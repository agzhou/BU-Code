%% Description

% GPU delay-and-sum on RF data collected from a row-column array probe (orthogonal plane wave compounding), using the
% ffdas library (https://github.com/BrainEchoLab/ffdas, MATLAB bindings). Counterpart of RCA_DAS.m: same inputs, same pixel
% grid, same delay model (plane wave TX time of flight referenced to the first-firing edge element, TW.peak, lens
% correction, startDepth).

% Input: RF data and various parameters (same as RCA_DAS), plus an optional Opts struct
% Output: IQ_CR (TX column, RX row) and IQ_RC (TX row, RX column) volumes, each coherently compounded over the angles.
%         Size (x, y, z, frames), or (x, y, z, angles, frames) if Opts.perAngle
%         grid: struct with the x, y, z pixel coordinates in m

% Notes:
%           - ffdas computes the receive delay from the 3D channel-to-voxel distance, but a row/column element is a
%             line: the correct delay only uses the distance in the plane perpendicular to the element. With
%             Opts.exactRx the volume is beamformed one slice at a time, and each channel's coordinate along its
%             length is set to the slice's coordinate, which makes the 3D distance equal the in-plane distance.
%           - ffdas sums over all channels for every sequence entry, so the two orientations (different receive
%             channels) are separate das calls. They are not summed here; add them (or not) after.
%           - ffdas has no per-channel weights. Use Opts.rxHalfAngle for a hard directivity cone, not the
%             Trans.ElementSens weighting that RCA_DAS uses
%           - Not yet checked against a real ffdas install: sample index base, the sign of the phase rotation (wavenum),
%             and out-of-range sample handling. test_RCA_DAS_ffdas.m checks them with a simulated point scatterer.

% Opts fields (all optional):
%   exactRx      (true)   per-slice in-plane receive distance. false: one call with 3D distance (faster, wrong off-center)
%   perAngle     (false)  keep the angle dimension instead of compounding
%   rxHalfAngle  ([])     half-angle in rad of the receive cone (ffdas srcdir), [] for no directivity mask
%   angleWeights ([])     1 x na transmit apodization over the angles, in the order of P.angles
%   frameChunk   (Inf)    # of frames beamformed per batch, lower this if the GPU runs out of memory
%   wavenumSign  (1)      set to -1 if ffdas' phase rotation convention is opposite
%   sampleShift  (0)      samples added to the transmit offsets (e.g., -1 if ffdas indexes samples from 1)
%   useGPU       (true)   false keeps the arrays on the host, only for testing with a CPU stand-in for ffdas.das

function [IQ_CR, IQ_RC, grid] = RCA_DAS_ffdas(RcvData, P, ypix_spacing, xpix_spacing, zpix_spacing, Opts)

    if nargin < 6
        Opts = struct();
    end
    Opts = setDefault(Opts, 'exactRx', true);
    Opts = setDefault(Opts, 'perAngle', false);
    Opts = setDefault(Opts, 'rxHalfAngle', []);
    Opts = setDefault(Opts, 'angleWeights', []);
    Opts = setDefault(Opts, 'frameChunk', Inf);
    Opts = setDefault(Opts, 'wavenumSign', 1);
    Opts = setDefault(Opts, 'sampleShift', 0);
    Opts = setDefault(Opts, 'useGPU', true);

    %% Parameters
    Trans = P.Trans;
    Receive = P.Receive;
    na = P.na;
    numElements = P.numElements;
    angles = reshape(P.angles, 1, []);
    startDepth = P.startDepth;
    endDepth = P.endDepth;
    if numel(angles) ~= na
        error('P.angles has %d entries but P.na = %d', numel(angles), na)
    end
    if isempty(Opts.angleWeights)
        angleW = ones(1, na);
    else
        angleW = reshape(Opts.angleWeights, 1, []);
    end

    c = P.Resource.Parameters.speedOfSound;
    fc = Trans.frequency * 1e6;
    wl = c / fc;
    spw = Receive(1).samplesPerWave;
    fs = fc * spw;
    ks = fs / c; % m -> sampling wavelengths (c / fs), the unit ffdas takes for positions and offsets

    d_correction_1way = Trans.lensCorrection * wl; % TX and RX each pass through the lens

    %% Reconstruction region (same pixel grid as RCA_DAS)
    xsize = numElements * Trans.spacingMm / 1e3;
    ysize = xsize;
    zsize = (endDepth - startDepth) * wl;

    znumpix = ceil(zsize / zpix_spacing);
    xnumpix = ceil(xsize / xpix_spacing);
    ynumpix = ceil(ysize / ypix_spacing);

    grid.z = startDepth * wl + ((1:znumpix) - 1/2) / znumpix * zsize;
    grid.x = ((1:xnumpix) - (xnumpix + 1)/2) / xnumpix * xsize;
    grid.y = ((1:ynumpix) - (ynumpix + 1)/2) / ynumpix * ysize;

    elementx = (wl .* Trans.ElementPos(1:numElements, 1))';
    elementy = (wl .* Trans.ElementPos(numElements + 1 : numElements*2, 2))';
    x_end = max(elementx);
    y_end = max(elementy);

    % Terms of the round trip path that do not depend on the pixel or the angle, in m. The padding RCA_DAS adds for a
    % nonzero start depth is replaced by this shift of the sample axis (first sample is at a round trip of 2*startDepth)
    pathConst = wl * P.TW.peak + 2 * d_correction_1way - 2 * startDepth * wl;

    G = struct('zv', grid.z, 'ks', ks, 'spw', spw, 'pathConst', pathConst);

    %% The two acquisition groups of each frame
    % Group 1 (acqs 1:na): TX columns steered in x-z, RX on the rows (channels numElements+1:2*numElements, positioned along y)
    % Group 2 (acqs na+1:2*na): TX rows steered in y-z, RX on the columns (channels 1:numElements, positioned along x).
    % RCA_DAS flips the angles for this group
    %   u: axis the receive elements are positioned along (resolved by receive beamforming)
    %   s: other lateral axis (resolved by the transmit steering). Slice axis for the in-plane receive distance
    %   sg: sign of s in the TX path, endPos: edge element position of the transmit aperture along s
    grp1 = struct('acq', 1:na, 'ch', numElements+1 : 2*numElements, 'uElem', elementy, 'uv', grid.y, 'sv', grid.x, ...
                  'sg', 1, 'endPos', x_end, 'angles', angles, 'w', angleW);
    grp2 = struct('acq', na+1 : 2*na, 'ch', 1 : numElements, 'uElem', elementx, 'uv', grid.x, 'sv', grid.y, ...
                  'sg', -1, 'endPos', y_end, 'angles', fliplr(angles), 'w', fliplr(angleW));

    %% Get raw data and remap according to element-channel pairings
    if iscell(RcvData)
        r = RcvData{1, 1};
    else
        r = RcvData;
    end
    r = r(:, Trans.Connector, :); % (# z samples * # acquisitions, # channels, # frames)
    if isfield(P, 'numSubFrames') % older P structures name the frame count differently, so take it from the data
        nf = P.numSubFrames;
    else
        nf = size(r, 3);
    end

    if Opts.perAngle
        nk = na;
    else
        nk = 1;
    end
    IQ_CR = complex(zeros(xnumpix, ynumpix, znumpix, nk, nf, 'single'));
    IQ_RC = complex(zeros(xnumpix, ynumpix, znumpix, nk, nf, 'single'));

    chunkSize = min(Opts.frameChunk, nf);
    for f0 = 1 : chunkSize : nf
        fi = f0 : min(f0 + chunkSize - 1, nf);

        % grp1: u = y, s = x -> (y, x, z, k, f)
        vol = groupDAS(r(:, :, fi), Receive, grp1, G, Opts);
        IQ_CR(:, :, :, :, fi) = permute(vol, [2, 1, 3, 4, 5]);

        % grp2: u = x, s = y -> (x, y, z, k, f)
        IQ_RC(:, :, :, :, fi) = groupDAS(r(:, :, fi), Receive, grp2, G, Opts);
    end

    if ~Opts.perAngle
        IQ_CR = reshape(IQ_CR, xnumpix, ynumpix, znumpix, nf);
        IQ_RC = reshape(IQ_RC, xnumpix, ynumpix, znumpix, nf);
    end

end


%% Beamform one acquisition group. Output is (u, s, z, k, frames)
function vol = groupDAS(rf, Receive, grp, G, Opts)

    na = numel(grp.angles);
    nu = numel(grp.uv);
    ns = numel(grp.sv);
    nz = numel(G.zv);
    numCh = numel(grp.uElem);
    nfc = size(rf, 3);
    nS = Receive(1).endSample; % # samples per acquisition

    % (# z samples, # angles, # channels, # frames)
    d = zeros(nS, na, numCh, nfc, 'single');
    for k = 1:na
        a = grp.acq(k);
        d(:, k, :, :) = permute(rf(Receive(a).startSample:Receive(a).endSample, grp.ch, :), [1, 4, 2, 3]);
    end

    % Analytic signal, then demodulate to baseband (ffdas' interpolation works on baseband IQ with a phase rotation)
    d = reshape(hilbert(reshape(d, nS, [])), size(d));
    d = d .* exp(-1j * 2*pi * single((0 : nS-1)') / G.spw);
    d = toDevice(d, Opts);

    wavenum = single(-Opts.wavenumSign * 2*pi / G.spw); % -2*pi*fc/fs

    if isempty(Opts.rxHalfAngle)
        srcdir = [];
    else
        srcdir = toDevice(single(repmat([0; 0; 1; cos(Opts.rxHalfAngle)], 1, numCh)), Opts); % cone around +z
    end

    if Opts.exactRx
        sliceSets = num2cell(1:ns); % one slice per call, channels put in the slice's plane
    else
        sliceSets = {1:ns};
    end
    if Opts.perAngle
        angSets = num2cell(1:na);
    else
        angSets = {1:na};
    end
    dSets = cell(size(angSets));
    for k = 1:numel(angSets)
        dSets{k} = d(:, angSets{k}, :, :);
    end

    th = reshape(grp.angles, 1, 1, 1, na);
    angleW = reshape(grp.w, 1, 1, 1, na);

    vol = complex(zeros(nu, ns, nz, numel(angSets), nfc, 'single'));

    for si = 1:numel(sliceSets)
        sidx = sliceSets{si};
        nsi = numel(sidx);
        sv = grp.sv(sidx);

        % Channel positions (3, # channels): along the element array axis, along the slice axis, depth
        if Opts.exactRx
            src = [grp.uElem; sv * ones(1, numCh); zeros(1, numCh)];
        else
            src = [grp.uElem; zeros(2, numCh)];
        end
        src = toDevice(single(src * G.ks), Opts);

        % Voxel positions (3, nu, nsi, nz)
        [U, S, Z] = ndgrid(grp.uv, sv, G.zv);
        dst = permute(cat(4, U, S, Z), [4, 1, 2, 3]);

        % TX one-way time of flight to each voxel, in samples (nu, nsi, nz, na). Same expression as RCA_DAS
        txPath = Z .* cos(th) + (grp.sg * S + sign(th) * grp.endPos) .* sin(th);
        offs = single((txPath + G.pathConst) * G.ks + Opts.sampleShift);
        wts = single(ones(size(offs)) .* angleW);

        if nsi == 1 % drop the singleton slice dim
            dst = reshape(dst, [3, nu, nz]);
        end
        dst = toDevice(single(dst * G.ks), Opts);

        for k = 1:numel(angSets)
            a = angSets{k};
            o = offs(:, :, :, a);
            w = wts(:, :, :, a);
            if nsi == 1
                o = reshape(o, [nu, nz, numel(a)]);
                w = reshape(w, [nu, nz, numel(a)]);
            end
            out = ffdas.das(dSets{k}, src, dst, toDevice(o, Opts), toDevice(w, Opts), srcdir, wavenum);
            vol(:, sidx, :, k, :) = reshape(toHost(out, Opts), [nu, nsi, nz, 1, nfc]);
        end
    end

end


function x = toDevice(x, Opts)
    if Opts.useGPU
        x = gpuArray(x);
    end
end


function x = toHost(x, Opts)
    if Opts.useGPU
        x = gather(x);
    end
end


function S = setDefault(S, name, val)
    if ~isfield(S, name)
        S.(name) = val;
    end
end
