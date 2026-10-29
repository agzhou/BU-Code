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
%             Opts.exactRx each voxel's coordinate along the element length is set to the channels' (0) for the
%             receive distance, which makes the 3D distance equal the in-plane distance. The voxel's real position
%             along that axis only enters through the transmit offsets, so each angle is one das call for the whole
%             volume (in depth chunks if it is large), not one call per slice.
%           - ffdas sums over all channels for every sequence entry, so the two orientations (different receive
%             channels) are separate das calls. They are not summed here; add them (or not) after.
%           - ffdas has no per-channel weights. Use Opts.rxHalfAngle for a hard directivity cone, not the
%             Trans.ElementSens weighting that RCA_DAS uses
%           - Checked with test_RCA_DAS_ffdas.m (simulated point scatterer): offsets index the samples from 0 (the
%             sampleShift default of 0 is right), and the data is demodulated with exp(-j*2*pi*n/samplesPerWave) so the
%             phase rotation is wavenum = +2*pi/samplesPerWave (wavenumSign = 1). The opposite sign loses most of the peak

% Opts fields (all optional):
%   exactRx      (true)   in-plane receive distance. false: 3D distance to the channel's midpoint (wrong off-center)
%   perAngle     (false)  keep the angle dimension instead of compounding
%   rxHalfAngle  (0.6)    half-angle in rad of the receive cone (ffdas srcdir), [] for no directivity mask. The default
%                         is the Verasonics Recon.senscutoff used in the Verasonics recon
%   grid         ('verasonics')  'verasonics': the pixels of the Verasonics recon's PData (pixel 1 at PData.Origin, y
%                         decreasing with the index). 'centered': the pixel centers of RCA_DAS's grid (y increasing)
%   angleWeights ([])     1 x na transmit apodization over the angles, in the order of P.angles
%   frameChunk   (Inf)    # of frames beamformed per batch, lower this if the GPU runs out of memory
%   maxElems     (6e7)    max # of elements of the (x, y, z, angles) offset/weight arrays per das call (splits in depth)
%   wavenumSign  (1)      sign of ffdas' phase rotation (see Notes)
%   sampleShift  (0)      samples added to the transmit offsets
%   useGPU       (true)   false keeps the arrays on the host, only for testing with a CPU stand-in for ffdas.das

function [IQ_CR, IQ_RC, grid] = RCA_DAS_ffdas(RcvData, P, ypix_spacing, xpix_spacing, zpix_spacing, Opts)

    if nargin < 6
        Opts = struct();
    end
    Opts = setDefault(Opts, 'exactRx', true);
    Opts = setDefault(Opts, 'perAngle', false);
    Opts = setDefault(Opts, 'rxHalfAngle', 0.6); % Verasonics' Recon.senscutoff
    Opts = setDefault(Opts, 'grid', 'verasonics');
    Opts = setDefault(Opts, 'angleWeights', []);
    Opts = setDefault(Opts, 'frameChunk', Inf);
    Opts = setDefault(Opts, 'wavenumSign', 1);
    Opts = setDefault(Opts, 'sampleShift', 0);
    Opts = setDefault(Opts, 'useGPU', true);
    Opts = setDefault(Opts, 'maxElems', 6e7);

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

    %% Reconstruction region
    xsize = numElements * Trans.spacingMm / 1e3;
    ysize = xsize;
    zsize = (endDepth - startDepth) * wl;

    if strcmpi(Opts.grid, 'verasonics')
        % Same pixels as the PData of the Verasonics recon (RCA_VerasonicsRecon_loop*.m): pixel (1, 1, 1) is at
        % PData.Origin = [-halfProbeDist, +halfProbeDist, startDepth], x increases and y DEcreases with the index
        if isfield(Trans, 'spacing')
            spacingWl = Trans.spacing;
        else
            spacingWl = Trans.spacingMm / 1e3 / wl;
        end
        L = Trans.spacingMm / 1e3 * numElements / wl; % probe width in wavelengths
        pdelta = [xpix_spacing, ypix_spacing, zpix_spacing] / wl; % [x, y, z] in wavelengths
        xnumpix = ceil(round(L / pdelta(1), 9));
        ynumpix = ceil(round(L / pdelta(2), 9));
        znumpix = ceil(round((endDepth - startDepth) / pdelta(3), 9));
        halfProbeDist = (numElements - 1) / 2 * spacingWl;

        grid.x = (-halfProbeDist + (0 : xnumpix-1) * pdelta(1)) * wl;
        grid.y = (halfProbeDist - (0 : ynumpix-1) * pdelta(2)) * wl;
        grid.z = (startDepth + (0 : znumpix-1) * pdelta(3)) * wl;
    else
        % Same pixel grid as RCA_DAS (pixel centers of a grid spanning the probe, y increasing)
        znumpix = ceil(zsize / zpix_spacing);
        xnumpix = ceil(xsize / xpix_spacing);
        ynumpix = ceil(ysize / ypix_spacing);

        grid.z = startDepth * wl + ((1:znumpix) - 1/2) / znumpix * zsize;
        grid.x = ((1:xnumpix) - (xnumpix + 1)/2) / xnumpix * xsize;
        grid.y = ((1:ynumpix) - (ynumpix + 1)/2) / ynumpix * ysize;
    end

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
    dDev = toDevice(d, Opts);
    clear d

    wavenum = single(Opts.wavenumSign * 2*pi / G.spw); % +2*pi*fc/fs for the exp(-j*2*pi*fc*t) demodulation above (verified, see the test)

    if isempty(Opts.rxHalfAngle)
        srcdir = [];
    else
        srcdir = toDevice(single(repmat([0; 0; 1; cos(Opts.rxHalfAngle)], 1, numCh)), Opts); % cone around +z
    end

    % Channel positions (3, # channels): along the element array axis only (zero along the slice axis and depth)
    src = toDevice(single([grp.uElem; zeros(2, numCh)] * G.ks), Opts);

    % Depths per call, limited so the (nu, ns, nzc) offsets/weights stay below Opts.maxElems elements
    zChunk = max(2, min(nz, floor(Opts.maxElems / (nu * ns))));

    if Opts.perAngle
        nk = na;
    else
        nk = 1;
    end
    vol = complex(zeros(nu, ns, nz, nk, nfc, 'single'));
    for z0 = 1 : zChunk : nz
        zi = z0 : min(z0 + zChunk - 1, nz);
        if numel(zi) < 2 % ffdas needs the voxel array to keep its depth dimension
            zi = max(1, zi(1) - 1) : zi(1);
        end
        nzc = numel(zi);

        % Voxel positions (3, nu, ns, nzc). With exactRx the slice coordinate is set to the channels' (0) for the
        % receive distance, so the 3D distance equals the in-plane distance, for the whole volume in one call. The
        % slice dependence of the transmit path is in the offsets below
        [U, S, Z] = ndgrid(grp.uv, grp.sv, G.zv(zi));
        if Opts.exactRx
            dst = permute(cat(4, U, zeros(size(S)), Z), [4, 1, 2, 3]);
        else
            dst = permute(cat(4, U, S, Z), [4, 1, 2, 3]);
        end
        dst = toDevice(single(dst * G.ks), Opts);
        wts = toDevice(ones(nu, ns, nzc, 'single'), Opts);

        acc = [];
        % One das call per angle: this ffdas build applies the first sequence entry's offsets/weights to every entry
        % of the sequence dimension (checked with impulses), so compounding all angles in one call is wrong
        for k = 1:na
            % TX one-way time of flight to each voxel, in samples (nu, ns, nzc). Same expression as RCA_DAS
            th = grp.angles(k);
            txPath = Z .* cos(th) + (grp.sg * S + sign(th) * grp.endPos) .* sin(th);
            offs = toDevice(single((txPath + G.pathConst) * G.ks + Opts.sampleShift), Opts);

            out = grp.w(k) * ffdas.das(dDev(:, k, :, :), src, dst, offs, wts, srcdir, wavenum);
            out = reshape(out, [nu, ns, nzc, 1, nfc]);
            if Opts.perAngle
                vol(:, :, zi, k, :) = toHost(out, Opts);
            elseif isempty(acc)
                acc = out;
            else
                acc = acc + out;
            end
        end
        if ~Opts.perAngle
            vol(:, :, zi, 1, :) = toHost(acc, Opts);
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
