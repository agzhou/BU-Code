%% Description

% GPU delay-and-sum on RF data collected from a row-column array probe (orthogonal plane wave compounding), using the
% rcabeam library (https://github.com/sdiebolt/rcabeam) through MATLAB's Python interface. Counterpart of RCA_DAS_ffdas.m
% with the same inputs, pixel grid (Verasonics PData by default) and delay model (plane wave TX time of flight referenced to the first-firing edge
% element, TW.peak, lens correction, startDepth), so it can replace it in the recon loop.

% Input: RF data and various parameters (same as RCA_DAS_ffdas), plus an optional Opts struct
% Output: IQ_CR (TX column, RX row) and IQ_RC (TX row, RX column) volumes, each coherently compounded over the angles.
%         Size (x, y, z, frames). With Opts.sumOnly, IQ_CR holds the coherent sum of both (OPW) and IQ_RC is empty
%         grid: struct with the x, y, z pixel coordinates in m
% Call RCA_DAS_rcabeam('clear') to release the GPU workspace that is kept between calls.

% Setup (once per MATLAB session):
%   - MATLAB R2026b+ with the Parallel Computing Toolbox, and a Python environment with cupy (CUDA 13) and rcabeam
%   - set the RCABEAM_PYTHON environment variable (or Opts.pythonExe) to that environment's python.exe
%   - run test_rcabeam_layout.m once (checks the MATLAB <-> Python array hand-over), then test_RCA_DAS_rcabeam.m

% Notes:
%           - rcabeam computes the receive delay from the in-plane distance sqrt(z^2 + (v - element)^2) and sums over
%             the angles inside the kernel, so no per-slice loop is needed (unlike RCA_DAS_ffdas). It does not do
%             per-angle output (no perAngle) or the Trans.ElementSens weighting that RCA_DAS uses. Opts.fNumber is a
%             hard receive aperture z/(2 fNumber)
%           - rcabeam shares one list of angles and first-sample times between RC and CR, and takes the same angles
%             for both. RCA_DAS flips the angles for the second group, which is the same list only when the angles are
%             symmetric about 0 and the array is square (the RC15gV setup). Otherwise this errors
%           - rcabeam needs IQ demodulated to absolute time, exp(-2j pi fc (t_start + n/fs)), so there is a constant
%             phase per angle on top of the usual demodulation. Without it the angles do not add coherently
%             (Opts.phaseReference = 'firstSample' turns it off to show the difference)
%           - Sample indices in rcabeam are 0-based and fractional (linear interpolation), unlike RCA_DAS' rounding
%           - Unchecked on real data: the Verasonics sample timing (see Opts.sampleShift) and the y index direction
%             relative to the Verasonics IQ. test_RCA_DAS_rcabeam.m checks the simulated case

% Opts fields (all optional):
%   frameChunk       (16)           # of frames per workspace call. The last chunk is zero padded to this size
%   sumOnly          (false)        one pass that returns the RC+CR sum in IQ_CR. false: two passes (RC, then CR)
%   fNumber          (0)            receive f-number, <= 0 for the full aperture
%   sampleShift      (0)            samples added to the sample index (e.g., -1 if the model is a sample off)
%   iqStorage        ('float32')    'float16' halves the packed input memory, lossy
%   phaseReference   ('absolute')   or 'firstSample' (see Notes)
%   backend          ('rcabeam')    'reference' runs rcabeam_reference_beamform on the CPU (for tests, slow)
%   pythonExe        ('')           python.exe for rcabeam_setup, if RCABEAM_PYTHON is not set
%   grid             ('verasonics') 'verasonics': the pixels of the Verasonics recon's PData (pixel 1 at PData.Origin, y
%                                   decreasing with the index). 'centered': the pixel centers of RCA_DAS's grid (y increasing)
%   verbose          (false)        print the adapter's timings for each chunk

function [IQ_CR, IQ_RC, grid] = RCA_DAS_rcabeam(RcvData, P, ypix_spacing, xpix_spacing, zpix_spacing, Opts)

    persistent cache % struct with the key of the configuration and the adapter session

    if ischar(RcvData) || isstring(RcvData)
        if strcmpi(RcvData, 'clear')
            if ~isempty(cache)
                try cache.sess.release(); catch, end
            end
            cache = [];
            IQ_CR = []; IQ_RC = []; grid = [];
            return
        end
        error('RCA_DAS_rcabeam:input', 'The only text input is "clear"')
    end

    if nargin < 6
        Opts = struct();
    end
    Opts = setDefault(Opts, 'frameChunk', 16);
    Opts = setDefault(Opts, 'sumOnly', false);
    Opts = setDefault(Opts, 'fNumber', 0);
    Opts = setDefault(Opts, 'sampleShift', 0);
    Opts = setDefault(Opts, 'iqStorage', 'float32');
    Opts = setDefault(Opts, 'phaseReference', 'absolute');
    Opts = setDefault(Opts, 'backend', 'rcabeam');
    Opts = setDefault(Opts, 'pythonExe', '');
    Opts = setDefault(Opts, 'grid', 'verasonics');
    Opts = setDefault(Opts, 'verbose', false);
    useRcabeam = strcmpi(Opts.backend, 'rcabeam');
    if ~useRcabeam && ~strcmpi(Opts.backend, 'reference')
        error('Unknown backend "%s"', Opts.backend)
    end
    useGPU = useRcabeam; % the reference runs on the host

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

    c = P.Resource.Parameters.speedOfSound;
    fc = Trans.frequency * 1e6;
    wl = c / fc;
    spw = Receive(1).samplesPerWave;
    fs = fc * spw;
    ns = Receive(1).endSample; % # samples per acquisition

    d_correction_1way = Trans.lensCorrection * wl; % TX and RX each pass through the lens

    %% Reconstruction region (same pixel grid as RCA_DAS and RCA_DAS_ffdas)
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
    elseif strcmpi(Opts.grid, 'centered')
        % Same pixel grid as RCA_DAS (pixel centers of a grid spanning the probe, y increasing)
        znumpix = ceil(zsize / zpix_spacing);
        xnumpix = ceil(xsize / xpix_spacing);
        ynumpix = ceil(ysize / ypix_spacing);

        grid.z = startDepth * wl + ((1:znumpix) - 1/2) / znumpix * zsize;
        grid.x = ((1:xnumpix) - (xnumpix + 1)/2) / xnumpix * xsize;
        grid.y = ((1:ynumpix) - (ynumpix + 1)/2) / ynumpix * ysize;
    else
        error('Unknown grid "%s". Use verasonics or centered.', Opts.grid)
    end

    elementx = (wl .* Trans.ElementPos(1:numElements, 1))';
    elementy = (wl .* Trans.ElementPos(numElements + 1 : numElements*2, 2))';
    x_end = max(elementx);
    y_end = max(elementy);

    %% rcabeam's angles and first-sample times, shared by both orientations
    % RCA_DAS: group 1 (TX columns steered in x-z, RX rows along y) uses angles; group 2 (TX rows steered in y-z, RX columns
    % along x) uses fliplr(angles) with the sign of the y term reversed, which is rcabeam's convention with the angle
    % -fliplr(angles). That equals angles when the list is symmetric about 0
    if max(abs(angles + fliplr(angles))) > 1e-6
        error('RCA_DAS_rcabeam:angles', 'The angles are not symmetric about 0, so RC and CR cannot share them in rcabeam')
    end
    if abs(x_end - y_end) > 1e-7
        error('RCA_DAS_rcabeam:array', 'x and y element extents differ (%.3g vs %.3g m); RC and CR cannot share t_start', x_end, y_end)
    end
    % rcabeam references the transmit to the origin, RCA_DAS to the first-firing edge element (+ x_end |sin|). The constants
    % of the round trip path go in with the (negative) first-sample time. First sample is at a round trip of 2*startDepth
    pathConst = wl * P.TW.peak + 2 * d_correction_1way - 2 * startDepth * wl;
    t0 = -(pathConst + x_end * abs(sin(angles))) / c - Opts.sampleShift / fs;

    %% The two acquisition groups of each frame (channels are positioned by RC15gV's connector map in r)
    grp1 = struct('acq', 1:na, 'ch', numElements+1 : 2*numElements); % CR: TX columns, RX rows (y)
    grp2 = struct('acq', na+1 : 2*na, 'ch', 1 : numElements);        % RC: TX rows, RX columns (x)
    G = struct('ns', ns, 'spw', spw, 'fc', fc, 't0', t0);

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
    nT = min(Opts.frameChunk, nf); % frames per workspace call

    %% Workspace (kept between calls with the same geometry)
    if useRcabeam
        mod = rcabeam_setup(Opts.pythonExe);
        key = string(sprintf('%d,', [ns, numElements, na, nT, xnumpix, ynumpix, znumpix])) + ...
              string(sprintf('%.12g,', [angles, t0, grid.x([1, end]), grid.y([1, end]), grid.z([1, end]), c, fs, fc, ...
                                        Opts.fNumber, elementx, elementy])) + string(Opts.iqStorage);
        if isempty(cache) || ~strcmp(cache.key, key)
            if ~isempty(cache)
                try cache.sess.release(); catch, end
                cache = [];
            end
            cache.sess = mod.Session(pyvec([ns, numElements, na, nT]), pyvec(angles), pyvec(t0), pyvec(grid.x), ...
                pyvec(grid.z), pyvec(grid.y), pyvec(elementx), pyvec(elementy), c, fs, fc, Opts.fNumber, char(Opts.iqStorage));
            cache.key = key;
        end
        sess = cache.sess;
    end

    %% Beamform in chunks of frames
    IQ_CR = complex(zeros(xnumpix, ynumpix, znumpix, nf, 'single'));
    if Opts.sumOnly
        IQ_RC = [];
    else
        IQ_RC = complex(zeros(xnumpix, ynumpix, znumpix, nf, 'single'));
    end
    want = [xnumpix, ynumpix, znumpix, nT];

    for f0 = 1 : nT : nf
        fi = f0 : min(f0 + nT - 1, nf);
        nfc = numel(fi);

        D_rc = prepGroup(r, fi, nT, Receive, grp2, G, Opts, useGPU);
        D_cr = prepGroup(r, fi, nT, Receive, grp1, G, Opts, useGPU);

        if useRcabeam
            pyRC = pycupyarray(D_rc);
            pyCR = pycupyarray(D_cr);
            clear D_rc D_cr
            if Opts.sumOnly
                volSum = rcabeam_from_py(sess.beamform_sum(pyRC, pyCR), want);
                IQ_CR(:, :, :, fi) = gather(volSum(:, :, :, 1:nfc));
                clear volSum
            else
                volCR = rcabeam_from_py(sess.beamform_cr(pyCR), want);
                IQ_CR(:, :, :, fi) = gather(volCR(:, :, :, 1:nfc));
                clear volCR
                volRC = rcabeam_from_py(sess.beamform_rc(pyRC), want);
                IQ_RC(:, :, :, fi) = gather(volRC(:, :, :, 1:nfc));
                clear volRC
            end
            if Opts.verbose
                fprintf('frames %d-%d: %s\n', fi(1), fi(end), char(sess.report()))
            end
        else
            args = {angles, t0, grid.x, grid.z, grid.y};
            volCR = rcabeam_reference_beamform(D_cr, 1, args{:}, elementy, c, fs, fc, Opts.fNumber);
            volRC = rcabeam_reference_beamform(D_rc, 0, args{:}, elementx, c, fs, fc, Opts.fNumber);
            if Opts.sumOnly
                IQ_CR(:, :, :, fi) = volCR(:, :, :, 1:nfc) + volRC(:, :, :, 1:nfc);
            else
                IQ_CR(:, :, :, fi) = volCR(:, :, :, 1:nfc);
                IQ_RC(:, :, :, fi) = volRC(:, :, :, 1:nfc);
            end
        end
    end

end


%% Channel data of one acquisition group as the IQ array the adapter takes: MATLAB [nT, na, nc, ns]
% (the same bytes as the workspace's C-contiguous (ns, nc, na, nT)). Frames beyond numel(fi) stay zero
function D = prepGroup(r, fi, nT, Receive, grp, G, Opts, useGPU)

    na = numel(grp.acq);
    nc = numel(grp.ch);
    ns = G.ns;
    nfc = numel(fi);
    if useGPU
        D = complex(zeros(nT, na, nc, ns, 'single', 'gpuArray'));
    else
        D = complex(zeros(nT, na, nc, ns, 'single'));
    end

    n = (0 : ns-1)';
    for k = 1:na
        a = grp.acq(k);
        x = r(Receive(a).startSample : Receive(a).endSample, grp.ch, fi); % (# samples, # channels, # frames), still int16 for the transfer
        if size(x, 1) ~= ns
            error('Acquisition %d has %d samples, expected %d', a, size(x, 1), ns)
        end
        if useGPU
            x = gpuArray(x);
        end
        x = analyticSignal(single(x));

        % Demodulate to absolute time: sample n is at t_start + n/fs, and the carrier phase must be referred to t = 0
        if strcmpi(Opts.phaseReference, 'absolute')
            ph = exp(-1j * 2*pi * (n / G.spw + G.fc * G.t0(k)));
        else
            ph = exp(-1j * 2*pi * n / G.spw);
        end
        x = x .* single(ph);

        D(1:nfc, k, :, :) = permute(x, [3, 4, 2, 1]);
    end

end


%% Analytic signal along the first dimension (same as hilbert, but runs on gpuArrays without a toolbox)
function y = analyticSignal(x)

    N = size(x, 1);
    h = zeros(N, 1, 'like', real(x));
    h(1) = 1;
    if mod(N, 2) == 0
        h(N/2 + 1) = 1;
        h(2 : N/2) = 2;
    else
        h(2 : (N+1)/2) = 2;
    end
    y = ifft(fft(x, [], 1) .* h, [], 1);

end


%% A numeric vector as a Python list of floats
function p = pyvec(v)

    p = py.list(num2cell(double(v(:)')));

end


function S = setDefault(S, name, val)

    if ~isfield(S, name)
        S.(name) = val;
    end

end
