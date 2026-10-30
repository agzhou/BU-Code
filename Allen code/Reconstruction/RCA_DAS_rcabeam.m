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
%           - Receive.sampleMode (all Receives must have the same one): NS200BW, BS67BW, BS100BW and BS50BW, the
%             Verasonics sampling schemes of Vantage Sequence Programming Manual 3.3.1 (NS200BWI and 'custom' are not
%             supported). The BS modes store fewer samples (BS67BW: 3/4 of a period apart, BS100BW and BS50BW: a pair a
%             quarter period apart every 1 or 2 periods) of a band that is narrow around Receive.demodFrequency Fd, so the
%             band is aliased and folded. Each acquisition's stored samples are turned into the complex baseband signal at
%             Fd and upsampled to 4*Fd (the rate of NS200BW) in the frequency domain, which is exact for a band limited
%             signal, and that is what the kernel interpolates (fs = 4*Fd, f_demod = Fd). Same conversion as RCA_DAS_ffdas.
%             The first stored sample is at t = 0 for all modes. The demodulation frequency is Receive.demodFrequency, not
%             Trans.frequency (13.39 vs. 13.6 MHz for the BS67BW RC15gV fUS runs)
%           - The IQ is returned as the complex CONJUGATE of the kernel's output (Opts.conjugateIQ), to follow the Verasonics
%             convention. Found on the in vivo BS67BW RC15gV data (test_rcabeam_vs_verasonics_fus.m): with the kernel's own sign
%             the clutter filtered IQ correlates 0.04 per voxel with the Verasonics IQ and the color Doppler has the opposite
%             sign (correlation -0.99), after conjugation 0.92 and +0.99. Magnitudes (PDI, g1) are not affected. Only
%             checked on BS67BW; the pair modes (BS100BW, BS50BW) assume the same convention
%           - rcabeam needs IQ demodulated to absolute time, exp(-2j pi Fd (t_start + n/fs)), so there is a constant
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
%   conjugateIQ      (true)         output the complex conjugate of what the rcabeam kernel computes, which is the Verasonics IQ
%                                   convention (see Notes). false: the kernel's own (a positive Doppler shift is a flow toward
%                                   the probe), which has the opposite sign of the Verasonics IQ's phase
%   yIndex           ([])           reconstruct only these rows of the y grid (indices into grid.y), e.g., one xz plane. The
%                                   output and grid.y have only these rows. [] for all
%   sampleMode       ([])           override P.Receive(1).sampleMode: NS200BW, BS67BW, BS100BW or BS50BW (see Notes)
%   demodFrequency   ([])           override P.Receive(1).demodFrequency, in MHz
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
    Opts = setDefault(Opts, 'conjugateIQ', true);
    Opts = setDefault(Opts, 'yIndex', []);
    Opts = setDefault(Opts, 'sampleMode', []);
    Opts = setDefault(Opts, 'demodFrequency', []);
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
    ns = Receive(1).endSample; % # stored samples per acquisition (fewer than the Nyquist count in the BS sample modes)

    % The kernel interpolates complex baseband IQ at one rate fs, so the stored samples of any sampleMode are first turned
    % into baseband at demodFrequency Fd, at 4*Fd (see basebandSetup/toBaseband), nsIQ samples per acquisition
    bb = basebandSetup(Receive(1), Opts, fc, ns);
    fs = bb.fsEff;
    Fd = bb.Fd;
    nsIQ = bb.nIQ;

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
    if ~isempty(Opts.yIndex) % only some rows of the grid, e.g., one xz plane
        grid.y = grid.y(Opts.yIndex);
        ynumpix = numel(grid.y);
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
    G = struct('ns', ns, 'bb', bb, 't0', t0);

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
        key = string(sprintf('%d,', [nsIQ, numElements, na, nT, xnumpix, ynumpix, znumpix])) + ...
              string(sprintf('%.12g,', [angles, t0, grid.x([1, end]), grid.y([1, end]), grid.z([1, end]), c, fs, Fd, ...
                                        Opts.fNumber, elementx, elementy])) + string(Opts.iqStorage);
        if isempty(cache) || ~strcmp(cache.key, key)
            if ~isempty(cache)
                try cache.sess.release(); catch, end
                cache = [];
            end
            cache.sess = mod.Session(pyvec([nsIQ, numElements, na, nT]), pyvec(angles), pyvec(t0), pyvec(grid.x), ...
                pyvec(grid.z), pyvec(grid.y), pyvec(elementx), pyvec(elementy), c, fs, Fd, Opts.fNumber, char(Opts.iqStorage));
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
    if Opts.conjugateIQ
        sgn = @(v) conj(v); % the Verasonics IQ is the complex conjugate of the rcabeam kernel's (see Notes)
    else
        sgn = @(v) v;
    end

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
                IQ_CR(:, :, :, fi) = gather(sgn(volSum(:, :, :, 1:nfc)));
                clear volSum
            else
                volCR = rcabeam_from_py(sess.beamform_cr(pyCR), want);
                IQ_CR(:, :, :, fi) = gather(sgn(volCR(:, :, :, 1:nfc)));
                clear volCR
                volRC = rcabeam_from_py(sess.beamform_rc(pyRC), want);
                IQ_RC(:, :, :, fi) = gather(sgn(volRC(:, :, :, 1:nfc)));
                clear volRC
            end
            if Opts.verbose
                fprintf('frames %d-%d: %s\n', fi(1), fi(end), char(sess.report()))
            end
        else
            args = {angles, t0, grid.x, grid.z, grid.y};
            volCR = rcabeam_reference_beamform(D_cr, 1, args{:}, elementy, c, fs, Fd, Opts.fNumber);
            volRC = rcabeam_reference_beamform(D_rc, 0, args{:}, elementx, c, fs, Fd, Opts.fNumber);
            if Opts.sumOnly
                IQ_CR(:, :, :, fi) = sgn(volCR(:, :, :, 1:nfc) + volRC(:, :, :, 1:nfc));
            else
                IQ_CR(:, :, :, fi) = sgn(volCR(:, :, :, 1:nfc));
                IQ_RC(:, :, :, fi) = sgn(volRC(:, :, :, 1:nfc));
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
        D = complex(zeros(nT, na, nc, G.bb.nIQ, 'single', 'gpuArray'));
    else
        D = complex(zeros(nT, na, nc, G.bb.nIQ, 'single'));
    end

    for k = 1:na
        a = grp.acq(k);
        x = r(Receive(a).startSample : Receive(a).endSample, grp.ch, fi); % (# samples, # channels, # frames), still int16 for the transfer
        if size(x, 1) ~= ns
            error('Acquisition %d has %d samples, expected %d', a, size(x, 1), ns)
        end
        if useGPU
            x = gpuArray(x);
        end

        % Stored samples (any sampleMode) -> complex baseband at Fd, sampled at 4*Fd: (G.bb.nIQ, # channels, # frames)
        x = toBaseband(single(x), G.bb, useGPU);

        % Demodulate to absolute time: sample n is at t_start + n/fs, and the carrier phase must be referred to t = 0.
        % toBaseband took out exp(-j 2 pi Fd n/fs), what is left is the constant phase of the first sample's time
        if strcmpi(Opts.phaseReference, 'absolute')
            x = x * single(exp(-1j * 2*pi * G.bb.Fd * G.t0(k)));
        end

        D(1:nfc, k, :, :) = permute(x, [3, 4, 2, 1]);
    end

end


%% How the stored samples of a Receive relate to the demodulation frequency Fd (Vantage Sequence Programming Manual, 3.3.1)
%   NS200BW: uniform, 4 samples per period of Fd (fsN = 4*Fd)
%   BS67BW:  uniform, 1 sample every 3/4 period (fsN = 4/3*Fd). The band of Fd is aliased and folded
%   BS100BW: pairs of samples, at t = m*Tp and m*Tp + 1/(4*Fd), with Tp = 1/Fd (the pair is a quadrature I, Q pair)
%   BS50BW:  same pairs with Tp = 2/Fd
% The complex baseband (at carrier Fd) is rebuilt from them and upsampled to fsEff = 4*Fd, the rate of NS200BW, which is
% fine enough for the kernel's linear interpolation. Same as RCA_DAS_ffdas (validated there, BS67BW on in vivo data)
function bb = basebandSetup(R, Opts, fc, nS)

    if ~isempty(Opts.sampleMode)
        mode = Opts.sampleMode;
    elseif isfield(R, 'sampleMode') && ~isempty(R.sampleMode)
        mode = R.sampleMode;
    else
        mode = 'NS200BW';
    end
    mode = upper(char(mode));

    if ~isempty(Opts.demodFrequency)
        Fd = Opts.demodFrequency * 1e6;
    elseif isfield(R, 'demodFrequency') && ~isempty(R.demodFrequency)
        Fd = R.demodFrequency * 1e6;
    else % older P structures: NS200BW, demodFrequency = decimSampleRate / 4
        Fd = fc * R.samplesPerWave / 4;
    end

    bb.mode = mode;
    bb.Fd = Fd;
    bb.fsEff = 4 * Fd;
    switch mode
        case 'NS200BW'
            bb.fsN = 4 * Fd; bb.pairs = false; bb.L = 1;
        case 'BS67BW'
            bb.fsN = 4/3 * Fd; bb.pairs = false; bb.L = 3;
        case 'BS100BW'
            bb.fsN = Fd; bb.pairs = true; bb.L = 4; % fsN is the rate of the baseband samples (1 pair per period)
        case 'BS50BW'
            bb.fsN = Fd / 2; bb.pairs = true; bb.L = 8;
        otherwise
            error('Receive.sampleMode %s is not supported. Use NS200BW, BS67BW, BS100BW or BS50BW.', mode)
    end
    if bb.pairs
        if mod(nS, 2) ~= 0
            error('Receive has %d samples, but %s samples come in pairs', nS, mode)
        end
        bb.nIQ = bb.L * nS / 2; % # of baseband samples after the upsampling to fsEff
    else
        bb.nIQ = bb.L * nS;
    end

    % The stored sample rate should be decimSampleRate / quadDecim (BS100BW and BS50BW are not checked)
    if ~bb.pairs && isfield(R, 'decimSampleRate') && isfield(R, 'quadDecim') && ...
            abs(R.decimSampleRate * 1e6 / R.quadDecim / bb.fsN - 1) > 0.01
        warning('RCA_DAS_rcabeam:sampleRate', ['Stored sample rate decimSampleRate/quadDecim = %.4g MHz does not match ' ...
                '%s with demodFrequency = %.4g MHz (%.4g MHz expected)'], R.decimSampleRate / R.quadDecim, mode, Fd / 1e6, bb.fsN / 1e6)
    end

end


%% Stored RF samples (nS, ...) -> complex baseband at carrier bb.Fd, bb.nIQ rows, at bb.fsEff
% The first stored sample is at t = 0. Exact for a band limited signal (frequency domain filtering and interpolation)
function a = toBaseband(d, bb, useGPU)

    sz = size(d);
    nS = sz(1);
    d = reshape(d, nS, []);

    if bb.pairs
        % x(t) = Re(a(t) exp(j 2 pi Fd t)) = I(t) cos(.) - Q(t) sin(.). At t = m*Tp x = I and at t = m*Tp + 1/(4 Fd) x = -Q.
        % I and Q are real baseband signals limited to |f| < 1/(2*Tp), sampled at 1/Tp, so they can be rebuilt exactly.
        % Q is sampled a quarter period late: Q(f) = Qs(f) exp(-j 2 pi f tau)
        M = nS / 2;
        tauSamples = (bb.fsN / (4 * bb.Fd)); % quarter period of Fd in baseband samples
        k = [0 : ceil(M/2) - 1, -floor(M/2) : -1]';
        shiftQ = toDevice(single(exp(-2j * pi * (k / M) * tauSamples)), useGPU);
        Q = real(ifft(fft(-d(2:2:end, :), [], 1) .* shiftQ, [], 1));
        a = complex(d(1:2:end, :), Q);
        Y = fft(a, [], 1);
        % zero padding of the spectrum, the Nyquist bin is split between the two sides
        h = M / 2;
        Yu = zeros(bb.L * M, size(Y, 2), 'like', Y);
        Yu(1:h, :) = Y(1:h, :);
        Yu(end-h+1 : end, :) = Y(h+1 : end, :);
        Yu(h+1, :) = 0.5 * Y(h+1, :);
        Yu(end-h+1, :) = 0.5 * Y(h+1, :);
    else
        % Demodulate at the sample times, then keep |f| < fsN/4 (the image of the other half of the spectrum is in
        % fsN/4 < |f| < 3*fsN/4 for both modes). For NS200BW this is the analytic signal
        M = nS;
        ph = (0 : nS-1)' * (bb.Fd / bb.fsN);
        ph = ph - floor(ph); % cycles, wrapped in double precision before the conversion
        y = 2 * d .* toDevice(single(exp(-2j * pi * ph)), useGPU);
        Y = fft(y, [], 1);
        K = floor(M / 4);
        w = ones(2*K + 1, 1);
        if mod(M, 4) == 0
            w([1, end]) = 0.5; % bins at exactly +-fsN/4 are shared with the image
        end
        Yu = zeros(bb.L * M, size(Y, 2), 'like', Y);
        Yu(1 : K+1, :) = w(K+1 : end) .* Y(1 : K+1, :);
        Yu(end-K+1 : end, :) = w(1 : K) .* Y(end-K+1 : end, :);
    end

    a = bb.L * ifft(Yu, [], 1);
    a = reshape(a, [size(a, 1), sz(2:end)]);

end


function x = toDevice(x, useGPU)

    if useGPU
        x = gpuArray(x);
    end

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
