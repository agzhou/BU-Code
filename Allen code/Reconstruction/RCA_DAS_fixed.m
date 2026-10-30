%% Description

% Copy of RCA_DAS.m (delay-and-sum on RF data from a row-column array probe, orthogonal plane wave compounding) with the
% phase and zero-sample problems fixed, and with the Verasonics sample modes handled. The delay model (TX/RX path lengths,
% edge-referenced TX, TW.peak, lens correction, ElementSens angles) is RCA_DAS's. Each change is a switch in Opts so that
% its effect can be tested (test_RCA_DAS_vs_rcabeam.m): the defaults are the fixed behavior, and
%   Opts = struct('iq', 'analytic', 'zeroSample', false, 'apod', 'sens', 'grid', 'centered', 'conjugateIQ', false)
% is RCA_DAS's behavior (apart from the sampleMode handling, which RCA_DAS does not have).

% What RCA_DAS did, and what is different here:
%   1. Zero sample. RCA_DAS rounds the 0-based sample position round(path / (c Ts)) and uses it as a 1-based MATLAB index,
%      so every channel is read one sample early (a 1/(4 Fd) = 19 ns error, a quarter period: ~90 degrees of carrier phase on
%      every channel, and a lower coherent gain). Opts.zeroSample (true) adds the 1.
%      The first-sample time also was startDepth*wl one way, the first sample of a Verasonics acquisition is at a round trip
%      of 2*startDepth*wl (the same for startDepth = 0)
%   2. Phase. RCA_DAS takes the analytic signal of the stored samples (hilbert) and reads the nearest sample, so the delay
%      is off by up to half a sample (+-1/(8 Fd), +-45 degrees of the carrier) in a different way for every channel, which
%      does not add up coherently. Opts.iq = 'baseband' (default) converts the data to complex baseband, reads it at the
%      fractional sample position with linear interpolation, and multiplies by exp(+j 2 pi Fd (n/fs)) there, so the carrier
%      phase is right for the exact delay (what the rcabeam kernel does). 'analytic' is RCA_DAS's nearest sample
%   3. Sample modes. hilbert() is only valid for data sampled at 4 samples per period (NS200BW). The BS67BW, BS100BW and BS50BW
%      samples (Vantage Sequence Programming Manual 3.3.1) are aliased and folded: they are converted to the complex baseband at
%      Receive.demodFrequency Fd and upsampled to fs = 4 Fd in the frequency domain (same conversion as RCA_DAS_ffdas and
%      RCA_DAS_rcabeam)
%   4. Apodization. RCA_DAS weights every channel by Trans.ElementSens(angle) ('sens'). The Verasonics recon includes a channel
%      or not (Recon.senscutoff, 0.6), which is a hard f-number (0.788 for the RC15gV): Opts.apod = 'cutoff' (default). 'sens+cutoff'
%      weights the included channels
%   5. Output sign. The Verasonics IQ is the complex conjugate of this phase convention (see RCA_DAS_rcabeam), so the output is
%      conjugated (Opts.conjugateIQ, default true)
%   6. Pixel grid: Opts.grid 'verasonics' (default, the pixels of the Verasonics PData: pixel 1 at [-halfProbeDist, +halfProbeDist,
%      startDepth], y decreasing) or RCA_DAS's 'centered'. Opts.yIndex limits the rows of the y grid
%   7. Speed only: no parfor and no per-pixel loops over x: for each (y, z) all x pixels, channels and acquisitions are
%      gathered at once, on the GPU if there is one, a chunk of frames at a time. Same sums as RCA_DAS's per pixel code

% Input: RF data (unstacked: (# samples * # acquisitions, # channels, # frames), or a cell with that) and various
%        parameters (as RCA_DAS)
% Output: IQ (x, y, z, # acquisitions (2 * # angles), # frames), or (x, y, z, 1, # frames) with Opts.sumAcq.
%         The acquisitions: 1:na TX columns / RX rows, na+1:2na TX rows / RX columns (RCA_DAS's order)
%         grid: pixel coordinates in m
% Opts fields (all optional):
%   iq             ('baseband')   'baseband' or 'analytic' (see 2.)
%   zeroSample     (true)         see 1.
%   apod           ('cutoff')     'cutoff', 'sens' or 'sens+cutoff' (see 4.)
%   senscutoff     (0.6)          Recon.senscutoff
%   grid           ('verasonics') or 'centered'
%   yIndex         ([])           rows of the y grid to reconstruct, [] for all
%   conjugateIQ    (true)
%   sumAcq         (false)        return the coherent sum over the acquisitions (OPW)
%   frameChunk     (64)           frames per batch
%   useGPU         (true if there is a GPU)
%   sampleMode, demodFrequency    override P.Receive(1).sampleMode, .demodFrequency (MHz)

function [IQ, grid] = RCA_DAS_fixed(RcvData, P, ypix_spacing, xpix_spacing, zpix_spacing, Opts)

    if nargin < 6
        Opts = struct();
    end
    Opts = setDefault(Opts, 'iq', 'baseband');
    Opts = setDefault(Opts, 'zeroSample', true);
    Opts = setDefault(Opts, 'apod', 'cutoff');
    Opts = setDefault(Opts, 'senscutoff', 0.6);
    Opts = setDefault(Opts, 'grid', 'verasonics');
    Opts = setDefault(Opts, 'yIndex', []);
    Opts = setDefault(Opts, 'conjugateIQ', true);
    Opts = setDefault(Opts, 'sumAcq', false);
    Opts = setDefault(Opts, 'frameChunk', 64);
    Opts = setDefault(Opts, 'useGPU', gpuDeviceCount > 0);
    Opts = setDefault(Opts, 'sampleMode', []);
    Opts = setDefault(Opts, 'demodFrequency', []);
    useGPU = Opts.useGPU;
    interp = strcmpi(Opts.iq, 'baseband');
    if ~interp && ~strcmpi(Opts.iq, 'analytic')
        error('Opts.iq must be baseband or analytic')
    end

    %% Parameters
    Trans = P.Trans;
    Receive = P.Receive;
    na = P.na;
    numElements = P.numElements;
    angles = reshape(P.angles, 1, []);
    startDepth = P.startDepth;
    endDepth = P.endDepth;
    c = P.Resource.Parameters.speedOfSound;
    fc = Trans.frequency * 1e6;
    wl = c / fc;
    tw_peak = P.TW.peak;
    d_correction_1way = Trans.lensCorrection * wl; % TX and RX each pass through the lens
    ns = Receive(1).endSample; % # stored samples per acquisition

    bb = basebandSetup(Receive(1), Opts, fc, ns);
    fs = bb.fsEff; Fd = bb.Fd; nIQ = bb.nIQ;
    cTs = c / fs; % distance per sample

    %% Pixel grid
    xsize = numElements * Trans.spacingMm / 1e3;
    ysize = xsize;
    zsize = (endDepth - startDepth) * wl;
    if strcmpi(Opts.grid, 'verasonics')
        spacingWl = Trans.spacing;
        L = Trans.spacingMm / 1e3 * numElements / wl;
        pdelta = [xpix_spacing, ypix_spacing, zpix_spacing] / wl;
        xnumpix = ceil(round(L / pdelta(1), 9));
        ynumpix = ceil(round(L / pdelta(2), 9));
        znumpix = ceil(round((endDepth - startDepth) / pdelta(3), 9));
        halfProbeDist = (numElements - 1) / 2 * spacingWl;
        grid.x = (-halfProbeDist + (0 : xnumpix-1) * pdelta(1)) * wl;
        grid.y = (halfProbeDist - (0 : ynumpix-1) * pdelta(2)) * wl;
        grid.z = (startDepth + (0 : znumpix-1) * pdelta(3)) * wl;
    elseif strcmpi(Opts.grid, 'centered')
        znumpix = ceil(zsize / zpix_spacing);
        xnumpix = ceil(xsize / xpix_spacing);
        ynumpix = ceil(ysize / ypix_spacing);
        grid.z = startDepth * wl + ((1:znumpix) - 1/2) / znumpix * zsize;
        grid.x = ((1:xnumpix) - (xnumpix + 1)/2) / xnumpix * xsize;
        grid.y = ((1:ynumpix) - (ynumpix + 1)/2) / ynumpix * ysize;
    else
        error('Unknown grid "%s"', Opts.grid)
    end
    if ~isempty(Opts.yIndex)
        grid.y = grid.y(Opts.yIndex);
        ynumpix = numel(grid.y);
    end

    elementx = (wl .* Trans.ElementPos(1:numElements, 1))';
    elementy = (wl .* Trans.ElementPos(numElements + 1 : numElements*2, 2))';
    x_end = max(elementx);
    y_end = max(elementy);

    %% Element sensitivity (RCA_DAS: ElementSens interpolated at the angle from the pixel to the element)
    esAngle = linspace(-pi/2, pi/2, numel(Trans.ElementSens));
    esAng = esAngle(esAngle >= -1e-9 & Trans.ElementSens(:)' < Opts.senscutoff);
    if isempty(esAng)
        thetaCut = pi/2;
    else
        thetaCut = esAng(1); % the aperture ends at the first tabulated angle that is below the cutoff
    end
    apodFun = @(ang) apodization(ang, Opts.apod, Trans.ElementSens, esAngle, thetaCut);

    %% Raw data and remap according to element-channel pairings
    if iscell(RcvData)
        r = RcvData{1, 1};
    else
        r = RcvData;
    end
    r = r(:, Trans.Connector, :);
    if isfield(P, 'numSubFrames')
        nf = P.numSubFrames;
    else
        nf = size(r, 3);
    end
    nAcq = 2 * na;
    grpCh = {numElements + 1 : 2*numElements, 1 : numElements}; % channels of the acquisitions 1:na and na+1:2na

    angCols = reshape(angles, 1, na); % angles of group 1 (acq a), group 2 uses fliplr(angles) (RCA_DAS)
    angG2 = fliplr(angCols);

    ex = reshape(elementx, [], 1); ey = reshape(elementy, [], 1);
    xg = reshape(grid.x, 1, 1, []); % x pixels along dim 3

    if Opts.sumAcq
        IQ = complex(zeros(xnumpix, ynumpix, znumpix, 1, nf, 'single'));
    else
        IQ = complex(zeros(xnumpix, ynumpix, znumpix, nAcq, nf, 'single'));
    end

    ee = (1 : numElements)';
    aa = reshape(1 : na, 1, na);

    %% Frames in chunks
    for f0 = 1 : Opts.frameChunk : nf
        fi = f0 : min(f0 + Opts.frameChunk - 1, nf);
        nfc = numel(fi);

        % Channel data of the chunk as baseband (interp) or analytic (nearest) IQ at fs: Dr2 (frames, elements * samples * acquisitions)
        Dr = complex(zeros(nfc, numElements, nIQ, nAcq, 'single'));
        if useGPU
            Dr = gpuArray(Dr);
        end
        m = (0 : nIQ-1)';
        carrier = single(exp(2j * pi * ((Fd / fs) * m - floor((Fd / fs) * m)))); % exp(+j 2 pi Fd m/fs), cycles wrapped
        for q = 1 : nAcq
            x = r(Receive(q).startSample : Receive(q).endSample, grpCh{1 + (q > na)}, fi);
            if size(x, 1) ~= ns
                error('Acquisition %d has %d samples, expected %d', q, size(x, 1), ns)
            end
            if useGPU, x = gpuArray(x); end
            a = toBaseband(single(x), bb, useGPU); % (nIQ, elements, frames)
            if ~interp
                a = a .* toDevice(carrier, useGPU); % analytic signal of the stored samples (what hilbert gave)
            end
            Dr(:, :, :, q) = permute(a, [3, 2, 1]);
        end
        Dr2 = reshape(Dr, nfc, []);
        clear Dr

        for iy = 1:ynumpix
            y = grid.y(iy);
            for iz = 1:znumpix
                z = grid.z(iz);

                % ---- group 1: TX columns steered in x-z, RX rows (y elements). dims (elements, angles, x)
                rx1 = sqrt(z^2 + (ey - y).^2) + d_correction_1way; % (E, 1)
                tx1 = z .* cos(angCols) + (xg + sign(angCols) .* x_end) .* sin(angCols) + d_correction_1way; % (1, na, nx)
                td1 = rx1 + tx1 + wl * tw_peak; % (E, na, nx)
                w1 = apodFun(atan((y - ey) ./ z)); % (E, 1)

                % ---- group 2: TX rows steered in y-z, RX columns (x elements)
                rx2 = sqrt(z^2 + (ex - reshape(xg, 1, 1, [])).^2) + d_correction_1way; % (E, 1, nx)
                tx2 = z .* cos(angG2) + (sign(angG2) .* y_end - y) .* sin(angG2) + d_correction_1way; % (1, na)
                td2 = rx2 + tx2 + wl * tw_peak; % (E, na, nx)
                w2 = apodFun(atan((reshape(xg, 1, 1, []) - ex) ./ z)); % (E, 1, nx)

                S = complex(zeros(nfc, nAcq, numel(xg), 'single')); % (frames, acquisitions, x)
                for g = 1:2
                    if g == 1, td = td1; w = w1; else, td = td2; w = w2; end
                    n = (td - 2 * startDepth * wl) / cTs; % 0-based fractional sample position
                    n = single(n);
                    if interp
                        i0 = floor(n);
                        wf = n - i0;
                        % the 1-based index of the sample at 0-based position i0 is i0 + 1 (zeroSample). RCA_DAS used i0
                        iA = i0 + double(Opts.zeroSample); iB = iA + 1;
                        ok = iA >= 1 & iB <= nIQ;
                        ph = single(exp(2j * pi * (0.25 * n - floor(0.25 * n)))) ; % carrier exp(+j 2 pi Fd n / fs), Fd/fs = 1/4
                        ca = (1 - wf) .* ph; cb = wf .* ph;
                    else
                        iA = round(n) + double(Opts.zeroSample);
                        ok = iA >= 1 & iA <= nIQ;
                        iB = iA; ca = ones(size(iA), 'single'); cb = zeros(size(iA), 'single');
                    end
                    ww = single(w) .* single(ok) + zeros(size(iA), 'single');
                    ca = ca .* ww; cb = cb .* ww;
                    iA(~ok) = 1; iB(~ok) = 1;
                    aq = (g - 1) * na + aa; % acquisition numbers (1, na)
                    base = ee + numElements * (nIQ * (reshape(aq, 1, na) - 1)); % (E, na)
                    LA = base + numElements * (iA - 1); LB = base + numElements * (iB - 1); % (E, na, nx)
                    if useGPU
                        LA = gpuArray(LA); LB = gpuArray(LB); ca = gpuArray(ca); cb = gpuArray(cb);
                    end
                    VA = reshape(Dr2(:, LA(:)), nfc, numElements, na, numel(xg));
                    val = VA .* reshape(ca, 1, numElements, na, numel(xg));
                    if interp
                        VB = reshape(Dr2(:, LB(:)), nfc, numElements, na, numel(xg));
                        val = val + VB .* reshape(cb, 1, numElements, na, numel(xg));
                    end
                    Sg = reshape(sum(val, 2), nfc, na, numel(xg)); % (frames, angles, x)
                    S(:, (g-1)*na + (1:na), :) = gather(Sg);
                end
                S = gather(S);
                if Opts.conjugateIQ, S = conj(S); end
                if Opts.sumAcq
                    IQ(:, iy, iz, 1, fi) = permute(sum(S, 2), [3, 1, 2]);
                else
                    IQ(:, iy, iz, :, fi) = permute(S, [3, 2, 1]);
                end
            end
        end
    end

end


%% Weight of a channel for the angle from the pixel to the element
function w = apodization(ang, mode, ES, esAngle, thetaCut)

    switch lower(mode)
        case 'sens' % RCA_DAS
            w = interp1(esAngle, ES, ang);
        case 'cutoff' % Verasonics' Recon.senscutoff
            w = double(abs(ang) <= thetaCut);
        case 'sens+cutoff'
            w = interp1(esAngle, ES, ang) .* (abs(ang) <= thetaCut);
        otherwise
            error('Unknown apod "%s"', mode)
    end

end


%% How the stored samples of a Receive relate to the demodulation frequency Fd (Vantage Sequence Programming Manual, 3.3.1)
% (copy of RCA_DAS_rcabeam's)
%   NS200BW: uniform, 4 samples per period of Fd (fsN = 4*Fd)
%   BS67BW:  uniform, 1 sample every 3/4 period (fsN = 4/3*Fd). The band of Fd is aliased and folded
%   BS100BW: pairs of samples, at t = m*Tp and m*Tp + 1/(4*Fd), with Tp = 1/Fd (the pair is a quadrature I, Q pair)
%   BS50BW:  same pairs with Tp = 2/Fd
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
    else
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
            bb.fsN = Fd; bb.pairs = true; bb.L = 4;
        case 'BS50BW'
            bb.fsN = Fd / 2; bb.pairs = true; bb.L = 8;
        otherwise
            error('Receive.sampleMode %s is not supported. Use NS200BW, BS67BW, BS100BW or BS50BW.', mode)
    end
    if bb.pairs
        if mod(nS, 2) ~= 0
            error('Receive has %d samples, but %s samples come in pairs', nS, mode)
        end
        bb.nIQ = bb.L * nS / 2;
    else
        bb.nIQ = bb.L * nS;
    end

end


%% Stored RF samples (nS, ...) -> complex baseband at carrier bb.Fd, bb.nIQ rows, at bb.fsEff (copy of RCA_DAS_rcabeam's)
% The first stored sample is at t = 0. Exact for a band limited signal (frequency domain filtering and interpolation)
function a = toBaseband(d, bb, useGPU)

    sz = size(d);
    nS = sz(1);
    d = reshape(d, nS, []);

    if bb.pairs
        M = nS / 2;
        tauSamples = (bb.fsN / (4 * bb.Fd));
        k = [0 : ceil(M/2) - 1, -floor(M/2) : -1]';
        shiftQ = toDevice(single(exp(-2j * pi * (k / M) * tauSamples)), useGPU);
        Q = real(ifft(fft(-d(2:2:end, :), [], 1) .* shiftQ, [], 1));
        a = complex(d(1:2:end, :), Q);
        Y = fft(a, [], 1);
        h = M / 2;
        Yu = zeros(bb.L * M, size(Y, 2), 'like', Y);
        Yu(1:h, :) = Y(1:h, :);
        Yu(end-h+1 : end, :) = Y(h+1 : end, :);
        Yu(h+1, :) = 0.5 * Y(h+1, :);
        Yu(end-h+1, :) = 0.5 * Y(h+1, :);
    else
        M = nS;
        ph = (0 : nS-1)' * (bb.Fd / bb.fsN);
        ph = ph - floor(ph);
        y = 2 * d .* toDevice(single(exp(-2j * pi * ph)), useGPU);
        Y = fft(y, [], 1);
        K = floor(M / 4);
        w = ones(2*K + 1, 1);
        if mod(M, 4) == 0
            w([1, end]) = 0.5;
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


function s = setDefault(s, name, value)

    if ~isfield(s, name)
        s.(name) = value;
    end

end
