%% Description

% Checks that RCA_DAS_rcabeam gives the same IQ from the compressed RF of compressRCA_RF as from the full unstacked RF.
% Random int16 RF with a scrambled Trans.Connector and a gap of unused samples between acquisitions, so a wrong channel
% map, acquisition order or sample offset shows up. Also checks a single frame file (MATLAB drops the 4th dim then).
% Set backend = 'rcabeam' beforehand to also run it on the GPU (default 'reference', CPU, no Python needed).

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'Ultrasound data acquisition', 'RFcompression'))
if exist('backend', 'var') ~= 1, backend = 'reference'; end

c = 1540; fc = 15.625e6; wl = c / fc; spw = 4;
N = 8; na = 3; nf = 3; pitch = 0.1e-3; pix = 0.2e-3;
startDepth = 10; endDepth = 20; % wavelengths
nS = 2 * (endDepth - startDepth) * spw; % samples per acquisition
gap = 16;                                % unused samples after each acquisition
nCh = 2*N + 4;                           % more receive channels than elements, like a system with spare channels

P.Trans.name = 'RC15gV';
P.Trans.frequency = fc / 1e6;
P.Trans.spacingMm = pitch * 1e3;
P.Trans.lensCorrection = 0;
rng(1)
P.Trans.Connector = randperm(nCh, 2*N)';
xe = ((0:N-1) - (N-1)/2) * pitch;
P.Trans.ElementPos = zeros(2*N, 4);
P.Trans.ElementPos(1:N, 1) = xe / wl;
P.Trans.ElementPos(N+1 : 2*N, 2) = xe / wl;
P.na = na; P.numElements = N; P.angles = linspace(-10, 10, na) * pi/180;
P.startDepth = startDepth; P.endDepth = endDepth;
P.TW.peak = 1;
P.Resource.Parameters.speedOfSound = c;
for a = 1 : 2*na
    P.Receive(a).startSample = (a-1) * (nS + gap) + 1;
    P.Receive(a).endSample = P.Receive(a).startSample + nS - 1;
    P.Receive(a).samplesPerWave = spw;
end
nspf = 2*na * (nS + gap);

r = int16(randi([-2000, 2000], nspf, nCh, nf));
[rc, layout] = compressRCA_RF(r, P);

opts = struct('backend', backend, 'grid', 'centered', 'fNumber', 0, 'frameChunk', 2);

%% Layout: each acquisition's slice is the right samples of the right channels
okLayout = isequal(size(rc), [nS, 2*na, N, nf]) && isa(rc, 'int16');
for a = 1 : 2*na
    if a <= na, el = N+1 : 2*N; else, el = 1 : N; end
    ref = r(P.Receive(a).startSample : P.Receive(a).endSample, P.Trans.Connector(el), :);
    okLayout = okLayout && isequal(reshape(rc(:, a, :, :), nS, N, nf), ref);
end
fprintf('Compressed layout matches the remapped slices of the full RF: %s\n', passFail(okLayout))

%% Full vs. compressed, RC and CR separately
[CR1, RC1] = RCA_DAS_rcabeam(r, P, pix, pix, pix, opts);
[CR2, RC2] = RCA_DAS_rcabeam(rc, P, pix, pix, pix, opts);
okRCCR = isequal(CR1, CR2) && isequal(RC1, RC2) && norm(CR1(:)) > 0;
fprintf('RC and CR volumes identical (full vs. compressed): %s\n', passFail(okRCCR))

%% sumOnly (the OPW output of the recon)
o = opts; o.sumOnly = true;
S1 = RCA_DAS_rcabeam(r, P, pix, pix, pix, o);
S2 = RCA_DAS_rcabeam(rc, P, pix, pix, pix, o);
okSum = isequal(S1, S2);
fprintf('OPW sum identical: %s\n', passFail(okSum))

%% One frame (rc is 3-D then)
S1 = RCA_DAS_rcabeam(r(:, :, 2), P, pix, pix, pix, o);
S2 = RCA_DAS_rcabeam(rc(:, :, :, 2), P, pix, pix, pix, o);
okOne = ndims(rc(:, :, :, 2)) == 3 && isequal(size(S1), size(S2)) && isequal(S1, S2);
fprintf('Single frame identical: %s\n', passFail(okOne))

if strcmpi(backend, 'rcabeam'), RCA_DAS_rcabeam('clear'), end
if ~(okLayout && okRCCR && okSum && okOne)
    error('test_compressRCA_RF failed')
end


function s = passFail(tf)
    if tf, s = 'PASS'; else, s = 'FAIL'; end
end
