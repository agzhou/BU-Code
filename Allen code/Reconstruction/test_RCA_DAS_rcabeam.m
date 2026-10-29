%% Description

% Checks RCA_DAS_rcabeam against a simulated point scatterer. The data is built directly from the RCA_DAS delay model
% (independent of rcabeam, same simulation as test_RCA_DAS_ffdas.m), so this checks the conventions the rcabeam interface
% leaves to the caller:
%   - angle and group mapping (RC/CR, flipped angles) and the shift of the transmit reference to the origin (t_start)
%   - sample index base and units (peak lands at the scatterer, within a pixel)
%   - demodulation to absolute time: the angles must add coherently ('absolute' vs. 'firstSample' phase reference)
%   - frame order and the zero padding of the last chunk (frames scaled 1, 2, 3 with frameChunk = 2)
%   - the sum of the RC and CR volumes equals the single pass sum (the sumOnly path)
%   - if the rcabeam backend is also run: the CUDA result against rcabeam_reference_beamform (a MATLAB port of its kernel)
% Set beforehand: backend = 'reference' (default, CPU, no Python needed) or 'rcabeam' (GPU, needs R2026b and the Python
% environment, see RCA_DAS_rcabeam). With 'rcabeam' the reference is run too and the two are compared.
% Optional presets: pythonExe, fNumber, numAngles, maxAngleDeg, N (# elements per side), pix (pixel size in m)

c = 1540; fc = 15.625e6; wl = c / fc; spw = 4; fs = fc * spw; dps = c / fs;
if exist('N', 'var') ~= 1, N = 32; end
pitch = 0.1e-3; nf = 3;
startDepth = 10; endDepth = 40; twPeak = 1; % wavelengths
if exist('numAngles', 'var'), na = numAngles; else, na = 11; end
if exist('maxAngleDeg', 'var'), maxAng = maxAngleDeg; else, maxAng = 15; end
if exist('backend', 'var') ~= 1, backend = 'reference'; end
if exist('fNumber', 'var') ~= 1, fNumber = 0; end
if exist('pythonExe', 'var') ~= 1, pythonExe = ''; end
if exist('pix', 'var') ~= 1, pix = 0.1e-3; end
angles = linspace(-maxAng, maxAng, na) * pi/180;

xe = ((0:N-1) - (N-1)/2) * pitch; ye = xe;
x_end = max(xe); y_end = max(ye);

% Minimal P with the fields RCA_DAS_rcabeam reads
P.Trans.frequency = fc / 1e6;
P.Trans.spacingMm = pitch * 1e3;
P.Trans.lensCorrection = 0;
P.Trans.Connector = (1 : 2*N)';
P.Trans.ElementPos = zeros(2*N, 4);
P.Trans.ElementPos(1:N, 1) = xe / wl;
P.Trans.ElementPos(N+1 : 2*N, 2) = ye / wl;
P.na = na; P.numElements = N; P.angles = angles;
P.startDepth = startDepth; P.endDepth = endDepth;
P.TW.peak = twPeak;
P.Resource.Parameters.speedOfSound = c;
P.numFramesPerBuffer = nf;

nS = 2 * (endDepth - startDepth) * spw; % samples per acquisition
for a = 1 : 2*na
    P.Receive(a).startSample = (a-1) * nS + 1;
    P.Receive(a).endSample = a * nS;
    P.Receive(a).samplesPerWave = spw;
end

%% Simulate a point scatterer (off-center in x and y), placed on a voxel center of the grid RCA_DAS_rcabeam builds so the
% peak is not reduced by being between voxels. x, y: centers are at (i - 0.5 - n/2) * pix, z: see the grid in RCA_DAS_rcabeam
zsize = (endDepth - startDepth) * wl; nzp = ceil(zsize / pix);
scat = [(25 - 0.5 - N/2) * pitch, (10 - 0.5 - N/2) * pitch, startDepth * wl + (16 - 1/2) / nzp * zsize]; % x, y, z in m
if pix ~= pitch
    warning('The scatterer is only on a voxel center for pix = pitch, so the peak checks below can fail')
end
j = (0 : nS-1)';
sigma = 1.5 * spw;
pulse = @(n0) exp(-0.5 * ((j - n0) / sigma).^2) .* cos(2*pi * (j - n0) / spw); % real RF, peak at (fractional) sample n0

rf = zeros(nS * 2*na, 2*N, nf);
for a = 1:na
    % Group 1: TX columns steered in x-z, RX rows (positioned along y)
    th = angles(a);
    tx = scat(3) * cos(th) + (scat(1) + sign(th) * x_end) * sin(th);
    rx = sqrt(scat(3)^2 + (scat(2) - ye).^2);
    n0 = (tx + rx + wl * twPeak - 2 * startDepth * wl) / dps; % 0-based sample index of the echo peak
    for ch = 1:N
        rf((a-1)*nS + (1:nS), N + ch, :) = repmat(pulse(n0(ch)), 1, 1, nf);
    end

    % Group 2: TX rows steered in y-z, RX columns (positioned along x). The angle order is flipped
    th = angles(na + 1 - a);
    tx = scat(3) * cos(th) + (sign(th) * y_end - scat(2)) * sin(th);
    rx = sqrt(scat(3)^2 + (scat(1) - xe).^2);
    n0 = (tx + rx + wl * twPeak - 2 * startDepth * wl) / dps;
    for ch = 1:N
        rf((na + a - 1)*nS + (1:nS), ch, :) = repmat(pulse(n0(ch)), 1, 1, nf);
    end
end
rf = rf .* reshape(1:nf, 1, 1, nf); % frame f is f times frame 1

%% Reconstruct and locate the peaks
base = struct('frameChunk', 2, 'fNumber', fNumber, 'pythonExe', pythonExe, 'backend', backend);
fprintf('Scatterer at x = %.0f, y = %.0f, z = %.0f um, backend %s\n', scat * 1e6, backend)
[IQ_CR, IQ_RC, grid] = RCA_DAS_rcabeam(rf, P, pix, pix, pix, base);

vols = {IQ_CR, IQ_RC};
names = {'C-R', 'R-C'};
fprintf('%-6s %10s %8s %8s %8s %12s\n', 'vol', 'peak', 'dx um', 'dy um', 'dz um', 'frame ratios')
posOK = true; ratioOK = true; gainOK = true; peak = zeros(1, 2); at = zeros(2, 3);
for v = 1:2
    [m, ind] = max(reshape(abs(vols{v}(:, :, :, 1)), [], 1));
    [ix, iy, iz] = ind2sub(size(vols{v}, 1:3), ind);
    err = [grid.x(ix) - scat(1), grid.y(iy) - scat(2), grid.z(iz) - scat(3)] * 1e6;
    pk = squeeze(max(max(max(abs(vols{v}), [], 1), [], 2), [], 3))';
    ratios = pk / pk(1);
    fprintf('%-6s %10.3g %8.1f %8.1f %8.1f %12s\n', names{v}, m, err, mat2str(ratios, 4))
    posOK = posOK && all(abs(err) <= pix * 1e6 / 2);
    ratioOK = ratioOK && all(abs(ratios - (1:nf)) < 1e-3 * (1:nf));
    gainOK = gainOK && m > 0.95 * na * N; % all angles and channels add in phase (pulse amplitude 1)
    peak(v) = m;
    at(v, :) = [ix, iy, iz];
end

%% Coherence: the demodulation phase reference, compared at the voxel of the scatterer
o = base; o.phaseReference = 'firstSample';
[a1, a2] = RCA_DAS_rcabeam(rf, P, pix, pix, pix, o);
peakFirst = [abs(a1(at(1, 1), at(1, 2), at(1, 3), 1)), abs(a2(at(2, 1), at(2, 2), at(2, 3), 1))];

%% The RC+CR sum in one pass
o = base; o.sumOnly = true;
S = RCA_DAS_rcabeam(rf, P, pix, pix, pix, o);
sumErr = norm(S(:) - IQ_CR(:) - IQ_RC(:)) / norm(S(:));

fprintf('\nPeak on the scatterer voxel (within half a pixel, %.0f um): %s\n', pix * 1e6 / 2, passFail(posOK))
fprintf('Peak above 95%% of the ideal coherent sum (%d angles x %d channels = %d), C-R %.3g, R-C %.3g: %s\n', ...
    na, N, na * N, peak, passFail(gainOK))
fprintf('Frame amplitudes scale as 1:2:3 (frame order, padded last chunk): %s\n', passFail(ratioOK))
fprintf('At that voxel, absolute vs. first sample phase reference (should be much higher for absolute): C-R %.3g vs. %.3g, R-C %.3g vs. %.3g -> %s\n', ...
    peak(1), peakFirst(1), peak(2), peakFirst(2), passFail(all(peak > 2 * peakFirst)))
fprintf('RC + CR vs. single pass sum, relative difference %.2g: %s\n', sumErr, passFail(sumErr < 1e-4))

%% Against the MATLAB port of the CUDA kernel, when the GPU backend was used
if strcmpi(backend, 'rcabeam')
    o = base; o.backend = 'reference';
    [R_CR, R_RC] = RCA_DAS_rcabeam(rf, P, pix, pix, pix, o);
    e1 = norm(IQ_CR(:) - R_CR(:)) / norm(R_CR(:));
    e2 = norm(IQ_RC(:) - R_RC(:)) / norm(R_RC(:));
    fprintf('rcabeam (CUDA) vs. rcabeam_reference_beamform, relative RMS difference C-R %.2g, R-C %.2g: %s\n', e1, e2, ...
        passFail(max(e1, e2) < 2e-3))
    RCA_DAS_rcabeam('clear')
end


function s = passFail(tf)
    if tf
        s = 'PASS';
    else
        s = 'FAIL';
    end
end
