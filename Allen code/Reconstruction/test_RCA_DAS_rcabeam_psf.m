%% Description

% Checks the GPU rcabeam backend of RCA_DAS_rcabeam on the real simulated point scatterer data (the same RF file that
% plot_PSF_rcabeam.m uses) against rcabeam_reference_beamform, the MATLAB port of the CUDA kernel. test_RCA_DAS_rcabeam.m does
% this on a synthetic scatterer; this adds the real P (Receive start samples, Trans.Connector, lens correction, TW.peak, the
% actual angle list and sampling) and a real-sized volume. The two backends get exactly the same data preparation, so a
% mismatch is in the GPU hand-over (layout, frame chunking, workspace caching) or the kernel, not in the recon setup.
% Checks:
%   - C-R and R-C volumes of the GPU against the reference, per frame (relative RMS and peak difference)
%   - the GPU volumes were made by the GPU (backend, and the second call is faster than the first, which builds the workspace)
%   - frameChunk = 1 (a different chunking and padding than the default) gives the same volumes
%   - the sumOnly pass equals C-R + R-C
%   - peak of the OPW sum is at the scatterer, within a voxel (informational: it also tests the recon model, not just the GPU)

% Set beforehand (all optional):
%   psfPath      folder with params.mat and the RF file(s)
%   pix          voxel size in m (default one wavelength: the reference runs on the CPU, and the half wavelength grid of
%                plot_PSF_rcabeam.m takes minutes)
%   nFrames      # frames to reconstruct (default 2, capped at what the file has. The RC15gV PSF file has one frame, so the
%                frame order and chunking are only exercised by test_RCA_DAS_rcabeam.m, or by a multi-frame file here)
%   fNumber      receive f-number (default 0, full aperture). Try 1.5 too, the aperture mask is a separate code path
%   pythonExe    python.exe with cupy and rcabeam, if RCABEAM_PYTHON is not set
%   gpuBackend   'rcabeam' (default). 'reference' makes both runs the CPU reference, to check this script without a GPU
%   tolerance    pass threshold of the relative RMS difference (default 2e-3, like test_RCA_DAS_rcabeam.m)

if exist('psfPath', 'var') ~= 1
    psfPath = 'G:\My Drive\Data\PSF Simulations\RC15gV PSF sim - 5 angles from -6 to 6 deg';
end
if exist('nFrames', 'var') ~= 1, nFrames = 2; end
if exist('fNumber', 'var') ~= 1, fNumber = 0; end
if exist('pythonExe', 'var') ~= 1, pythonExe = ''; end
if exist('gpuBackend', 'var') ~= 1, gpuBackend = 'rcabeam'; end
if exist('tolerance', 'var') ~= 1, tolerance = 2e-3; end

S = load(fullfile(psfPath, 'params.mat'), 'P'); P = S.P;
if ~isfield(P, 'numFramesPerBuffer') && isfield(P, 'numFramesPerBuff'), P.numFramesPerBuffer = P.numFramesPerBuff; end
wl = P.Resource.Parameters.speedOfSound / (P.Trans.frequency * 1e6);
if exist('pix', 'var') ~= 1, pix = wl; end
rfFile = dir(fullfile(psfPath, 'RF-*.mat')); rfFile = rfFile(1).name;

R = load(fullfile(psfPath, rfFile), 'RcvData');
r = unstackFrames(R.RcvData, P);
clear R
nFrames = min(nFrames, size(r, 3));
r = r(:, :, 1 : nFrames);
if isfield(P, 'numSubFrames'), P.numSubFrames = nFrames; end % RCA_DAS_rcabeam takes the frame count from here if it exists
fprintf('%s, %d frame(s), voxel %.0f um (wavelength %.0f um), fNumber %g, GPU backend "%s"\n', ...
    rfFile, nFrames, pix * 1e6, wl * 1e6, fNumber, gpuBackend)

base = struct('fNumber', fNumber, 'pythonExe', pythonExe);

%% Reference (CPU)
o = base; o.backend = 'reference'; o.frameChunk = nFrames;
tic
[R_CR, R_RC, grid] = RCA_DAS_rcabeam(r, P, pix, pix, pix, o);
tRef = toc;
fprintf('Reference: %.0f s, volume %s\n', tRef, mat2str(size(R_CR)))

%% GPU, default chunking. The first call builds the workspace, the second reuses it
o = base; o.backend = gpuBackend; o.frameChunk = nFrames;
RCA_DAS_rcabeam('clear');
tic; [G_CR, G_RC] = RCA_DAS_rcabeam(r, P, pix, pix, pix, o); tFirst = toc;
tic; [G_CR2, G_RC2] = RCA_DAS_rcabeam(r, P, pix, pix, pix, o); tSecond = toc;
fprintf('%s: first call %.1f s, second call %.1f s\n', gpuBackend, tFirst, tSecond)

%% GPU, one frame per chunk
o.frameChunk = 1;
[C_CR, C_RC] = RCA_DAS_rcabeam(r, P, pix, pix, pix, o);

%% GPU, one pass sum
o = base; o.backend = gpuBackend; o.frameChunk = nFrames; o.sumOnly = true;
S_sum = RCA_DAS_rcabeam(r, P, pix, pix, pix, o);

%% Compare
rel = @(a, b) norm(a(:) - b(:)) / norm(b(:));
fprintf('\n%-34s %12s %12s\n', 'relative RMS difference to the reference', 'C-R', 'R-C')
fprintf('%-34s %12.2g %12.2g\n', 'default chunking', rel(G_CR, R_CR), rel(G_RC, R_RC))
fprintf('%-34s %12.2g %12.2g\n', 'frameChunk = 1', rel(C_CR, R_CR), rel(C_RC, R_RC))
for f = 1 : nFrames
    fprintf('%-34s %12.2g %12.2g\n', sprintf('  frame %d (default chunking)', f), ...
        rel(G_CR(:, :, :, f), R_CR(:, :, :, f)), rel(G_RC(:, :, :, f), R_RC(:, :, :, f)))
end
pk = @(v) max(abs(v(:)));
fprintf('\nPeak |IQ|, reference vs. GPU: C-R %.4g vs. %.4g, R-C %.4g vs. %.4g\n', pk(R_CR), pk(G_CR), pk(R_RC), pk(G_RC))

e = max([rel(G_CR, R_CR), rel(G_RC, R_RC), rel(C_CR, R_CR), rel(C_RC, R_RC)]);
okGPU = e < tolerance;
okRepeat = isequal(G_CR, G_CR2) && isequal(G_RC, G_RC2);
okNonzero = pk(R_CR) > 0 && pk(G_CR) > 0 && pk(R_RC) > 0 && pk(G_RC) > 0; % an empty-vs-empty comparison would pass
okSum = rel(S_sum, G_CR + G_RC) < 1e-4;
okSpeed = strcmpi(gpuBackend, 'reference') || tSecond < tFirst;

%% Peak of the OPW sum against the scatterer (x = 0, y = 0, z = P.Media.MP(3) wavelengths)
okPeak = true;
if isfield(P, 'Media') && isfield(P.Media, 'MP')
    scat = [0, 0, P.Media.MP(3) * wl];
    A = abs(G_CR(:, :, :, 1) + G_RC(:, :, :, 1));
    [~, i] = max(A(:));
    [ix, iy, iz] = ind2sub(size(A), i);
    err = [grid.x(ix) - scat(1), grid.y(iy) - scat(2), grid.z(iz) - scat(3)];
    okPeak = all(abs(err) <= pix);
    fprintf('Peak of the GPU OPW sum is %.0f, %.0f, %.0f um from the scatterer (limit one voxel, %.0f um)\n', err * 1e6, pix * 1e6)
end

fprintf('\nGPU vs. reference, relative RMS %.2g (limit %.0e), all outputs: %s\n', e, tolerance, passFail(okGPU))
fprintf('Volumes are not empty: %s\n', passFail(okNonzero))
fprintf('Second call gives the same volumes as the first (workspace reuse): %s\n', passFail(okRepeat))
fprintf('Second call faster than the first (workspace was kept): %s\n', passFail(okSpeed))
fprintf('sumOnly pass equals C-R + R-C: %s\n', passFail(okSum))
fprintf('Peak within one voxel of the scatterer: %s\n', passFail(okPeak))

RCA_DAS_rcabeam('clear');


function s = passFail(tf)
    if tf
        s = 'PASS';
    else
        s = 'FAIL';
    end
end
