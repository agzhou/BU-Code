%% Description

% Compares RCA_DAS_fixed (a copy of the handwritten RCA_DAS with the zero sample, phase and sample mode fixes) and rcabeam
% (RCA_DAS_rcabeam) with the Verasonics recon, on a real RC15gV acquisition. RCA_DAS_fixed is run with the fixes turned on one at
% a time, starting from RCA_DAS's behavior (Opts of RCA_DAS_fixed, rows A to D), and rcabeam is run with its defaults. For each
% the IQ of a few xz planes goes through the clutter filter of IQ2g1T_3D_loop.m (SVD cut, 25 Hz high pass) and is compared with
% the Verasonics IQ of the same plane: complex correlation of the filtered IQ per voxel (median), correlation of the PDI
% divided by the median PDI of its depth, and the vessel to background contrast of that PDI. The last rows compare rcabeam and
% the fully fixed RCA_DAS_fixed with each other.
% The Verasonics recon is the folder 'IQ Data - Verasonics recon' (IData/QData, [y x z 1 frames], double) of runPath.
% Note that the correlation of filtered IQ cannot reach 1 even for identical recons: it drops to ~0.98 with 1e-6 relative noise
% on the IQ, since the SVD filter's high order components are sensitive to the noise level

% Set beforehand (all optional):
%   runPath    folder with params.mat, the RF files and the Verasonics IQ
%   fileNum    RF/IQ file number (1)
%   rows       rows of the PData y grid to compare (41, the middle; 25 and 60 have more vessels)
%   nFrames    # of frames (all of the file)
%   svRange    [lower, upper] singular values of the clutter filter ([50, nFrames])
%   zRange     z voxels to compare (20:142)

if exist('runPath', 'var') ~= 1
    runPath = 'E:\Ultrasound data from 7-21-2025\08-28-2025 AZ01 RCA fUS air puff\run 1';
end
if exist('fileNum', 'var') ~= 1, fileNum = 1; end
if exist('rows', 'var') ~= 1, rows = 41; end
if exist('zRange', 'var') ~= 1, zRange = 20:142; end

S = load(fullfile(runPath, 'params.mat'), 'P'); P = S.P;
if exist('nFrames', 'var') ~= 1, nFrames = P.numFramesPerBuffer; end
if exist('svRange', 'var') ~= 1, svRange = [50, nFrames]; end
wl = P.Resource.Parameters.speedOfSound / (P.Trans.frequency * 1e6);
pix = [P.Trans.spacing, P.Trans.spacing, 0.5] * wl; % the Verasonics PDelta of the fUS runs, [x, y, z] in m
nS = P.Receive(1).endSample; nspf = nS * P.na * 2;
name = sprintf('%d-%d-%d-%d-1-%d.mat', round(P.maxAngle), P.na, P.frameRate, P.numFramesPerBuffer, fileNum);

fprintf('Loading the Verasonics IQ and the RF data (%d frames)...\n', nFrames)
mv = matfile(fullfile(runPath, 'IQ Data - Verasonics recon', ['IQ-', name]));
V = cell(1, numel(rows));
for q = 1:numel(rows)
    V{q} = single(squeeze(complex(mv.IData(rows(q), :, :, 1, 1:nFrames), mv.QData(rows(q), :, :, 1, 1:nFrames)))); % [x, z, frames]
end
m = matfile(fullfile(runPath, ['RF-', name]));
r = permute(reshape(m.RcvData(1 : nspf * nFrames, :), nspf, nFrames, []), [1, 3, 2]); % unstacked

%% Reconstructions of the rows
cfg = {
    'A  RCA_DAS logic: analytic signal, nearest sample, 1 sample early, ElementSens weights', ...
        struct('iq', 'analytic', 'zeroSample', false, 'apod', 'sens')
    'B  A + zero sample fix', ...
        struct('iq', 'analytic', 'zeroSample', true, 'apod', 'sens')
    'C  B + baseband, fractional sample, carrier phase of the exact delay', ...
        struct('iq', 'baseband', 'zeroSample', true, 'apod', 'sens')
    'D  C + Verasonics aperture (Recon.senscutoff), the RCA_DAS_fixed defaults', ...
        struct()
};
IQ = cell(size(cfg, 1) + 1, 1);
for k = 1:size(cfg, 1)
    o = cfg{k, 2}; o.yIndex = rows; o.sumAcq = true;
    tic; X = RCA_DAS_fixed(r, P, pix(2), pix(1), pix(3), o);
    fprintf('%s: %.0f s\n', cfg{k, 1}(1:1), toc)
    IQ{k} = squeeze(X(:, :, :, 1, :)); % [x, rows, z, frames]
end
o = struct('sumOnly', true, 'yIndex', rows, 'frameChunk', 16);
tic; X = RCA_DAS_rcabeam(r, P, pix(2), pix(1), pix(3), o); fprintf('rcabeam: %.0f s\n', toc)
RCA_DAS_rcabeam('clear');
IQ{end} = squeeze(X); clear X r

%% Comparison with the Verasonics IQ, per row
names = [cellfun(@(c) c(1:1), cfg(:, 1), 'UniformOutput', false); {'rcabeam'}];
addpath(genpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'Processing'))) % getSVs2D, applySVs2D_accel
for q = 1:numel(rows)
    fprintf('\nrow %d, filtered IQ (SVD %d-%d, 25 Hz high pass) vs Verasonics:\n', rows(q), svRange)
    fprintf('  %-9s %-26s %-26s %s\n', '', 'filtered IQ corr (top 10%)', 'flat PDI corr', 'vessel/background (Verasonics)')
    pv = clutterFilter(V{q}, P, svRange, zRange);
    fv = pv.PDI ./ median(pv.PDI, 1);
    top = fv > prctile(fv(:), 95); bg = fv < prctile(fv(:), 50);
    P_ = cell(1, numel(IQ));
    for k = 1:numel(IQ)
        if ndims(IQ{k}) == 4, Xk = squeeze(IQ{k}(:, q, :, :)); else, Xk = IQ{k}; end
        P_{k} = clutterFilter(Xk, P, svRange, zRange);
        [rho, rhoTop] = filtCorr(pv.Xf, P_{k}.Xf);
        fk = P_{k}.PDI ./ median(P_{k}.PDI, 1);
        fprintf('  %-9s %.3f (%.3f)%14s %.3f%19s %.2f (%.2f)\n', names{k}, rho, rhoTop, '', corr(fv(:), fk(:)), '', ...
            mean(fk(top)) / mean(fk(bg)), mean(fv(top)) / mean(fv(bg)))
    end
    fprintf('  between the recons: filtered IQ correlation of rcabeam with\n')
    for k = 1:numel(IQ) - 1
        [rho, rhoTop] = filtCorr(P_{end}.Xf, P_{k}.Xf);
        fprintf('    %-7s %.3f (%.3f)\n', names{k}, rho, rhoTop)
    end
end


%% Clutter filter as in IQ2g1T_3D_loop.m: [x, z, frames] -> filtered IQ and PDI (x, z)
function out = clutterFilter(IQ, P, sv, zr)

    X = reshape(IQ(:, zr, :), size(IQ, 1), 1, numel(zr), size(IQ, 3));
    [b, a] = butter(4, 25 / (P.frameRate / 2), 'high');
    [CM, EVs, Vs] = getSVs2D(X);
    Xf = filter(b, a, applySVs2D_accel(X, CM, EVs, Vs, sv(1), min(sv(2), size(IQ, 3))), [], 4);
    out.Xf = squeeze(Xf);
    out.PDI = mean(abs(out.Xf) .^ 2, 3);

end


%% Median complex correlation over time per voxel (and the median over the 10% of voxels with the most energy of a)
function [rho, rhoTop] = filtCorr(a, b)

    r = abs(sum(a .* conj(b), 3)) ./ sqrt(sum(abs(a) .^ 2, 3) .* sum(abs(b) .^ 2, 3));
    e = sum(abs(a) .^ 2, 3);
    rho = median(r(:));
    rhoTop = median(r(e > prctile(e(:), 90)));

end
