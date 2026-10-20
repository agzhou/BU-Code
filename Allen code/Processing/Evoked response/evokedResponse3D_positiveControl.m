%% Description: positive control (synthetic injection) for the 3D evoked response processing in evokedResponse3D*.m
%   A known, HRF-shaped, stimulus-locked change is injected into ROI voxels of the REAL PDI data (so the real noise is
%   kept), then the same analysis functions used in evokedResponse3D.m (periStimulusAverage, glmActivationMap, and the legacy
%   interp1 + 'extrap' trial-average path) are run on it. If the injected response is recovered, the processing chain works and a
%   missing response in the real data means the response is below the noise floor; if it is not recovered, the chain is broken.
%
%   Tests:
%     A. Peri-stimulus averaging: recovered ROI-mean evoked response vs the expected one (injected - uninjected, so the noise cancels)
%     B. Voxelwise GLM: noise-free fidelity (beta/amp), then detection vs injected amplitude (raw and spatially smoothed data),
%        sensitivity to the stim timing offset and to the analysis HRF peak. Also reports the noise floor of the real data
%        (beta standard error in % PDI) = the smallest response you can expect to detect.
%     C. Legacy path (evokedResponse3D.m "Resample the data": interp1 with 'extrap' -> trial average -> correlation), with and
%        without extrapolation.
%
%   Injection model: PDI_inj(t) = PDI(t) * (1 + amp * reg(t)), reg = the exact GLM regressor (HRF-convolved stim, averaged over each
%   superframe), so amp = the GLM beta / 100 for a sustained stimulus (e.g. 0.05 -> +5 %). NOTE: this scales the TOTAL PDI. A real
%   CBV change scales only the blood signal, so if a fraction f of PDI is blood signal (see 1 - noise fraction), the equivalent
%   physiological change is amp / f.
%
%   Timing is derived from the measured data (RF time tags, stim trace), not from the nominal acquisition parameters.
%
%   To override options, define a struct `opt` with the fields below before running (e.g. opt.nROIs = 4).

%% 0. Options
def = struct();
def.jUse         = 3;               % Which PDI cell to use if fUSdata's PDI is a cell (1: neg, 2: pos, 3: all freq.). Ignored if PDI is numeric (e.g. XDoppler)
def.amps         = [0.01, 0.02, 0.05, 0.10];  % Injected amplitudes (fractional change of PDI for a sustained stim)
def.ampSweep     = 0.05;            % Amplitude used for the timing / HRF sweeps, the evoked plots and the legacy test
def.injHRFPeak   = 3;               % HRF peak time of the injected response [s]
def.hrfPeaks     = [1, 2, 3, 4, 6]; % Analysis HRF peak times to sweep [s]
def.timeOffsets  = [-2, -1, 0, 1, 2]; % Offsets of the analysis stim relative to the injected one [s]
def.smoothSigma  = 1.5;             % Gaussian smoothing sigma for the 'smoothed' analysis [voxels]
def.roiMode      = 'random';        % 'random': nROIs random brain ROIs. 'click': also click a barrel-cortex ROI on the MIPs. 'preset': use roiCenterVox as the first ROI(s)
def.nROIs        = 8;               % Number of random ROIs (injected simultaneously, non-overlapping)
def.roiSemi_mm   = [0.5, 0.5, 0.4]; % ROI ellipsoid semi-axes along dims 1, 2, 3 [mm]
def.roiCenterVox = [];              % (n x 3) ROI centres in voxels (dims 1-3), for 'preset'
def.brainThrFrac = 0.3;             % Brain mask: PDI average > this fraction of its 99th percentile
def.seed         = 1;               % RNG seed (ROI placement, background subsets)
def.runLegacy    = true;            % Run test C (legacy resampling path)
def.legacyNBg    = 10000;          % # of background voxels used in test C
def.showFigs     = true;
def.savePath     = '';              % Folder to save the results struct in ('' = don't save)
if exist('opt', 'var') && isstruct(opt)
    fn = fieldnames(opt);
    for k = 1:numel(fn), def.(fn{k}) = opt.(fn{k}); end
end
opt = def; clearvars def fn k
rng(opt.seed)
if ~opt.showFigs, set(0, 'DefaultFigureVisible', 'off'); end

%% 1. Paths and data
thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir), thisDir = cd; end
addpath(fileparts(thisDir), thisDir) % Processing folder (getReconCoords3D) and the Evoked response folder

if ~exist('PDpath', 'var')
    PDpath = uigetdir('D:\Allen\Data\', 'Select the Processed data (PDI) path');
    PDpath = [PDpath, filesep];
end
runPath = [fileparts(PDpath(1:end-1)), filesep]; % Run folder (params.mat, RFTimeTagsPerChannel.mat, triggerData.mat, PData.mat)

if ~exist('P', 'var'),          load(findFile({runPath, PDpath}, 'params.mat'), 'P'), end
if ~exist('RFtimeTags', 'var'), load(findFile({runPath, PDpath}, 'RFTimeTagsPerChannel.mat'), 'RFtimeTags', 'RFcount'), end
if ~exist('inScanData', 'var'), load(findFile({runPath, PDpath}, 'triggerData.mat'), 'inScanData'), end
if ~exist('PData', 'var'),      load(findFile({PDpath, runPath}, 'PData.mat'), 'PData'), end
if ~exist('RFcount', 'var'),    RFcount = size(RFtimeTags, 3); end

%% 2. Timing from the measured data
analogHigh = 2; % Trigger value past which we treat a signal as 'on'
fs = P.daqrate;

% -- Superframe times from the RF time tags (the frame rate and superframe period are NOT the nominal values) --
%   median over channels: robust to channels without a valid time tag
tt = squeeze(median(RFtimeTags, 2));  % (# frames, # superframes)
RFTT = tt - tt(1, 1);                 % t = 0 at the first acquired frame
chSpread = max(max(RFtimeTags(1, :, :), [], 2) - min(RFtimeTags(1, :, :), [], 2), [], 'all');
if chSpread > 1e-3, warning('Time tags differ across channels by up to %.3g s; check RFtimeTags', chSpread), end
sfStarts = RFTT(1, :).';              % Superframe start times [s]
frameDt = median(diff(RFTT, 1, 1), 'all'); % Real frame spacing within a superframe [s]
sfWidth = size(RFTT, 1) * frameDt;    % Real superframe acquisition duration [s]
sfCenters = sfStarts + sfWidth / 2;
nSF = numel(sfStarts);
sfPeriod = median(diff(sfStarts));

% -- Stim onsets from the recorded stim trace, relative to the Verasonics trigger (= acquisition start) --
if isfield(P, 'Mcr_fcp') && isfield(P.Mcr_fcp, 'vts') && isfield(P.Mcr_fcp.vts, 'signal')
    indTrig = find(P.Mcr_fcp.vts.signal > analogHigh, 1, 'first'); % Rising edge of the programmed Verasonics trigger [DAQ sample]
else
    indTrig = round(P.vts.delay_s * fs) + 1;
end
stimTrace = inScanData(:, 1);
edgesOn = find(diff([false; stimTrace > analogHigh]) == 1);          % Rising edges of the (3 Hz) puff train
trialStartIdx = edgesOn([true; diff(edgesOn) > 1 * fs]);              % First puff of each trial
if numel(trialStartIdx) ~= P.numTrials
    warning('Found %d trials in the stim trace, expected %d; using the first onset + the nominal trial period', numel(trialStartIdx), P.numTrials)
    trialStartIdx = trialStartIdx(1) + round((0:P.numTrials-1).' * P.apis.seq_length_s * fs);
end
stimOnsetTimestamps = (trialStartIdx(:) - indTrig) / fs;             % [s] in the acquisition time base
stimEndTimestamps = stimOnsetTimestamps + P.apis.stim_length_s;
onsetInTrial_s = median(mod(stimOnsetTimestamps, P.apis.seq_length_s)); % Stim onset within a trial window that starts at t = 0 [s]
fprintf('Measured timing: first stim onset %.3f s after the trigger (nominal %.3f s); trial period %.4f s (nominal %.1f s)\n', ...
    stimOnsetTimestamps(1), P.apis.delay_time_ms/1e3, median(diff(stimOnsetTimestamps)), P.apis.seq_length_s)
fprintf('Real RF timing: superframe period %.4f s (nominal %.4f s), frame spacing %.4f ms (nominal %.4f ms), superframe width %.4f s (nominal %.4f s)\n', ...
    sfPeriod, 1/P.sfRate, frameDt*1e3, 1e3/P.frameRate, sfWidth, P.numFramesPerBuffer/P.frameRate)

% Full-rate stim waveform (0/1) in the acquisition time base; sample n is at t = (n-1)/fs
nStim = ceil((sfStarts(end) + sfWidth) * fs) + 1;
stimBox = @(onsets) boxcarStim(onsets, P.apis.stim_length_s, fs, nStim);
cleanStim = stimBox(stimOnsetTimestamps);

if opt.showFigs
    figure
    tDAQ = ((1:numel(stimTrace)).' - indTrig) / fs;
    plot(tDAQ, stimTrace); hold on
    yl = [min(stimTrace), max(stimTrace)];
    for k = 1:nSF
        patch([sfStarts(k), sfStarts(k) + sfWidth, sfStarts(k) + sfWidth, sfStarts(k)], yl([2 2 1 1]), 'g', 'FaceAlpha', .3, 'EdgeColor', 'none')
    end
    plot(stimOnsetTimestamps, repmat(yl(2), size(stimOnsetTimestamps)), 'rv', 'MarkerFaceColor', 'r')
    xlim([-1, sfStarts(end) + 5]); xlabel('Time from the Verasonics trigger [s]'); ylabel('Stim trace [V]')
    title('Stim trace, superframes (green) and detected trial onsets (red)')
end

%% 3. Load the PDI of every superframe
d1 = load(fullfile(PDpath, 'fUSdata-1.mat'), 'PDI');
v1 = pdiVolume(d1.PDI, opt.jUse);
PDIallSF = zeros([size(v1), RFcount], 'single'); % (x, y, z, superframe); single to limit memory
okSF = false(RFcount, 1);
for fi = 1:RFcount
    fp = fullfile(PDpath, ['fUSdata-', num2str(fi), '.mat']);
    if ~isfile(fp), PDIallSF(:, :, :, fi) = NaN; continue, end
    d = load(fp, 'PDI');
    PDIallSF(:, :, :, fi) = pdiVolume(d.PDI, opt.jUse);
    okSF(fi) = true;
end
clearvars d d1 v1 fp
if ~all(okSF), warning('%d of %d superframe files are missing (set to NaN)', nnz(~okSF), RFcount), end
assert(size(PDIallSF, 4) == nSF, 'Number of superframes in the data (%d) doesn''t match the RF time tags (%d)', size(PDIallSF, 4), nSF)
PDIA = mean(PDIallSF, 4, 'omitnan');
sz3 = size(PDIA);
vox_mm = [PData.PDelta(2), PData.PDelta(1), PData.PDelta(3)] * P.wl * 1e3; % mm per voxel along dims 1-3 (dim 1 = rows = y, dim 2 = columns = x)

%% 4. Brain mask and ROIs to inject into
brain = PDIA > opt.brainThrFrac * pctl(PDIA(isfinite(PDIA)), 99);
brain = bwareaopen(brain, 500);
presetCenters = zeros(0, 3);
switch opt.roiMode
    case 'preset'
        presetCenters = opt.roiCenterVox;
    case 'click'
        V = PDIA .^ 0.5;
        figure; imagesc(squeeze(max(V, [], 1)).'); axis image; title('Coronal MIP (x vs z): click barrel cortex')
        [xi, zi] = ginput(1);
        figure; imagesc(squeeze(max(V, [], 2)).'); axis image; title('Sagittal MIP (y vs z): click AP position')
        [yi, ~] = ginput(1);
        presetCenters = round([yi, xi, zi]);
end
[roiLabel, roiCenters, blocked] = makeROIs(brain, vox_mm, opt.roiSemi_mm, opt.nROIs, presetCenters);
nROI = size(roiCenters, 1);
roiIdx = arrayfun(@(r) find(roiLabel == r), 1:nROI, 'UniformOutput', false);
roiAny = find(roiLabel > 0);
bgIdx = find(brain & ~blocked); % Background (brain voxels away from all ROIs)
fprintf('Brain mask: %d voxels (%.0f%% of the volume); %d ROIs of ~%d voxels each; %d background voxels\n', ...
    nnz(brain), 100 * nnz(brain) / numel(brain), nROI, round(mean(cellfun(@numel, roiIdx))), numel(bgIdx))

if opt.showFigs
    figure
    subplot(1, 2, 1); imagesc(squeeze(max(PDIA, [], 1)).' .^ 0.5); axis image; colormap hot; hold on
    plot(roiCenters(:, 2), roiCenters(:, 3), 'co', 'MarkerSize', 8, 'LineWidth', 1.5); title('Coronal MIP (dim 2 vs dim 3) + ROI centres')
    subplot(1, 2, 2); imagesc(squeeze(max(PDIA, [], 2)).' .^ 0.5); axis image; colormap hot; hold on
    plot(roiCenters(:, 1), roiCenters(:, 3), 'co', 'MarkerSize', 8, 'LineWidth', 1.5); title('Sagittal MIP (dim 1 vs dim 3) + ROI centres')
end

%% 5. The injected response: the exact regressor the GLM uses (via a dummy 1 x nSF fit)
%   glmActivationMap builds reg = (stim * HRF), averaged over each superframe window; info.X(:, 1) is that regressor
[~, ginfo] = glmActivationMap(rand(1, nSF) + 1, sfStarts, cleanStim, 'Fs', fs, 'SFWidth', sfWidth, ...
    'HRFPeak', opt.injHRFPeak, 'DriftOrder', 0, 'AR', 0);
regInj = ginfo.X(:, 1);
assert(all(isfinite(regInj)), 'Injected regressor has non-finite values (are the superframes outside the stim recording?)')
fprintf('Injected regressor: peak %.3f (a sustained-stim response of amp gives a peak of amp*%.3f)\n', max(regInj), max(regInj))

V0 = reshape(PDIallSF, [], nSF); % (voxels x superframes), no copy
preTime = min(5, onsetInTrial_s);
postTime = P.apis.seq_length_s - onsetInTrial_s;

%% 6. Test A: peri-stimulus averaging recovers the injected ROI-mean response
amp = opt.ampSweep;
Vinj = V0; Vinj(roiAny, :) = V0(roiAny, :) .* (1 + amp * single(regInj.')); % Inject into every ROI voxel
% Expected evoked curve: peri-stimulus average of the noise-free synthetic trace 1 + amp * reg (same settings)
[evExp, tG] = periStimulusAverage(1 + amp * regInj.', sfStarts, stimOnsetTimestamps, 'SFWidth', sfWidth, ...
    'Pre', preTime, 'Post', postTime, 'Normalize', 'percent');
evExp = evExp(:);
A.gain = zeros(nROI, 1); A.corr = zeros(nROI, 1); A.peakRec = zeros(nROI, 1); A.peakOrig = zeros(nROI, 1);
evRecAll = zeros(nROI, numel(tG)); evOrigAll = evRecAll; evSemAll = evRecAll;
for r = 1:nROI
    ts0 = mean(V0(roiIdx{r}, :), 1, 'omitnan');           % ROI-mean PDI, uninjected
    ts1 = mean(Vinj(roiIdx{r}, :), 1, 'omitnan');         % ROI-mean PDI, injected
    [ev0, ~, ~] = periStimulusAverage(ts0, sfStarts, stimOnsetTimestamps, 'SFWidth', sfWidth, 'Pre', preTime, 'Post', postTime, 'Normalize', 'percent');
    [ev1, ~, inf1] = periStimulusAverage(ts1, sfStarts, stimOnsetTimestamps, 'SFWidth', sfWidth, 'Pre', preTime, 'Post', postTime, 'Normalize', 'percent');
    ev0 = ev0(:); ev1 = ev1(:); dEv = ev1 - ev0;          % The difference is the injected part after processing (noise cancels)
    ok = isfinite(dEv) & isfinite(evExp);
    A.gain(r) = (evExp(ok)' * dEv(ok)) / (evExp(ok)' * evExp(ok)); % Least-squares gain of recovered vs expected (1 = perfect)
    c = corrcoef(evExp(ok), dEv(ok)); A.corr(r) = c(1, 2);
    A.peakRec(r) = max(ev1); A.peakOrig(r) = max(ev0);
    evRecAll(r, :) = ev1.'; evOrigAll(r, :) = ev0.'; evSemAll(r, :) = inf1.sem(:).';
end
fprintf('\n== Test A: peri-stimulus average, amp = %.1f%% ==\n', 100 * amp)
fprintf('Recovered/expected gain (injected part): median %.3f [min %.3f, max %.3f]; correlation: median %.4f\n', ...
    median(A.gain), min(A.gain), max(A.gain), median(A.corr))
fprintf('Expected peak %.2f%%; ROI-mean peak of the evoked curve: uninjected median %.2f%%, injected median %.2f%%\n', ...
    max(evExp), median(A.peakOrig), median(A.peakRec))
fprintf('Trials with usable data: %d of %d (peri-stimulus window -%.1f to +%.1f s)\n', numel(inf1.trialsUsed), P.numTrials, preTime, postTime)
clearvars Vinj

if opt.showFigs
    figure; hold on
    errorbar(tG, mean(evOrigAll, 1), mean(evSemAll, 1) / sqrt(nROI), 'k')
    errorbar(tG, mean(evRecAll, 1), mean(evSemAll, 1) / sqrt(nROI), 'r')
    plot(tG, mean(evOrigAll, 1) + evExp.', 'b--', 'LineWidth', 1.5)
    xline(0); xline(P.apis.stim_length_s); xlabel('Time from stim onset [s]'); ylabel('\DeltaPDI [%]')
    legend('uninjected (mean over ROIs)', 'injected', 'uninjected + expected', 'Location', 'best')
    title(sprintf('Test A: ROI-mean evoked response, injected %.1f%%', 100 * amp))
end

%% 7. Test B: voxelwise GLM
%% 7a. Noise-free fidelity: does the GLM return beta = 100 * amp?
stimFor = @(dt, peak) deal(stimBox(stimOnsetTimestamps + dt), peak);
nfRows = 3;
B.fidelity = struct('label', {}, 'betaOverTruth', {});
rowsDef = {};
for peak = opt.hrfPeaks, rowsDef(end+1, :) = {sprintf('HRF peak %g s (stim offset 0)', peak), 0, peak}; end %#ok<SAGROW>
for dt = opt.timeOffsets, rowsDef(end+1, :) = {sprintf('stim offset %+g s (HRF peak %g s)', dt, opt.injHRFPeak), dt, opt.injHRFPeak}; end %#ok<SAGROW>
fprintf('\n== Test B1: noise-free GLM fidelity (recovered beta / injected beta; injection HRF peak %g s, offset 0) ==\n', opt.injHRFPeak)
syn = repmat(1 + amp * regInj.', nfRows, 1) + 1e-6 * randn(nfRows, nSF);
for k = 1:size(rowsDef, 1)
    [stimK, peakK] = stimFor(rowsDef{k, 2}, rowsDef{k, 3});
    rk = glmActivationMap(syn, sfStarts, stimK, 'Fs', fs, 'SFWidth', sfWidth, 'HRFPeak', peakK, 'Normalize', 'percent');
    B.fidelity(k).label = rowsDef{k, 1};
    B.fidelity(k).betaOverTruth = median(rk.beta) / (100 * amp);
    fprintf('  %-40s beta / (100*amp) = %.3f\n', B.fidelity(k).label, B.fidelity(k).betaOverTruth)
end

%% 7b. Real-data noise floor, and detection vs injected amplitude (raw and spatially smoothed)
fprintf('\n== Test B2: GLM on real data with injected responses (hrf peak %g s, offset 0) ==\n', opt.injHRFPeak)
variants = {'raw'; 'smoothed'};
ampList = [0, opt.amps(:).'];
B.det = zeros(numel(ampList), 2); B.betaRatio = B.det; B.bgFP = B.det; B.roiT = B.det; B.roiDet = B.det;
for vi = 1:2
    for ai = 1:numel(ampList)
        a = ampList(ai);
        vol = PDIallSF;
        if a > 0
            Vi = reshape(vol, [], nSF); Vi(roiAny, :) = Vi(roiAny, :) .* (1 + a * single(regInj.')); vol = reshape(Vi, size(PDIallSF)); clear Vi
        end
        if vi == 2, vol = smoothVol(vol, opt.smoothSigma); end
        res = glmActivationMap(vol, sfStarts, cleanStim, 'Fs', fs, 'SFWidth', sfWidth, 'HRFPeak', opt.injHRFPeak, 'Normalize', 'percent');
        rm = roiMeanGLM(vol, roiIdx, sfStarts, cleanStim, fs, sfWidth, opt.injHRFPeak);
        detR = cellfun(@(ix) mean(res.q(ix) < 0.05 & res.beta(ix) > 0, 'omitnan'), roiIdx);
        betaR = cellfun(@(ix) median(res.beta(ix), 'omitnan'), roiIdx);
        B.det(ai, vi) = median(detR);
        if a > 0, B.betaRatio(ai, vi) = median(betaR) / (100 * a); else, B.betaRatio(ai, vi) = NaN; end
        B.bgFP(ai, vi) = mean(res.q(bgIdx) < 0.05 & res.beta(bgIdx) > 0, 'omitnan');
        B.roiT(ai, vi) = median(rm.t); B.roiDet(ai, vi) = mean(rm.p < 0.05 & rm.beta > 0);
        if a == 0
            B.seVox(vi) = median(cellfun(@(ix) median(abs(res.beta(ix) ./ res.t(ix)), 'omitnan'), roiIdx)); % Voxelwise beta standard error [%]
            B.seROI(vi) = median(abs(rm.beta ./ rm.t));                                                   % ROI-mean beta standard error [%]
            B.realRoiBeta(vi) = median(rm.beta);
        end
        if a == opt.ampSweep && vi == 2, B.zSmoothed = res.z; end
        if a == opt.ampSweep && vi == 1, B.zRaw = res.z; end
        fprintf('  %-8s amp %5.1f%%: voxels detected (q<0.05) %5.1f%% (median over ROIs), median beta/truth %.2f, ROI-mean t %.1f (p<.05 in %d/%d ROIs), background false-positive rate %.2f%%\n', ...
            variants{vi}, 100 * a, 100 * B.det(ai, vi), B.betaRatio(ai, vi), B.roiT(ai, vi), round(B.roiDet(ai, vi) * nROI), nROI, 100 * B.bgFP(ai, vi))
    end
end
fprintf('\nNoise floor of the REAL data (GLM beta standard error, in %% PDI):\n')
fprintf('  raw:      voxelwise %.3f%%  -> detectable at q<0.05 (z~4.5) only above ~%.2f%%;  ROI-mean %.3f%% -> 3-sigma at ~%.2f%%\n', B.seVox(1), 4.5 * B.seVox(1), B.seROI(1), 3 * B.seROI(1))
fprintf('  smoothed: voxelwise %.3f%%  -> detectable at q<0.05 (z~4.5) only above ~%.2f%%;  ROI-mean %.3f%% -> 3-sigma at ~%.2f%%\n', B.seVox(2), 4.5 * B.seVox(2), B.seROI(2), 3 * B.seROI(2))
fprintf('  (compare with the response you expect, and with the ~0.3%% ROI response seen so far)\n')

%% 7c. Sensitivity to the stim timing offset and to the analysis HRF peak (smoothed data)
amp = opt.ampSweep;
Vi = reshape(PDIallSF, [], nSF); Vi(roiAny, :) = Vi(roiAny, :) .* (1 + amp * single(regInj.')); volI = smoothVol(reshape(Vi, size(PDIallSF)), opt.smoothSigma); clear Vi
fprintf('\n== Test B3: analysis timing offset / HRF mismatch (smoothed data, injected %.1f%%, injection HRF peak %g s) ==\n', 100 * amp, opt.injHRFPeak)
B.sweep = struct('label', {}, 'det', {}, 'roiT', {}, 'betaRatio', {});
for k = 1:size(rowsDef, 1)
    [stimK, peakK] = stimFor(rowsDef{k, 2}, rowsDef{k, 3});
    res = glmActivationMap(volI, sfStarts, stimK, 'Fs', fs, 'SFWidth', sfWidth, 'HRFPeak', peakK, 'Normalize', 'percent');
    rm = roiMeanGLM(volI, roiIdx, sfStarts, stimK, fs, sfWidth, peakK);
    B.sweep(k).label = rowsDef{k, 1};
    B.sweep(k).det = median(cellfun(@(ix) mean(res.q(ix) < 0.05 & res.beta(ix) > 0, 'omitnan'), roiIdx));
    B.sweep(k).roiT = median(rm.t);
    B.sweep(k).betaRatio = median(cellfun(@(ix) median(res.beta(ix), 'omitnan'), roiIdx)) / (100 * amp);
    fprintf('  %-40s voxels detected %5.1f%%, ROI-mean t %5.1f, beta/truth %.2f\n', B.sweep(k).label, 100 * B.sweep(k).det, B.sweep(k).roiT, B.sweep(k).betaRatio)
end
clearvars volI res rm

if opt.showFigs
    figure
    subplot(1, 3, 1); hold on
    plot(100 * ampList, 100 * B.det(:, 1), '-o'); plot(100 * ampList, 100 * B.det(:, 2), '-s')
    xlabel('Injected amplitude [%]'); ylabel('ROI voxels detected [%] (q<0.05)'); legend('raw', 'smoothed', 'Location', 'best'); title('B2: voxelwise detection')
    subplot(1, 3, 2); hold on
    plot(100 * ampList, B.roiT(:, 1), '-o'); plot(100 * ampList, B.roiT(:, 2), '-s'); yline(2, 'k--')
    xlabel('Injected amplitude [%]'); ylabel('ROI-mean GLM t'); title('B2: ROI-mean statistic')
    subplot(1, 3, 3)
    bar(categorical({B.sweep.label}, {B.sweep.label}), 100 * [B.sweep.det]); ylabel('ROI voxels detected [%]'); title(sprintf('B3: mismatch (%.0f%% injected)', 100 * amp))
    figure
    subplot(1, 2, 1); imagesc(squeeze(max(B.zRaw, [], 1)).'); axis image; colorbar; title(sprintf('z (raw), %.0f%% injected: max over dim 1', 100 * amp))
    subplot(1, 2, 2); imagesc(squeeze(max(B.zSmoothed, [], 1)).'); axis image; colorbar; title('z (smoothed)')
end

%% 8. Test C: legacy resampling path (interp1 + 'extrap' -> trial average -> correlation)
C = struct();
if opt.runLegacy
    amp = opt.ampSweep; rr = 10;
    nBg = min(opt.legacyNBg, numel(bgIdx));
    bgSub = bgIdx(randperm(numel(bgIdx), nBg));
    sub = [roiAny(:); bgSub(:)]; isRoi = [true(numel(roiAny), 1); false(nBg, 1)];
    Vs0 = double(V0(sub, :));
    Vs1 = Vs0; Vs1(isRoi, :) = Vs0(isRoi, :) .* (1 + amp * regInj.');
    fprintf('\n== Test C: legacy resampling path (rr = %d Hz), injected %.1f%% ==\n', rr, 100 * amp)
    fprintf('  (r = Pearson correlation of the trial-averaged ratio with the stim; "box" = undelayed boxcar as in evokedResponse3D.m, "HRF" = boxcar * gamma HRF, peak %g s)\n', opt.injHRFPeak)
    for fld = {'extrap', 'no_extrap'}
        f = fld{1}; ex = strcmp(f, 'extrap');
        o0 = legacyPath(Vs0, sfCenters, P.apis.seq_length_s, P.numTrials, onsetInTrial_s, P.apis.stim_length_s, rr, ex, opt.injHRFPeak);
        o1 = legacyPath(Vs1, sfCenters, P.apis.seq_length_s, P.numTrials, onsetInTrial_s, P.apis.stim_length_s, rr, ex, opt.injHRFPeak);
        oe = legacyPath(1 + amp * regInj.', sfCenters, P.apis.seq_length_s, P.numTrials, onsetInTrial_s, P.apis.stim_length_s, rr, ex, opt.injHRFPeak); % noise-free
        C.(f).TAroiInj = mean(o1.TA(:, isRoi), 2);
        C.(f).TAroiOrig = mean(o0.TA(:, isRoi), 2);
        C.(f).nExtrap = o1.nExtrap;
        fprintf('  %s: %d of %d resampled points extrapolated (%.2f s)\n', f, o1.nExtrap, numel(o1.tq), o1.nExtrap / rr)
        for ti = {'r', 'rT'}
            t = ti{1}; lab = 'box'; if strcmp(t, 'rT'), lab = 'HRF'; end
            thr = pctl(o0.(t)(~isRoi), 99);   % Empirical null: 99th percentile of r in the (uninjected) background
            C.(f).(t).ceiling = oe.(t)(1);    % r of the noise-free injected response: the best any voxel could do
            C.(f).(t).rRoiOrig = median(o0.(t)(isRoi));
            C.(f).(t).rRoiInj = median(o1.(t)(isRoi));
            C.(f).(t).thr = thr;
            C.(f).(t).detOrig = mean(o0.(t)(isRoi) > thr);
            C.(f).(t).detInj = mean(o1.(t)(isRoi) > thr);
            fprintf('    %-3s: noise-free ceiling r = %.2f; median ROI r %.3f (uninjected) -> %.3f (injected); ROI voxels above the background 99th percentile (r > %.2f): %.1f%% -> %.1f%%\n', ...
                lab, C.(f).(t).ceiling, C.(f).(t).rRoiOrig, C.(f).(t).rRoiInj, thr, 100 * C.(f).(t).detOrig, 100 * C.(f).(t).detInj)
        end
    end
    clearvars o0 o1 oe Vs0 Vs1
    dTA = max(abs(C.extrap.TAroiInj - C.no_extrap.TAroiInj));
    fprintf('  Max difference of the ROI-mean trial-average curve between extrap and no extrap: %.2e (ratio units; %.2f%% of baseline)\n', dTA, 100 * dTA)
end

%% 9. Summary
fprintf('\n================ SUMMARY ================\n')
fprintf('A  peri-stimulus average: recovered/expected gain %.2f (median over ROIs), correlation %.3f   -> %s\n', median(A.gain), median(A.corr), passStr(abs(median(A.gain) - 1) < 0.1 && median(A.corr) > 0.95))
iM = find(opt.hrfPeaks == opt.injHRFPeak, 1);
if isempty(iM), fidM = NaN; else, fidM = B.fidelity(iM).betaOverTruth; end
fprintf('B1 GLM noise-free fidelity at the matched HRF/timing: beta/truth %.2f                   -> %s\n', fidM, passStr(abs(fidM - 1) < 0.05))
mdaIdx = find(B.det(2:end, 2) >= 0.5, 1);
if isempty(mdaIdx), mdaStr = sprintf('> %.1f%%', 100 * max(opt.amps)); else, mdaStr = sprintf('%.1f%%', 100 * opt.amps(mdaIdx)); end
fprintf('B2 smallest injected amplitude detected in >=50%% of ROI voxels (smoothed, q<0.05): %s\n', mdaStr)
mdaRoi = find(B.roiT(2:end, 1) > 3, 1);
if isempty(mdaRoi), mdaRoiStr = sprintf('> %.1f%%', 100 * max(opt.amps)); else, mdaRoiStr = sprintf('%.1f%%', 100 * opt.amps(mdaRoi)); end
fprintf('   smallest injected amplitude with ROI-mean GLM t > 3 (raw data): %s\n', mdaRoiStr)
fprintf('   noise floor (beta SE): voxelwise %.3f%% (raw) / %.3f%% (smoothed); ROI-mean %.3f%% / %.3f%%\n', B.seVox(1), B.seVox(2), B.seROI(1), B.seROI(2))
fprintf('If A and B1 pass, the processing recovers a known response; compare the noise floor with the response you are looking for.\n')
fprintf('If they fail, the chain (timing, normalization, regressor) is at fault: inspect the Test B1 rows and the timing printout above.\n')

PC = struct('opt', opt, 'timing', struct('sfStarts', sfStarts, 'sfWidth', sfWidth, 'stimOnsetTimestamps', stimOnsetTimestamps), ...
    'roiCenters', roiCenters, 'tG', tG, 'A', A, 'B', rmfield(B, intersect(fieldnames(B), {'zRaw', 'zSmoothed'})), 'C', C);
if ~isempty(opt.savePath)
    save(fullfile(opt.savePath, 'positiveControl_results.mat'), 'PC')
end

%% Local functions
function p = passStr(tf)
    if tf, p = 'PASS'; else, p = 'CHECK'; end
end

% Look for a file in a list of folders; ask the user if it isn't in any of them
function f = findFile(folders, name)
    for k = 1:numel(folders)
        f = [folders{k}, name];
        if isfile(f), return, end
    end
    [fn, fp] = uigetfile('*.mat', ['Select ', name], folders{1});
    f = [fp, fn];
end

function v = pdiVolume(PDI, jUse)
    if iscell(PDI), v = PDI{jUse}; else, v = PDI; end
end

function p = pctl(x, pc)
    x = sort(x(:)); n = numel(x);
    p = x(min(n, max(1, round(pc / 100 * n))));
end

function s = boxcarStim(onsets, len, fs, nStim)
    s = zeros(nStim, 1);
    for k = 1:numel(onsets)
        i0 = round(onsets(k) * fs) + 1; i1 = round((onsets(k) + len) * fs);
        s(max(i0, 1):min(i1, nStim)) = 1;
    end
end

function vol = smoothVol(vol, sigma)
    if sigma <= 0, return, end
    for t = 1:size(vol, 4)
        vol(:, :, :, t) = imgaussfilt3(vol(:, :, :, t), sigma);
    end
end

function rm = roiMeanGLM(vol, roiIdx, sfStarts, stim, fs, sfWidth, hrfPeak)
    V = reshape(vol, [], size(vol, 4));
    S = zeros(numel(roiIdx), size(vol, 4));
    for r = 1:numel(roiIdx), S(r, :) = mean(V(roiIdx{r}, :), 1, 'omitnan'); end
    rm = glmActivationMap(S, sfStarts, stim, 'Fs', fs, 'SFWidth', sfWidth, 'HRFPeak', hrfPeak, 'Normalize', 'percent');
end

% Place non-overlapping ellipsoidal ROIs in the brain mask (preset centres first, then random ones)
function [roiLabel, centers, blocked] = makeROIs(brain, vox_mm, semi_mm, nRandom, presetCenters)
    sz = size(brain);
    roiLabel = zeros(sz, 'uint8'); blocked = false(sz); centers = zeros(0, 3);
    hw = ceil(semi_mm ./ vox_mm);
    [di, dj, dk] = ndgrid(-hw(1):hw(1), -hw(2):hw(2), -hw(3):hw(3));
    ell = ((di * vox_mm(1) / semi_mm(1)).^2 + (dj * vox_mm(2) / semi_mm(2)).^2 + (dk * vox_mm(3) / semi_mm(3)).^2) <= 1;
    cand = find(brain); cand = cand(randperm(numel(cand)));
    queue = [presetCenters; zeros(0, 3)];
    nPreset = size(presetCenters, 1);
    ci = 0; nPlaced = 0;
    while nPlaced < nPreset + nRandom
        if nPlaced < nPreset
            c = round(queue(nPlaced + 1, :)); isPre = true;
        else
            ci = ci + 1;
            if ci > numel(cand), warning('Could only place %d of %d ROIs', nPlaced, nPreset + nRandom), break, end
            [a, b, d] = ind2sub(sz, cand(ci)); c = [a, b, d]; isPre = false;
        end
        lo = c - hw; hi = c + hw;
        if any(lo < 1) || any(hi > sz)
            if isPre, nPlaced = nPlaced + 1; end
            continue
        end
        if ~isPre
            if blocked(c(1), c(2), c(3)), continue, end
            cov = brain(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3));
            if mean(cov(ell)) < 0.95, continue, end
        end
        nPlaced = nPlaced + 1;
        centers(end + 1, :) = c; %#ok<AGROW>
        sub = roiLabel(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3));
        sub(ell) = size(centers, 1);
        roiLabel(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)) = sub;
        l2 = max(c - 2 * hw, 1); h2 = min(c + 2 * hw, sz);
        blocked(l2(1):h2(1), l2(2):h2(2), l2(3):h2(3)) = true; % Keep a gap of >= 1 ROI half-width between ROIs
    end
end

% Legacy path from evokedResponse3D.m: interp1 onto a regular grid, per-trial baseline ratio, trial average, correlation with the stim.
%   Vs: (voxels x superframes)
%   out.r: correlation with the boxcar stim (what evokedResponse3D.m does); out.rT: correlation with an HRF-convolved stim template
function out = legacyPath(Vs, sfCenters, seqLen, nTrials, onsetInTrial, stimLen, rr, useExtrap, hrfPeak)
    tq = (0:1/rr:(nTrials * seqLen - 1/rr)).';
    if useExtrap
        Y = interp1(sfCenters, Vs.', tq, 'linear', 'extrap');
    else
        Y = interp1(sfCenters, Vs.', tq, 'linear');
    end
    nW = seqLen * rr; nb = round(onsetInTrial * rr);
    S = zeros(nW, size(Vs, 1)); N = S;
    for k = 1:nTrials
        W = Y((k - 1) * nW + 1 : k * nW, :);
        R = W ./ mean(W(1:nb, :), 1, 'omitnan');
        ok = ~isnan(R); R(~ok) = 0;
        S = S + R; N = N + ok;
    end
    TA = S ./ N; TA(N == 0) = 0;
    stimrs = zeros(nW, 1); stimrs(nb + 1 : nb + round(stimLen * rr)) = 1;
    Z = TA - mean(TA, 1); s = stimrs - mean(stimrs);
    out.r = (s.' * Z) ./ (sqrt(sum(Z.^2, 1)) * norm(s));
    out.r(~isfinite(out.r)) = 0;
    th = (0:1/rr:30).'; h = th.^3 .* exp(-th / (hrfPeak / 3)); h = h / sum(h); % Gamma HRF (shape 4), as in glmActivationMap
    tmpl = conv(stimrs, h); tmpl = tmpl(1:nW); tmpl = tmpl - mean(tmpl);
    out.rT = (tmpl.' * Z) ./ (sqrt(sum(Z.^2, 1)) * norm(tmpl));
    out.rT(~isfinite(out.rT)) = 0;
    out.TA = TA; out.tq = tq;
    out.nExtrap = nnz(tq < sfCenters(1) | tq > sfCenters(end));
end
