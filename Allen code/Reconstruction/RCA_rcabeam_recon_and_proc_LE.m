%% Description

% Block ("long ensemble") version of RCA_rcabeam_recon_and_proc.m.
%
% Reconstruct the raw RF data of a row-column array (RC15gV) acquisition folder with the GPU rcabeam library
% (RCA_DAS_rcabeam) and process it as in IQ2g1T_3D_loop.m, but not one RF file ("superframe") at a time. The selected RF files
% are cut into contiguous, non-overlapping BLOCKS of Opts.blockSize superframes, and for each block:
%   1. every superframe of the block is beamformed on the GPU, and its IQ gathered to RAM (z cropped to zstart:zend)
%   2. the block's IQ is joined in time into one long ensemble, and the SVD clutter filter (and the high pass filter) is run on
%      the WHOLE long ensemble (the temporal covariance is accumulated over the whole block, and one filter matrix is applied)
%   3. the filtered IQ is re-sliced in time into new "superframes" of sliceLen frames (non-overlapping, any length)
%   4. each new superframe gives PDI and g1, as in IQ2g1T_3D_nofreqsep_fast (mean power, and g1 for nTau lags)
% Nothing is carried from one block to the next (no overlap), so blocks can be processed, skipped or redone independently.

% Input folder: params.mat (structure P) and either
%           RF-<maxAngle>-<na>-<frameRate>-<numFrames>-1-<file #>.mat   (RcvData)         Verasonics .mat saves, or
%           <P.fsFileName><file #, 5 digits>.rf                          (binary, int16)   Verasonics FastSave
%         (see Opts.fileFormat)
% Output: one file per block in the processing save folder, fUSdata-block-<first file #>-<last file #>.mat, with
%           PDI              [y, x, z, slice]        power Doppler image of each re-sliced superframe
%           g1               [y, x, z, lag, slice]   temporal field autocorrelation of each re-sliced superframe
%           fileNumbers      the RF file numbers in the block
%           nfBlock, sliceLen, nSlices, sliceStartFrame   (first frame of each slice, counted from the start of the block)
%         and PData.mat, rcabeam_recon_params.mat (the pixel grid in m, and the options used) and fUS_proc_params.mat in the
%         save folder. The IQ is not saved

% Notes:
%           - Setup is the same as RCA_DAS_rcabeam: MATLAB R2026b+ with the Parallel Computing Toolbox, and the python.exe
%             of an environment with cupy and rcabeam in Opts.pythonExe or the RCABEAM_PYTHON environment variable
%           - Blocks are made inside runs of consecutive file numbers, so a block never spans a gap in the selected files. Files
%             left over at the end of a run (fewer than blockSize) are dropped unless Opts.keepPartialBlock = true, which
%             processes them as a shorter block. All files of a block must have the same number of frames, else the block is skipped
%           - blockSize = 1 and sliceLen = frames per file is the same as RCA_rcabeam_recon_and_proc.m (one file per ensemble)
%           - RAM: the whole block's IQ (cropped to zstart:zend) is held in complex single, the whole IQ volume x frames x 8 bytes
%             x blockSize. It is filtered in place, in chunks of voxels that go through the GPU
%           - The SVD thresholds count singular values of the WHOLE block, so they do not mean the same as for a single file.
%             The upper threshold can be 'end' (all frames of the block). With more frames per ensemble, the lower one that
%             removes the same clutter is usually larger
%           - The high pass filter is applied to the whole block together with the SVD filter (Opts.hpfOnBlock = true, default),
%             so only the start of the block has a filter transient, not the start of every slice. false filters each slice
%             separately instead (the transient is in the start of every slice, as when each file is processed alone)
%           - Frames at the end of a block that do not fill a whole slice are not used (printed)
%           - The pixel grid, the channel remapping, the .rf/.mat reading and the IQ conventions are the same as
%             RCA_rcabeam_recon.m. Volumes are [y, x, z] as in the Verasonics IQ files
%           - The slices have no overlap. Slices of one block share the SVD filter, so they are not independent of each other

% Inputs (asked with dialogs):
%   raw data folder, save folder for the processed data, then a dialog with the file numbers, voxel spacing, output and
%   stacking (same as RCA_rcabeam_recon.m), and a dialog with the processing parameters:
%     SVD lower / upper bound      kept singular values of the whole block (upper: a number, or 'end')
%     z start / z end              z index range kept (as in RCA_rcabeam_recon_and_proc.m)
%     High pass filter cutoff      [Hz]
%     Time lag window              [ms], g1 is calculated for that many lags
%     Block size                   # of superframes (RF files) in one block
%     Re-sliced superframe length  # of frames in each new superframe
% Opts fields that are not asked for (set them in the script below): fileFormat, overwrite, keepPartialBlock, hpfOnBlock,
%   das (options for RCA_DAS_rcabeam: frameChunk, fNumber, pythonExe, backend, ...)

%% Add the rcabeam/processing folders (and subfolders) to path
% reconDir = fileparts(mfilename('fullpath'));
reconDir = cd;
if isempty(reconDir), reconDir = pwd; end
reconDir_split = split(string(reconDir), filesep);
iGitHub = find(reconDir_split == "GitHub", 1);
if ~isempty(iGitHub)
    rcabeamCodePath = fullfile(strjoin(reconDir_split(1:iGitHub), filesep), 'rcabeam');
    if isfolder(rcabeamCodePath), addpath(genpath(rcabeamCodePath)), end
end
addpath(genpath(fullfile(fileparts(reconDir), 'Processing')))
addpath(reconDir) % after the Processing folder, which has its own copies of some of the functions

datapath = uigetdir(pwd, 'Select the raw data path');
if isequal(datapath, 0)
    error('No raw data path selected.')
end
PDsavepath = uigetdir(datapath, 'Select the save path for processed data');
if isequal(PDsavepath, 0)
    error('No save path selected.')
end

Opts = struct();
showGUI = true;
Opts = setDefault(Opts, 'fileRange', [0, Inf]);
Opts = setDefault(Opts, 'fileNumbers', []);
Opts = setDefault(Opts, 'pixelSpacingWl', [0.5, 0.5, 0.5]);
Opts = setDefault(Opts, 'output', 'OPW');
Opts = setDefault(Opts, 'unstackedFrames', true);
Opts = setDefault(Opts, 'fileFormat', 'auto');
Opts = setDefault(Opts, 'overwrite', false);
Opts = setDefault(Opts, 'keepPartialBlock', false);
Opts = setDefault(Opts, 'hpfOnBlock', true);
Opts = setDefault(Opts, 'das', struct());
Opts.das.gpuOutput = true; % Keep the IQ on the GPU until each chunk of frames is gathered to RAM (false if the GPU memory is short)

%% Load the acquisition parameters
S = load(fullfile(datapath, 'params.mat'), 'P');
P = S.P;
if ~isfield(P, 'numFramesPerBuffer') && isfield(P, 'numFramesPerBuff') % older P structures
    P.numFramesPerBuffer = P.numFramesPerBuff;
end
if ~isequal(P.Trans.name, 'RC15gV')
    error('This is for the RC15gV row-column array, not %s', P.Trans.name)
end

c = P.Resource.Parameters.speedOfSound;
wl = c / (P.Trans.frequency * 1e6);

%% Prompt for processing params
procParamsPrompt = {'SVD lower bound (of the whole block)', 'SVD upper bound (number, or end = all frames of the block)', 'z start', 'z end', ...
    'High pass filter cutoff [Hz]', 'Time lag window [ms]', 'Block size [# of superframes (RF files)]', 'Re-sliced superframe length [# of frames]'};
procParamsDefaults = {'50', 'end', '1', '110', '25', '20', '2', num2str(P.numFramesPerBuffer)};
procParamsUserInput = inputdlg(procParamsPrompt, 'Input Parameters', 1, procParamsDefaults);
if isempty(procParamsUserInput)
    fprintf('Cancelled.\n')
    return
end

sv_threshold_lower = str2double(procParamsUserInput{1});
svUpperText = lower(strtrim(procParamsUserInput{2}));
if isempty(svUpperText) || ismember(svUpperText, {'end', 'inf'})
    sv_threshold_upper = Inf; % all frames of the block, set per block
else
    sv_threshold_upper = str2double(svUpperText);
end
zstart = str2double(procParamsUserInput{3});
zend = str2double(procParamsUserInput{4});
tau_window_s = str2double(procParamsUserInput{6}) / 1e3;
blockSize = str2double(procParamsUserInput{7});
sliceLen = str2double(procParamsUserInput{8});
if any(isnan([sv_threshold_lower, sv_threshold_upper, zstart, zend, tau_window_s, blockSize, sliceLen])) || ...
        blockSize < 1 || sliceLen < 2 || blockSize ~= round(blockSize) || sliceLen ~= round(sliceLen)
    error('Can''t use these processing parameters. The block size and the re-sliced superframe length must be whole numbers')
end

% Set up the High Pass Filter (parameters from the 2020 vUS paper)
HPF.fc = str2double(procParamsUserInput{5}); % Cutoff frequency [Hz]
% 25 Hz corresponds to 1 mm/s
HPF.order = 4; % Butterworth filter order
frameRate = frameRateFromTimeTags(P, datapath); % like checkFrameRate(P), without its uiopen dialog
HPF.fs = frameRate; % Sampling frequency [Hz]
[HPF.b, HPF.a] = butter(HPF.order, HPF.fc/(HPF.fs/2), 'high');

nTau = ceil(tau_window_s * frameRate); % # of time lags to consider; empirically set by assuming all g1 for voxels containing actual flow decay within 10 ms
if nTau >= sliceLen
    error('The time lag window (%d lags) must be shorter than the re-sliced superframe (%d frames)', nTau, sliceLen)
end
tau = (0:nTau - 1)' ./ frameRate; % Time lag vector [s]
clearvars procParamsPrompt procParamsDefaults procParamsUserInput svUpperText

%% Find the RF files, in order of file number
[f, fileNum, fmt] = findRFFiles(datapath, P, Opts.fileFormat);
fprintf('%d RF files (.%s) in %s\n', numel(f), fmt, datapath)

%% Options from the dialog
if showGUI
    Opts = reconDialog(Opts, fileNum, P, wl, datapath, PDsavepath, fmt);
    if isempty(Opts)
        fprintf('Reconstruction cancelled.\n')
        return
    end
end

%% Choose the files
if ~isempty(Opts.fileNumbers)
    missing = setdiff(Opts.fileNumbers, fileNum);
    if ~isempty(missing)
        warning('No RF file for file number(s) %s, skipping them', mat2str(missing(:)'))
    end
    sel = ismember(fileNum, Opts.fileNumbers);
    what = sprintf('the file numbers %s', mat2str(Opts.fileNumbers(:)'));
else
    sel = fileNum >= Opts.fileRange(1) & fileNum <= Opts.fileRange(2);
    what = sprintf('the file range [%g, %g]', Opts.fileRange(1), Opts.fileRange(2));
end
f = f(sel);
fileNum = fileNum(sel);
if isempty(f)
    error('No RF files in %s in %s', datapath, what)
end

outputMode = lower(Opts.output);
switch outputMode
    case 'opw'
        Opts.das.sumOnly = true;
    case 'rccr'
        error('Downstream processing is currently not built for RCCR')
    case 'allangles'
        error('rcabeam has no per-angle output. Use RCA_ffdas_recon with Opts.output = ''AllAngles''.')
    otherwise
        error('Unknown output type "%s". Use OPW or RCCR.', Opts.output)
end

useGPU = ~isfield(Opts.das, 'backend') || strcmpi(Opts.das.backend, 'rcabeam');
if useGPU && gpuDeviceCount == 0
    error('No supported GPU found. The rcabeam backend needs an NVIDIA GPU.')
end
if gpuDeviceCount == 0
    error('No supported GPU found. The block filtering also runs on the GPU.')
end

pix = Opts.pixelSpacingWl * wl; % [x, y, z] in m

if ~exist(PDsavepath, 'dir')
    mkdir(PDsavepath)
end

%% Plan the blocks: contiguous, non-overlapping groups of blockSize files, inside runs of consecutive file numbers
runStart = [1; find(diff(fileNum(:)) > 1) + 1];
runEnd = [runStart(2:end) - 1; numel(fileNum)];
if numel(runStart) > 1
    fprintf('The selected file numbers have %d runs of consecutive files, blocks are made inside each run\n', numel(runStart))
end
blocks = {}; % each: indices into f / fileNum
for r = 1:numel(runStart)
    idx = runStart(r):runEnd(r);
    nFull = floor(numel(idx) / blockSize);
    for b = 1:nFull
        blocks{end+1} = idx((b-1) * blockSize + (1:blockSize)); %#ok<SAGROW>
    end
    leftover = idx(nFull * blockSize + 1 : end);
    if ~isempty(leftover)
        if Opts.keepPartialBlock
            blocks{end+1} = leftover; %#ok<SAGROW>
            fprintf('Processing the %d leftover file(s) %s as a shorter block\n', numel(leftover), mat2str(reshape(fileNum(leftover), 1, [])))
        else
            fprintf('Not processing the %d leftover file(s) %s (fewer than the block size %d, Opts.keepPartialBlock = true to use them)\n', ...
                numel(leftover), mat2str(reshape(fileNum(leftover), 1, [])), blockSize)
        end
    end
end
if isempty(blocks)
    error('No complete block of %d files in %s. Use a smaller block size or Opts.keepPartialBlock = true', blockSize, what)
end
fprintf('%d block(s) of up to %d superframes, re-sliced into superframes of %d frames\n', numel(blocks), blockSize, sliceLen)

prm = struct('zstart', zstart, 'zend', zend, 'svLower', sv_threshold_lower, 'svUpper', sv_threshold_upper, 'HPF', HPF, ...
    'hpfOnBlock', Opts.hpfOnBlock, 'nTau', nTau, 'sliceLen', sliceLen, 'gpuChunkBytes', 0.5e9);

%% Reconstruct and process, block by block
cleanupGPU = onCleanup(@() RCA_DAS_rcabeam('clear')); % release the GPU workspace, also after an error or Ctrl+C
PDfiles = {};
grid = [];
for b = 1:numel(blocks)
% for b = 1
    idx = blocks{b};
    fnb = fileNum(idx);
    outFile = fullfile(PDsavepath, sprintf('fUSdata-block-%d-%d.mat', fnb(1), fnb(end)));
    if exist(outFile, 'file') == 2 && ~Opts.overwrite
        fprintf('Block %d/%d (files %d to %d): %s already exists, skipping.\n', b, numel(blocks), fnb(1), fnb(end), outFile)
        PDfiles{end+1} = outFile; %#ok<SAGROW>
        continue
    end

    fprintf('Block %d/%d: files %s\n', b, numel(blocks), mat2str(fnb(:)'))
    tBlock = tic;
    [PDI, g1, info] = processBlock(f(idx), fnb, fmt, P, pix, Opts, prm);
    if ~info.ok
        warning('Block %d (files %s) was skipped: %s', b, mat2str(fnb(:)'), info.reason)
        continue
    end
    grid = info.grid;
    voxelRange = {1:info.nx, 1:info.ny, info.zrange}; % x, y, z indices kept (of the whole reconstructed volume)
    xp = info.nx; yp = info.ny; zp = numel(info.zrange);

    fileNumbers = fnb; %#ok<NASGU>
    nfBlock = info.nfBlock; %#ok<NASGU>
    nSlices = info.nSlices; %#ok<NASGU>
    sliceStartFrame = info.sliceStartFrame; %#ok<NASGU>
    save(outFile, 'PDI', 'g1', 'fileNumbers', 'nfBlock', 'sliceLen', 'nSlices', 'sliceStartFrame', '-v7.3', '-nocompression')
    PDfiles{end+1} = outFile; %#ok<SAGROW>
    fprintf('Block %d/%d done and saved as %s (%.1f s)\n', b, numel(blocks), outFile, toc(tBlock))

    % clearvars PDI g1 info
end

if ~isempty(grid) % empty if every block was skipped
    % PData like the Verasonics recon's: size is [# rows (y), # columns (x), # sections (z)], PDelta and Origin in wavelengths
    if ~isfield(Opts.das, 'grid') || strcmpi(Opts.das.grid, 'verasonics')
        PData.Coord = 'rectangular';
        PData.PDelta = Opts.pixelSpacingWl;
        PData.Size = [numel(grid.y), numel(grid.x), numel(grid.z)];
        PData.Origin = [grid.x(1), grid.y(1), grid.z(1)] / wl;
        save(fullfile(PDsavepath, 'PData.mat'), 'PData')
    end

    ReconOpts = Opts; %#ok<NASGU>
    save(fullfile(PDsavepath, 'rcabeam_recon_params.mat'), 'grid', 'ReconOpts') % grid (m) is the same for every file
    save(fullfile(PDsavepath, 'fUS_proc_params.mat'), 'sv_threshold_lower', 'sv_threshold_upper', 'HPF', 'voxelRange', 'zstart', 'zend', ...
        'tau', 'nTau', 'tau_window_s', 'xp', 'yp', 'zp', 'frameRate', 'blockSize', 'sliceLen')
end


%% Beamform the files of one block, filter the whole block as one ensemble, re-slice it and get PDI and g1 of each slice.
% fb: file struct array (name, folder) of the block, fnb: their file numbers. prm: processing parameters (see above).
% Returns PDI [y, x, z, slice], g1 [y, x, z, lag, slice] and info (ok, reason, grid, nx, ny, zrange, nfBlock, nSlices, sliceStartFrame).
% The block's IQ is a [voxel, frame] matrix IQblock (y fastest, then x, then z), which is only modified by indexed assignments
% inside this function, so MATLAB does not copy it
function [PDI, g1, info] = processBlock(fb, fnb, fmt, P, pix, Opts, prm)

    PDI = []; g1 = [];
    info = struct('ok', false, 'reason', '');
    nFiles = numel(fb);
    gatherChunk = 50; % frames gathered from the GPU at a time

    %% 1. Beamform each superframe on the GPU and gather it to RAM
    IQblock = [];
    nfFile = 0; c0 = 0;
    for k = 1:nFiles
        tFile = tic;
        if strcmp(fmt, 'rf')
            RcvData = loadBinaryRF(fullfile(fb(k).folder, fb(k).name), P); % stacked [rowsPerFrame, colsPerFrame], int16
        else
            S = load(fullfile(fb(k).folder, fb(k).name), 'RcvData');
            RcvData = S.RcvData;
            clear S
        end
        if Opts.unstackedFrames
            r = unstackFrames(RcvData, P); % (# samples per frame, # channels, # frames)
        else
            r = RcvData;
            if iscell(r), r = r{1}; end
        end
        clear RcvData
        nfk = size(r, 3);

        if k == 1
            nfFile = nfk;
        elseif nfk ~= nfFile
            RCA_DAS_rcabeam('clear');
            info.reason = sprintf('file %d has %d frames, but file %d has %d (a file that is not complete?)', fnb(k), nfk, fnb(1), nfFile);
            return
        end

        % RCA_DAS_rcabeam takes the frame count from P.numSubFrames when that exists, so make sure it matches this file
        Pk = P;
        if isfield(Pk, 'numSubFrames')
            Pk.numSubFrames = nfk;
        end

        [IQ_CR, ~, grid] = RCA_DAS_rcabeam(r, Pk, pix(2), pix(1), pix(3), Opts.das); % [x, y, z, frame], on the GPU if Opts.das.gpuOutput
        clear r
        zr = prm.zstart : min(prm.zend, size(IQ_CR, 3));

        if k == 1
            nx = size(IQ_CR, 1); ny = size(IQ_CR, 2); nz = numel(zr);
            nvox = ny * nx * nz;
            IQblock = zeros(nvox, nfFile * nFiles, 'like', complex(single(0))); % complex single, voxel x frame
            fprintf('  The block''s IQ is %d x %d x %d voxels x %d frames, %.1f GB in RAM\n', ny, nx, nz, nfFile * nFiles, nvox * nfFile * nFiles * 8 / 1e9)
        end

        for a = 1:gatherChunk:nfk
            bb = min(a + gatherChunk - 1, nfk);
            tmp = gather(IQ_CR(:, :, zr, a:bb));
            IQblock(:, c0 + (a:bb)) = reshape(permute(tmp, [2, 1, 3, 4]), nvox, []); % Verasonics order: [y, x, z, frame]
        end
        c0 = c0 + nfk;
        clear IQ_CR tmp
        fprintf('  RF file %d beamformed and gathered to RAM (%.1f s)\n', fnb(k), toc(tFile))
    end
    RCA_DAS_rcabeam('clear'); % free the GPU memory for the filtering
    nfB = size(IQblock, 2);

    %% 2. SVD clutter filter (and HPF) on the whole long ensemble
    tFilt = tic;
    svLo = prm.svLower;
    svHi = min(prm.svUpper, nfB);
    if svLo > svHi
        error('The SVD lower bound (%d) is above the upper bound (%d frames in the block)', svLo, svHi)
    end

    % Temporal covariance of the whole block, summed over voxel chunks on the GPU
    vc = max(1, floor(prm.gpuChunkBytes / (8 * nfB)));
    CM = zeros(nfB, nfB);
    for v0 = 1:vc:nvox
        v1 = min(v0 + vc - 1, nvox);
        Pc = gpuArray(IQblock(v0:v1, :));
        CM = CM + double(gather(Pc' * Pc));
    end
    clear Pc
    CM = (CM + CM') / 2; % exactly Hermitian, so eig takes the Hermitian path
    [V, ev] = eig(CM, 'vector'); % ascending, real
    clear CM
    [~, order] = sort(real(ev), 'descend');
    V = V(:, order);

    % Clutter-filter projector Q = V_sel*V_sel'. If most components are kept, build it as I - V_rm*V_rm' instead
    nKeep = svHi - svLo + 1;
    if nfB - nKeep < nKeep
        V_rm = V(:, [1:svLo-1, svHi+1:nfB]);
        Q = eye(nfB) - V_rm * V_rm';
    else
        V_sel = V(:, svLo:svHi);
        Q = V_sel * V_sel';
    end
    clear V

    if prm.hpfOnBlock
        % causal HPF as a matrix (row j = response to an impulse at frame j), combined with Q, so the data is multiplied once
        Hm = filter(prm.HPF.b, prm.HPF.a, eye(nfB), [], 2);
        T = Q * Hm;
        clear Hm
    else
        T = Q;
    end
    clear Q
    T = gpuArray(single(T));

    for v0 = 1:vc:nvox % in place, voxel chunk by voxel chunk
        v1 = min(v0 + vc - 1, nvox);
        IQblock(v0:v1, :) = gather(gpuArray(IQblock(v0:v1, :)) * T);
    end
    clear T Pc
    fprintf('  SVD clutter filter on the whole %d-frame block (%.1f s)\n', nfB, toc(tFilt))

    %% 3-4. Re-slice into new superframes, PDI and g1 of each
    tProc = tic;
    sliceLen = prm.sliceLen;
    nSlices = floor(nfB / sliceLen);
    if nSlices < 1
        info.reason = sprintf('the block has %d frames, fewer than the re-sliced superframe length %d', nfB, sliceLen);
        return
    end
    if nfB > nSlices * sliceLen
        fprintf('  %d frame(s) at the end of the block do not fill a whole slice and are not used\n', nfB - nSlices * sliceLen)
    end
    nTau = prm.nTau;
    PDIv = zeros(nvox, nSlices, 'single');
    g1v = zeros(nvox, nTau, nSlices, 'like', complex(single(0)));

    if ~prm.hpfOnBlock % HPF of each slice alone, as a matrix too
        Hs = gpuArray(single(filter(prm.HPF.b, prm.HPF.a, eye(sliceLen), [], 2)));
    end
    L = 2^nextpow2(sliceLen + nTau - 1); % FFT length in g1T_fft
    vc2 = max(1, floor(prm.gpuChunkBytes / (8 * L)));
    for s = 1:nSlices
        fr = (s - 1) * sliceLen + (1:sliceLen);
        for v0 = 1:vc2:nvox
            v1 = min(v0 + vc2 - 1, nvox);
            Pc = gpuArray(IQblock(v0:v1, fr));
            if ~prm.hpfOnBlock
                Pc = Pc * Hs;
            end
            PDIv(v0:v1, s) = gather(sum(abs(Pc).^2, 2) ./ sliceLen);
            g1v(v0:v1, :, s) = gather(g1T_fft(Pc, nTau)); % [voxel, lag]
        end
    end
    fprintf('  PDI and g1 of %d slices of %d frames (%.1f s)\n', nSlices, sliceLen, toc(tProc))

    PDI = reshape(PDIv, [ny, nx, nz, nSlices]);
    g1 = reshape(g1v, [ny, nx, nz, nTau, nSlices]);

    info.ok = true;
    info.grid = grid;
    info.nx = nx; info.ny = ny; info.zrange = zr;
    info.nfBlock = nfB;
    info.nSlices = nSlices;
    info.sliceStartFrame = (0:nSlices-1)' * sliceLen + 1;

end


%% The RF files of the folder, sorted by file number: struct array f (name, folder), the file numbers, and the format ('mat' or 'rf')
function [f, fileNum, fmt] = findRFFiles(datapath, P, fileFormat)

    fileFormat = lower(fileFormat);
    if ~ismember(fileFormat, {'auto', 'mat', 'rf'})
        error('Opts.fileFormat must be auto, mat or rf')
    end

    % .mat saves: RF-...-<file #>.mat
    fm = dir(fullfile(datapath, 'RF-*.mat'));
    tokm = regexp({fm.name}, '^RF-.*-(\d+)\.mat$', 'tokens', 'once');
    isM = ~cellfun(@isempty, tokm);
    fm = fm(isM);
    numM = cellfun(@(t) str2double(t{1}), tokm(isM));

    % binary FastSave files: <P.fsFileName><5 digit file #>.rf
    fr = struct('name', {}, 'folder', {});
    numR = [];
    if isfield(P, 'fsFileName') && ~isempty(P.fsFileName)
        fr = dir(fullfile(datapath, '*.rf'));
        tokr = regexp({fr.name}, ['^', regexptranslate('escape', P.fsFileName), '(\d+)\.rf$'], 'tokens', 'once');
        isR = ~cellfun(@isempty, tokr);
        fr = fr(isR);
        numR = cellfun(@(t) str2double(t{1}), tokr(isR));
    end

    switch fileFormat
        case 'auto'
            if ~isempty(fm)
                fmt = 'mat';
                if ~isempty(fr)
                    fprintf('The folder has both .mat and .rf RF files, using the .mat files (Opts.fileFormat = ''rf'' for the others)\n')
                end
            elseif ~isempty(fr)
                fmt = 'rf';
            else
                error('No RF files in %s (RF-*-<file #>.mat, or <P.fsFileName><file #>.rf)', datapath)
            end
        case 'mat'
            fmt = 'mat';
        case 'rf'
            fmt = 'rf';
            if ~isfield(P, 'fsFileName')
                error('P has no fsFileName, so the .rf file names are unknown')
            end
    end
    if strcmp(fmt, 'mat'), f = fm; fileNum = numM; else, f = fr; fileNum = numR; end
    if isempty(f)
        error('No .%s RF files in %s', fmt, datapath)
    end
    [fileNum, order] = sort(fileNum);
    f = f(order);

end


%% Options dialog (modal). Shows the options of Opts as the defaults. Returns the updated Opts, or [] if cancelled
function Opts = reconDialog(Opts, fileNum, P, wl, datapath, savepath, fmt)

    result = [];

    % Geometry for the voxel count (same formulas as RCA_DAS_rcabeam's Verasonics grid)
    L = P.Trans.spacingMm / 1e3 * P.numElements / wl; % probe width in wavelengths
    depthWl = P.endDepth - P.startDepth;
    if isfield(P.Trans, 'spacing')
        sp = P.Trans.spacing;
    else
        sp = P.Trans.spacingMm / 1e3 / wl;
    end
    presetNames = {'Half wavelength (0.5, 0.5, 0.5)', sprintf('Verasonics fUS PData (%.4g, %.4g, 0.5)', sp, sp), 'Custom'};
    presetVals = {[0.5, 0.5, 0.5], [sp, sp, 0.5]};
    pd0 = reshape(Opts.pixelSpacingWl, 1, 3);
    preset0 = 3;
    for i = 1:2
        if max(abs(pd0 - presetVals{i})) < 1e-6, preset0 = i; end
    end
    out0 = upper(Opts.output);
    if ~ismember(out0, {'OPW', 'RCCR'}), out0 = 'OPW'; end

    fig = uifigure('Name', 'rcabeam reconstruction', 'Position', [100, 100, 620, 400], ...
        'WindowStyle', 'modal', 'CloseRequestFcn', @(src, ~) uiresume(src));
    movegui(fig, 'center')
    g = uigridlayout(fig, [11, 2]);
    g.ColumnWidth = {190, '1x'};
    g.RowHeight = {24, 24, 24, 28, 28, 28, 24, 28, 28, '1x', 32};
    g.Padding = [16, 16, 16, 16];

    addLabel(g, 1, 'Raw data folder'); lbl(g, 1, datapath);
    addLabel(g, 2, 'Save folder'); lbl(g, 2, savepath);
    gaps = '';
    if numel(fileNum) < max(fileNum) - min(fileNum) + 1, gaps = ' (with gaps)'; end
    addLabel(g, 3, 'RF files found'); lbl(g, 3, sprintf('%d .%s files, file numbers %d to %d%s', numel(fileNum), fmt, min(fileNum), max(fileNum), gaps));

    addLabel(g, 4, 'File numbers to reconstruct');
    eFiles = uieditfield(g, 'text', 'Value', fileNumbersText(Opts, fileNum), ...
        'Tooltip', 'e.g. 1:20, or 3, 7, 10:12, or 1:2:9 (first:step:last), or "all"');
    eFiles.Layout.Row = 4; eFiles.Layout.Column = 2;

    addLabel(g, 5, 'Voxel spacing preset');
    ddPreset = uidropdown(g, 'Items', presetNames, 'ItemsData', 1:3, 'Value', preset0, 'ValueChangedFcn', @onPreset);
    ddPreset.Layout.Row = 5; ddPreset.Layout.Column = 2;

    addLabel(g, 6, 'Spacing x, y, z (wavelengths)');
    gs = uigridlayout(g, [1, 3]); gs.Layout.Row = 6; gs.Layout.Column = 2;
    gs.Padding = [0, 0, 0, 0]; gs.ColumnSpacing = 8;
    eSp = gobjects(1, 3);
    for i = 1:3
        eSp(i) = uieditfield(gs, 'numeric', 'Value', pd0(i), 'Limits', [0.01, Inf], 'ValueChangedFcn', @onSpacing);
    end

    addLabel(g, 7, 'Volume'); lVol = lbl(g, 7, '');

    addLabel(g, 8, 'Output');
    ddOut = uidropdown(g, 'Items', {'OPW sum (single IQ volume)', 'Separate C-R and R-C volumes (RCCR, for XDoppler)'}, ...
        'ItemsData', {'OPW', 'RCCR'}, 'Value', out0);
    ddOut.Layout.Row = 8; ddOut.Layout.Column = 2;

    addLabel(g, 9, 'Frames in the RF files');
    cbStack = uicheckbox(g, 'Text', 'Stacked in one buffer frame (apply unstackFrames)', 'Value', logical(Opts.unstackedFrames), ...
        'Tooltip', 'On: like RCA_VerasonicsRecon_loop_unstackedFrames (always for binary .rf files). Off: RcvData is already [samples, channels, frames]');
    cbStack.Layout.Row = 9; cbStack.Layout.Column = 2;

    gb = uigridlayout(g, [1, 3]); gb.Layout.Row = 11; gb.Layout.Column = [1, 2];
    gb.Padding = [0, 0, 0, 0]; gb.ColumnWidth = {'1x', 100, 120};
    bCancel = uibutton(gb, 'Text', 'Cancel', 'ButtonPushedFcn', @(~, ~) uiresume(fig), 'Tag', 'cancel');
    bCancel.Layout.Column = 2;
    bOK = uibutton(gb, 'Text', 'Reconstruct', 'ButtonPushedFcn', @onOK, 'Tag', 'ok', 'FontWeight', 'bold');
    bOK.Layout.Column = 3;

    updateVolume()
    uiwait(fig)
    if isvalid(fig), delete(fig), end
    Opts = result;


    function onPreset(~, ~)
        if ddPreset.Value <= 2
            for k = 1:3, eSp(k).Value = presetVals{ddPreset.Value}(k); end
        end
        updateVolume()
    end

    function onSpacing(~, ~)
        ddPreset.Value = 3;
        for k = 1:2
            if max(abs([eSp.Value] - presetVals{k})) < 1e-6, ddPreset.Value = k; end
        end
        updateVolume()
    end

    function updateVolume()
        v = [eSp.Value];
        n = [ceil(round(L / v(1), 9)), ceil(round(L / v(2), 9)), ceil(round(depthWl / v(3), 9))]; % x, y, z
        lVol.Text = sprintf('%d x %d x %d voxels (x, y, z), %.1f MB per volume (complex single)', n, prod(n) * 8 / 1e6);
    end

    function onOK(~, ~)
        try
            nums = parseFileNumbers(eFiles.Value, fileNum);
        catch err
            uialert(fig, err.message, 'File numbers')
            return
        end
        if isempty(intersect(nums, fileNum))
            uialert(fig, sprintf('None of these file numbers has an RF file (found %d to %d).', min(fileNum), max(fileNum)), 'File numbers')
            return
        end
        result = Opts;
        result.fileNumbers = nums;
        result.pixelSpacingWl = [eSp.Value];
        result.output = ddOut.Value;
        result.unstackedFrames = cbStack.Value;
        uiresume(fig)
    end

end


function addLabel(g, row, txt)

    l = uilabel(g, 'Text', txt, 'HorizontalAlignment', 'right');
    l.Layout.Row = row; l.Layout.Column = 1;

end


function l = lbl(g, row, txt)

    l = uilabel(g, 'Text', txt, 'Tooltip', txt);
    l.Layout.Row = row; l.Layout.Column = 2;

end


%% The file numbers shown at first: the ones set in Opts, or the files inside Opts.fileRange
function s = fileNumbersText(Opts, fileNum)

    if ~isempty(Opts.fileNumbers)
        n = Opts.fileNumbers(:)';
    else
        n = fileNum(fileNum >= Opts.fileRange(1) & fileNum <= Opts.fileRange(2));
        n = n(:)';
    end
    if isempty(n)
        n = fileNum(:)';
    end
    if numel(n) > 1 && all(diff(n) == 1)
        s = sprintf('%d:%d', n(1), n(end));
    else
        s = char(strjoin(string(n), ', '));
    end

end


%% File numbers from text like '1:20', '3, 7, 10:12' or '1:2:9' (first:step:last), or 'all'. Errors on anything else
function n = parseFileNumbers(txt, available)

    txt = lower(strtrim(txt));
    if isempty(txt) || strcmp(txt, 'all')
        n = available(:)';
        return
    end
    n = [];
    parts = strsplit(txt, {',', ';', ' '}, 'CollapseDelimiters', true);
    for i = 1:numel(parts)
        t = parts{i};
        if isempty(t), continue, end
        % (a regexp with optional groups and 'tokens' only returns the first number here, so split at the colons instead)
        if isempty(regexp(t, '^\d+(:\d+){0,2}$', 'once'))
            error('Can''t read "%s". Use file numbers like 1:20 or 3, 7, 10:12.', t)
        end
        v = str2double(strsplit(t, ':'));
        switch numel(v)
            case 1
                n = [n, v]; %#ok<AGROW>
            case 2
                n = [n, v(1) : v(2)]; %#ok<AGROW>
            case 3
                n = [n, v(1) : v(2) : v(3)]; %#ok<AGROW>
        end
    end
    n = unique(n);

end


%% Actual frame rate from the RF time tags of the acquisition (same calculation as checkFrameRate.m, which opens the file with uiopen
% into the base workspace). Falls back to the nominal P.frameRate if the time tag file is not in the data folder
function frameRate = frameRateFromTimeTags(P, datapath)

    ttFile = fullfile(datapath, 'RFTimeTagsPerChannel.mat');
    if ~isfile(ttFile)
        warning('%s not found, using the nominal frame rate of %g Hz', ttFile, P.frameRate)
        frameRate = P.frameRate;
        return
    end
    S = load(ttFile, 'RFTimeTags_allSFStacked');
    dt_all = diff(squeeze(mean(S.RFTimeTags_allSFStacked, 1)), 1); % Difference in frame timetags
    dt_subframe = dt_all(dt_all < 0.8*(1/P.sfRate - P.numFramesPerBuffer/P.frameRate)); % Remove the large inter-superframe jumps
    frameRate_actual = 1 / mean(dt_subframe);
    if abs(P.frameRate - frameRate_actual) > 10 % hardware could lag behind the nominal frame rate
        warning('Nominal frame rate (%g Hz) was not supported by the hardware during data acquisition. Using the actual frame rate of %g Hz instead.', P.frameRate, frameRate_actual)
        frameRate = frameRate_actual;
    else
        frameRate = P.frameRate;
    end

end


function S = setDefault(S, name, val)
    if ~isfield(S, name)
        S.(name) = val;
    end
end
