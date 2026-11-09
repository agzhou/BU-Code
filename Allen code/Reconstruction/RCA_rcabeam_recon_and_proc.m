%% Description

% Reconstruct the raw RF data of a row-column array (RC15gV) acquisition folder with the GPU rcabeam library
% (RCA_DAS_rcabeam), one IQ file per RF file. Counterpart of RCA_ffdas_recon.m and of the Verasonics recon loops
% (RCA_VerasonicsRecon_loop.m, RCA_VerasonicsRecon_loop_unstackedFrames.m), without needing the Verasonics software.

% AND, process the data as in IQ2g1T_3D_loop.m

% Input folder: params.mat (structure P) and either
%           RF-<maxAngle>-<na>-<frameRate>-<numFrames>-1-<file #>.mat   (RcvData)         Verasonics .mat saves, or
%           <P.fsFileName><file #, 5 digits>.rf                          (binary, int16)   Verasonics FastSave
%         (see Opts.fileFormat)
% Output: one file per RF file in the save folder, named like the Verasonics recon output (the 'RF-' prefix is replaced; the
%         .rf files get IQ-<maxAngle>-<na>-<frameRate>-<numFrames>-1-<file #>.mat, as in the Verasonics loop for them):
%           'OPW'    IQ-...       IQ [y, x, z, frame], all 2*na acquisitions summed coherently (one pass on the GPU)
%           'RCCR'   IQ-RCCR-...  IQ_CR and IQ_RC [y, x, z, frame], each compounded over the angles
%         and PData.mat and rcabeam_recon_params.mat (the pixel grid in m, and the options used) in the save folder

% Notes:
%           - Setup is the same as RCA_DAS_rcabeam: MATLAB R2026b+ with the Parallel Computing Toolbox, and the python.exe
%             of an environment with cupy and rcabeam in Opts.pythonExe or the RCABEAM_PYTHON environment variable. Run
%             test_rcabeam_layout.m and test_RCA_DAS_rcabeam_psf.m once before trusting a recon
%           - Binary .rf files (FastSave) are read with loadBinaryRF, which gives the stacked [rowsPerFrame, colsPerFrame]
%             array, then unstackFrames, as in RCA_VerasonicsRecon_loop_unstackedFrames_fastsave.m. They are always stacked, so
%             Opts.unstackedFrames should stay true for them (false uses the array as is). The file numbers are the numbers in
%             the names (FastSave starts at 0). The channels are remapped with P.Trans.Connector inside RCA_DAS_rcabeam, as for
%             the .mat files
%           - Opts.unstackedFrames chooses how the .mat RF files are read (like the two Verasonics loops):
%               true (default)  each file holds the frames stacked in one buffer frame (RcvData is 2-D, or has a few
%                               stacked chunks), and is split with unstackFrames, as in RCA_VerasonicsRecon_loop_unstackedFrames.m
%               false           each file holds RcvData as [samples per frame, channels, frames] already, as in
%                               RCA_VerasonicsRecon_loop.m. Used as is
%           - rcabeam has no per-angle output, so there is no 'AllAngles' output (use RCA_ffdas_recon.m for that), and
%             no Trans.ElementSens weighting. See RCA_DAS_rcabeam for the other differences to RCA_DAS
%           - The pixel grid is the Verasonics recon's PData (pixel 1 at PData.Origin = [-halfProbeDist, +halfProbeDist,
%             startDepth], x increasing and y decreasing with the index), saved as PData.mat like the Verasonics loops, so
%             the volumes line up voxel by voxel with the Verasonics recon (use pixelSpacingWl = [Trans.spacing,
%             Trans.spacing, 0.5], 0.9714 wavelengths for the RC15gV, for the PData of the fUS runs). The y direction
%             follows from the PData definition, and is only checked by the simulated PSF at x = y = 0. Opts.das.grid =
%             'centered' gives RCA_DAS's centered grid (y increasing) instead, and then no PData.mat is written
%           - The IQ is complex single, demodulated to absolute time (see RCA_DAS_rcabeam), so its phase reference is not the
%             Verasonics one. Its magnitude is not scaled like the Verasonics IQ
%           - The GPU workspace is kept between files (same geometry), and released at the end

% Inputs (all optional, prompts for the folders if they are empty):
%   datapath   folder with params.mat and the RF files
%   savepath   folder to save to (default: 'IQ Data - rcabeam' inside datapath)
%   Opts fields:
%     gui              a dialog asks for the file numbers, voxel spacing, output type and whether the frames are stacked.
%                      It is shown when Opts is not given at all. With an Opts it is shown only if gui is true, and then the
%                      values in Opts are its defaults. Cancelling returns without reconstructing
%     fileRange        [first, last] RF file numbers to reconstruct (default [0, Inf]: all files, including file 0)
%     fileNumbers      list of RF file numbers to reconstruct, e.g. [1:10, 15] (default [], which uses fileRange)
%     pixelSpacingWl   [x, y, z] pixel spacing in wavelengths (default [0.5, 0.5, 0.5])
%     output           'OPW' (default) or 'RCCR'
%     unstackedFrames  true (default) or false, see the Notes
%     fileFormat       'auto' (default), 'mat' or 'rf'. auto: the .mat files if the folder has RF-*.mat files, else the .rf files
%                      named by P.fsFileName. If it has both, auto uses the .mat files ('rf' to use the others)
%     overwrite        false: skip files that already have an output (default)
%     das              struct of options passed to RCA_DAS_rcabeam (frameChunk, fNumber, pythonExe, backend, ...)
%
% Output: cell array with the paths of the saved files

%% Add the rcabeam/processing folders (and subfolders) to path
codeDir = cd;
codeDir_split = split(string(codeDir), filesep);
rcabeamCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "GitHub"))), '\') + "\rcabeam\");
addpath(genpath(rcabeamCodePath))
processingbeamCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "Allen code"))), '\') + "\Processing\");
addpath(genpath(processingbeamCodePath))

% function IQfiles = RCA_rcabeam_recon_and_proc(datapath, PDsavepath, Opts)

    % if nargin < 1 || isempty(datapath)
        datapath = uigetdir(pwd, 'Select the raw data path');
        if isequal(datapath, 0)
            error('No raw data path selected.')
        end
    % end
    % if nargin < 2 || isempty(IQsavepath)
    %     % IQsavepath = fullfile(datapath, 'IQ Data - rcabeam');
    %     IQsavepath = uigetdir(pwd, 'Select the IQ save path');
    %     if isequal(IQsavepath, 0)
    %         error('No IQ save path selected.')
    %     end
    % end
    % if nargin < 2 || isempty(PDsavepath)
        % PDsavepath = fullfile(datapath, 'Processed Data - rcabeam');
        PDsavepath = uigetdir(datapath, 'Select the save path for processed data');
        if isequal(PDsavepath, 0)
            error('No save path selected.')
        end
    % end
    % if nargin < 3
        Opts = struct();
        showGUI = true; % called without options: ask for them
    % else
    %     showGUI = isfield(Opts, 'gui') && Opts.gui;
    % end
    Opts = setDefault(Opts, 'fileRange', [0, Inf]);
    Opts = setDefault(Opts, 'fileNumbers', []);
    Opts = setDefault(Opts, 'pixelSpacingWl', [0.5, 0.5, 0.5]);
    Opts = setDefault(Opts, 'output', 'OPW');
    Opts = setDefault(Opts, 'unstackedFrames', true);
    Opts = setDefault(Opts, 'fileFormat', 'auto');
    Opts = setDefault(Opts, 'overwrite', false);
    Opts = setDefault(Opts, 'das', struct());
    % Opts = setDefault(Opts, 'grid', 'verasonics');
    % Opts = setDefault(Opts, 'gpuOutput', true);
    Opts.das.gpuOutput = true; % Default is to keep IQ on the GPU

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
    procParamsPrompt = {'SVD lower bound', 'SVD upper bound', 'z start', 'z end', 'High pass filter cutoff [Hz]', 'Time lag window [ms]'};
    procParamsDefaults = {'50', num2str(P.numFramesPerBuffer), '1', '110', '25', '20'};
    procParamsUserInput = inputdlg(procParamsPrompt, 'Input Parameters', 1, procParamsDefaults);

    % define # of files manually for now
    % str2double(parameterUserInput{});
    sv_threshold_lower = str2double(procParamsUserInput{1});
    sv_threshold_upper = str2double(procParamsUserInput{2});
    zstart = str2double(procParamsUserInput{3});
    zend = str2double(procParamsUserInput{4});
    tau_window_s = str2double(procParamsUserInput{6}) / 1e3;

    % Set up the High Pass Filter (parameters from the 2020 vUS paper)
    HPF.fc = str2double(procParamsUserInput{5}); % Cutoff frequency [Hz]
    % 25 Hz corresponds to 1 mm/s
    HPF.order = 4; % Butterworth filter order
    [frameRate] = checkFrameRate(P);
    HPF.fs = frameRate; % Sampling frequency [Hz]
    [HPF.b, HPF.a] = butter(HPF.order, HPF.fc/(HPF.fs/2), 'high');

    % zstart = 1
    % zend = 110
    % sv_threshold_lower = 50;
    % sv_threshold_upper = size(IQ, 4);

    

    nTau = ceil(tau_window_s *frameRate); % # of time lags to consider; empirically set by assuming all g1 for voxels containing actual flow decay within 10 ms
    % nTau = ceil(100e-3 *P.frameRate); % # of time lags to consider; empirically set by assuming all g1 for voxels containing actual flow decay within 10 ms
    tau = (0:nTau - 1)' ./ frameRate; % Time lag vector [s]
    clearvars procParamsPrompt procParamsDefaults procParamsUserInput

    %% Find the RF files, in order of file number
    [f, fileNum, fmt] = findRFFiles(datapath, P, Opts.fileFormat);
    fprintf('%d RF files (.%s) in %s\n', numel(f), fmt, datapath)

    %% Options from the dialog
    if showGUI
        Opts = reconDialog(Opts, fileNum, P, wl, datapath, PDsavepath, fmt);
        if isempty(Opts)
            fprintf('Reconstruction cancelled.\n')
            IQfiles = {};
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
            namePrefix = 'IQ-';
            Opts.das.sumOnly = true;
        case 'rccr'
            error('Downstream processing is currently not built for RCCR')
            namePrefix = 'IQ-RCCR-';
            Opts.das.sumOnly = false;
        case 'allangles'
            error('rcabeam has no per-angle output. Use RCA_ffdas_recon with Opts.output = ''AllAngles''.')
        otherwise
            error('Unknown output type "%s". Use OPW or RCCR.', Opts.output)
    end

    useGPU = ~isfield(Opts.das, 'backend') || strcmpi(Opts.das.backend, 'rcabeam');
    if useGPU && gpuDeviceCount == 0
        error('No supported GPU found. The rcabeam backend needs an NVIDIA GPU.')
    end

    pix = Opts.pixelSpacingWl * wl; % [x, y, z] in m

    if ~exist(PDsavepath, 'dir')
        mkdir(PDsavepath)
    end

    %% Reconstruct
    cleanup = onCleanup(@() RCA_DAS_rcabeam('clear')); % release the GPU workspace, also after an error or Ctrl+C
    IQfiles = {};
    for k = 1:numel(f) % Loop through files
        if strcmp(fmt, 'rf')
            % like the Verasonics loops: IQ-<maxAngle>-<na>-<frameRate>-<numFramesPerBuffer>-1-<file #>
            outName = sprintf('%s%d-%d-%d-%d-1-%d.mat', namePrefix, round(P.maxAngle), P.na, round(P.frameRate), P.numFramesPerBuffer, fileNum(k));
        else
            outName = [namePrefix, f(k).name(4:end)]; % replace the 'RF-' prefix
        end
        outFile = fullfile(PDsavepath, outName);
        if exist(outFile, 'file') && ~Opts.overwrite
            fprintf('%s already exists, skipping.\n', outName)
            IQfiles{end+1} = outFile; %#ok<AGROW>
            continue
        end

        tic
        if strcmp(fmt, 'rf')
            RcvData = loadBinaryRF(fullfile(f(k).folder, f(k).name), P); % stacked [rowsPerFrame, colsPerFrame], int16
        else
            S = load(fullfile(f(k).folder, f(k).name), 'RcvData');
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

        % RCA_DAS_rcabeam takes the frame count from P.numSubFrames when that exists, so make sure it matches this file
        Pk = P;
        if isfield(Pk, 'numSubFrames')
            Pk.numSubFrames = size(r, 3);
        end

        [IQ_CR, IQ_RC, grid] = RCA_DAS_rcabeam(r, Pk, pix(2), pix(1), pix(3), Opts.das);
        clear r

        % Verasonics IQ files are [y, x, z, frame]; RCA_DAS_rcabeam is [x, y, z, frame]
        switch outputMode
            case 'opw' % with sumOnly, IQ_CR holds the coherent sum of both
                IQ = permute(IQ_CR, [2, 1, 3, 4]);
                % save(outFile, 'IQ', '-v7.3', '-nocompression') % Don't save IQ
            case 'rccr'
                IQ_CR = permute(IQ_CR, [2, 1, 3, 4]);
                IQ_RC = permute(IQ_RC, [2, 1, 3, 4]);
                % save(outFile, 'IQ_CR', 'IQ_RC', '-v7.3', '-nocompression')
        end
        clear IQ_CR IQ_RC

        IQfiles{end+1} = outFile; %#ok<AGROW>
        % fprintf('RF file %d reconstructed and saved as %s (%.1f s)\n', fileNum(k), outName, toc)
        fprintf('RF file %d reconstructed (%.1f s)\n', fileNum(k), toc)

        % ==== TESTING ==== %
        
        [xpo, ypo, zpo, nfo] = size(IQ); % Original sizes
        % figure; imagesc(squeeze(max(abs(IQ(:, :, :, 1)), [], 1)).')
        xrange = 1:xpo;
        yrange = 1:ypo;
        % zrange = 1:100;
        % zrange = 1:zpo;
        zrange = zstart:zend;
        voxelRange = {xrange, yrange, zrange};
        xp = length(xrange); yp = length(yrange); zp = length(zrange); nf = size(IQ, 4);

        tic
        [PDI, g1] = IQ2g1T_3D_nofreqsep(IQ, voxelRange, sv_threshold_lower, sv_threshold_upper, HPF, nTau);
        % toc
    
        % voxelTimeseriesGUI(gather(g1), gather(abs(g1(:, :, :, 2))), 'ComplexMode', 'abs', 'Colormap', 'turbo')
        % save(fullfile(PDsavepath, 'fUSdata-', num2str(k)', '.mat'), 'PDI', 'g1', '-v7.3', '-nocompression')
        
        % Get the file index from the RF name
        if strcmp(fmt, 'rf')
            temp = strsplit(f(k).name, '-');
            temp2 = strsplit(temp{end}, '.');
            if strcmp(temp2{1}, '00000')
                fileInd = 0;
            else
                fileInd = strip(temp2{1}, 'left', '0');
            end
        else
            error('Code is not set up to get the file index properly for non FastSave data')
        end
        save([PDsavepath, '\fUSdata-', num2str(fileInd), '.mat'], 'PDI', 'g1', '-v7.3', '-nocompression')
        fprintf('Processed data file %d saved (%.1f s)\n', fileNum(k), toc)

        clearvars IQ PDI g1
    end

    if exist('grid', 'var') % not set if every file was skipped
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
        save(fullfile(PDsavepath, 'fUS_proc_params.mat'), 'sv_threshold_lower', 'sv_threshold_upper', 'HPF', 'voxelRange', 'zstart', 'zend', 'tau', 'nTau', 'tau_window_s', 'xp', 'yp', 'zp', 'frameRate') % grid (m) is the same for every file
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


function S = setDefault(S, name, val)
    if ~isfield(S, name)
        S.(name) = val;
    end
end
