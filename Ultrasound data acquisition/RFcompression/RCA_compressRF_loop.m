%% Description

% Compress the binary RF files (Verasonics FastSave, <P.fsFileName><5 digit file #>.rf) of a row-column array (RC15gV)
% acquisition folder into .mat files that keep only the receive channels the RCA beamformers use (see compressRCA_RF).
% About half the size of the .rf files, before any HDF5 compression.

% Input folder: params.mat (structure P) and the .rf files
% Output folder: RFc-<file #>.mat for each .rf file (same file # as in the .rf name, FastSave starts at 0), with
%                  RcvData  int16, (# samples per acquisition, 2*na, numElements, # frames), see compressRCA_RF
%                  layout   struct describing RcvData
%                  srcFile  name of the .rf file it came from
%                and a copy of params.mat
% Read by RCA_rcabeam_recon_and_proc.m (Opts.fileFormat 'rfc', or 'auto' picks it up) and RCA_DAS_rcabeam directly.
% NOT readable by the Verasonics recon loops (RCA_VerasonicsRecon_loop*.m), which need the full buffer.

clearvars
close all

%% Options
overwrite = false;      % false: skip files that already have an RFc file
hdf5Compression = false; % true: gzip inside the -v7.3 file (smaller, slower to save and load). false: '-nocompression'
fileNumbers = [];       % default shown in the file selection dialog: [] for all files, or a list, e.g. 0:99
deleteSource = false;   % true: delete each .rf file after its RFc file is saved and read back. Not undoable

%% Folders
RFPath = uigetdir('C:\', 'Select the raw .rf data folder (with params.mat)');
if isequal(RFPath, 0), error('No RF folder selected.'), end
outPath = uigetdir(RFPath, 'Select the folder to save the compressed RF');
if isequal(outPath, 0), error('No save folder selected.'), end
if strcmp(RFPath, outPath)
    warning('The save folder is the RF folder. The recon''s auto file format picks the RFc files there over the .rf files.')
end

%% Parameters and files
S = load(fullfile(RFPath, 'params.mat'), 'P');
P = S.P;
clear S
if ~isequal(P.Trans.name, 'RC15gV')
    error('This is for the RC15gV row-column array, not %s', P.Trans.name)
end
if ~isfield(P, 'numFramesPerBuffer') && isfield(P, 'numFramesPerBuff') % older P structures
    P.numFramesPerBuffer = P.numFramesPerBuff;
end

fr = dir(fullfile(RFPath, '*.rf'));
tok = regexp({fr.name}, ['^', regexptranslate('escape', P.fsFileName), '(\d+)\.rf$'], 'tokens', 'once');
isR = ~cellfun(@isempty, tok);
fr = fr(isR);
fileNum = cellfun(@(t) str2double(t{1}), tok(isR));
[fileNum, order] = sort(fileNum);
fr = fr(order);
if isempty(fr)
    error('No %s<file #>.rf files in %s', P.fsFileName, RFPath)
end
fprintf('%d .rf files in %s, file numbers %d to %d\n', numel(fr), RFPath, min(fileNum), max(fileNum))

%% Choose the files to compress
if isempty(fileNumbers)
    defaultText = sprintf('%d:%d', min(fileNum), max(fileNum));
else
    defaultText = strjoin(string(fileNumbers(:)'), ', ');
end
gaps = '';
if numel(fileNum) < max(fileNum) - min(fileNum) + 1, gaps = ', with gaps'; end
prompt = sprintf(['Found %d .rf files, file numbers %d to %d%s.\n' ...
    'File numbers to compress (e.g. 0:99, or 3, 7, 10:12, or 0:2:20, or all):'], numel(fileNum), min(fileNum), max(fileNum), gaps);
while true
    answer = inputdlg(prompt, 'Files to compress', [1, 70], {char(defaultText)});
    if isempty(answer)
        error('Compression cancelled.')
    end
    try
        fileNumbers = parseFileNumbers(answer{1}, fileNum);
    catch err
        uiwait(errordlg(err.message, 'File numbers'))
        defaultText = answer{1};
        continue
    end
    missing = setdiff(fileNumbers, fileNum);
    sel = ismember(fileNum, fileNumbers);
    if ~any(sel)
        uiwait(errordlg(sprintf('None of these file numbers has an .rf file (found %d to %d).', min(fileNum), max(fileNum)), 'File numbers'))
        defaultText = answer{1};
        continue
    end
    break
end
if ~isempty(missing)
    warning('No .rf file for file number(s) %s, skipping them', mat2str(missing(:)'))
end
fr = fr(sel);
fileNum = fileNum(sel);
fprintf('Compressing %d files: %s\n', numel(fileNum), answer{1})

if ~exist(fullfile(outPath, 'params.mat'), 'file')
    copyfile(fullfile(RFPath, 'params.mat'), fullfile(outPath, 'params.mat'))
end

if hdf5Compression
    saveArgs = {'-v7.3'};
else
    saveArgs = {'-v7.3', '-nocompression'};
end

%% Compress
bytesIn = 0;
bytesOut = 0;
for k = 1:numel(fr)
    outFile = fullfile(outPath, sprintf('RFc-%d.mat', fileNum(k)));
    if exist(outFile, 'file') == 2 && ~overwrite
        fprintf('RFc-%d.mat already exists, skipping.\n', fileNum(k))
        continue
    end

    tic
    srcFile = fr(k).name;
    RcvData = loadBinaryRF(fullfile(fr(k).folder, srcFile), P); % stacked [rowsPerFrame, colsPerFrame], int16
    RcvData = unstackFrames(RcvData, P);                          % (# samples per frame, # channels, # frames)
    [RcvData, layout] = compressRCA_RF(RcvData, P);               % (# samples per acquisition, 2*na, numElements, # frames)
    save(outFile, 'RcvData', 'layout', 'srcFile', saveArgs{:})

    inInfo = dir(fullfile(fr(k).folder, srcFile));
    outInfo = dir(outFile);
    bytesIn = bytesIn + inInfo.bytes;
    bytesOut = bytesOut + outInfo.bytes;
    fprintf('%s -> RFc-%d.mat, %.0f -> %.0f MB (%.2fx), %.1f s\n', srcFile, fileNum(k), inInfo.bytes / 1e6, ...
        outInfo.bytes / 1e6, inInfo.bytes / outInfo.bytes, toc)

    if deleteSource
        % Only delete once the saved file reads back identical
        chk = load(outFile, 'RcvData');
        if isequal(chk.RcvData, RcvData)
            delete(fullfile(fr(k).folder, srcFile))
        else
            warning('RFc-%d.mat does not read back identical, keeping %s', fileNum(k), srcFile)
        end
        clear chk
    end
    clear RcvData
end
if bytesOut > 0
    fprintf('Total: %.2f -> %.2f GB (%.2fx)\n', bytesIn / 1e9, bytesOut / 1e9, bytesIn / bytesOut)
end


%% File numbers from text like '0:99', '3, 7, 10:12' or '0:2:20' (first:step:last), or 'all'. Errors on anything else
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
        if isempty(regexp(t, '^\d+(:\d+){0,2}$', 'once'))
            error('Can''t read "%s". Use file numbers like 0:99 or 3, 7, 10:12.', t)
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
