%% Description

% Reconstruct the raw RF data of a row-column array (RC15gV) acquisition folder with ffdas (RCA_DAS_ffdas), one IQ file
% per RF file. GPU counterpart of the Verasonics recon loop (RCA_VerasonicsRecon_loop_unstackedFrames*.m), without
% needing the Verasonics software.

% Input folder: params.mat (structure P) and RF-<maxAngle>-<na>-<frameRate>-<numFrames>-1-<file #>.mat (RcvData),
%               e.g., the folders in G:\My Drive\Data\PSF Simulations
% Output: one file per RF file in the save folder, named like the Verasonics recon output (the 'RF-' prefix is replaced):
%           'OPW'        IQ-...       IQ [y, x, z, frame], all 2*na acquisitions summed coherently
%           'RCCR'       IQ-RCCR-...  IQ_CR and IQ_RC [y, x, z, frame], each compounded over the angles
%           'AllAngles'  IQ-...       IQ [y, x, z, page, frame]. Pages 1:na are C-R and na+1:2*na are R-C, in acquisition order
%         and ffdas_recon_params.mat with the pixel grid (m) and the options used

% Notes:
%           - The IQ is single precision and is the analytic signal of the RF at each voxel's round trip time (see
%             RCA_DAS_ffdas), so its magnitude is comparable to the Verasonics IQ but its phase reference is not the
%             Verasonics demodulation reference. Phase differences between frames/volumes are not affected
%           - Pixel (x, y, z) positions are the pixel centers of the RCA_DAS grid, not the PData corners
%             Verasonics uses, so there can be a sub-pixel shift relative to the Verasonics IQ
%           - Not checked: whether the y index runs in the same direction as in the Verasonics IQ. Compare on a
%             scatterer that is off-center in y before using the two interchangeably

% Inputs (all optional, prompts for the folders if they are empty):
%   datapath   folder with params.mat and the RF files
%   savepath   folder to save to (default: 'IQ Data - ffdas' inside datapath)
%   Opts fields:
%     fileRange       [first, last] RF file numbers to reconstruct (default [1, Inf])
%     pixelSpacingWl  [x, y, z] pixel spacing in wavelengths (default [0.5, 0.5, 0.5])
%     output          'OPW' (default), 'RCCR' or 'AllAngles'
%     overwrite       false: skip files that already have an output (default)
%     das             struct of options passed to RCA_DAS_ffdas (exactRx, rxHalfAngle, frameChunk, ...)
%
% Output: cell array with the paths of the saved files

function IQfiles = RCA_ffdas_recon(datapath, savepath, Opts)

    if nargin < 1 || isempty(datapath)
        datapath = uigetdir(pwd, 'Select the raw data path');
        if isequal(datapath, 0)
            error('No raw data path selected.')
        end
    end
    if nargin < 2 || isempty(savepath)
        savepath = fullfile(datapath, 'IQ Data - ffdas');
    end
    if nargin < 3
        Opts = struct();
    end
    Opts = setDefault(Opts, 'fileRange', [1, Inf]);
    Opts = setDefault(Opts, 'pixelSpacingWl', [0.5, 0.5, 0.5]);
    Opts = setDefault(Opts, 'output', 'OPW');
    Opts = setDefault(Opts, 'overwrite', false);
    Opts = setDefault(Opts, 'das', struct());
    Opts.das = setDefault(Opts.das, 'frameChunk', 16); % frames beamformed per batch

    outputMode = lower(Opts.output);
    switch outputMode
        case 'opw'
            namePrefix = 'IQ-';
            Opts.das.perAngle = false;
        case 'rccr'
            namePrefix = 'IQ-RCCR-';
            Opts.das.perAngle = false;
        case 'allangles'
            namePrefix = 'IQ-';
            Opts.das.perAngle = true;
        otherwise
            error('Unknown output type "%s". Use OPW, RCCR or AllAngles.', Opts.output)
    end

    useGPU = ~isfield(Opts.das, 'useGPU') || Opts.das.useGPU;
    if useGPU && gpuDeviceCount == 0
        error('No supported GPU found. ffdas needs an NVIDIA GPU.')
    end

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
    pix = Opts.pixelSpacingWl * wl; % [x, y, z] in m

    %% Find the RF files, in order of file number
    f = dir(fullfile(datapath, 'RF-*.mat'));
    tok = regexp({f.name}, '^RF-.*-(\d+)\.mat$', 'tokens', 'once');
    isRF = ~cellfun(@isempty, tok);
    f = f(isRF);
    fileNum = cellfun(@(t) str2double(t{1}), tok(isRF));
    [fileNum, order] = sort(fileNum);
    f = f(order);
    inRange = fileNum >= Opts.fileRange(1) & fileNum <= Opts.fileRange(2);
    f = f(inRange);
    if isempty(f)
        error('No RF files in %s in the file range [%g, %g]', datapath, Opts.fileRange(1), Opts.fileRange(2))
    end

    if ~exist(savepath, 'dir')
        mkdir(savepath)
    end

    %% Reconstruct
    IQfiles = {};
    for k = 1:numel(f)
        outName = [namePrefix, f(k).name(4:end)]; % replace the 'RF-' prefix
        outFile = fullfile(savepath, outName);
        if exist(outFile, 'file') && ~Opts.overwrite
            fprintf('%s already exists, skipping.\n', outName)
            IQfiles{end+1} = outFile; %#ok<AGROW>
            continue
        end

        tic
        S = load(fullfile(f(k).folder, f(k).name), 'RcvData');
        r = unstackFrames(S.RcvData, P); % (# samples per frame, # channels, # frames)
        clear S

        [IQ_CR, IQ_RC, grid] = RCA_DAS_ffdas(r, P, pix(2), pix(1), pix(3), Opts.das);
        clear r

        % Verasonics IQ files are [y, x, z, (page,) frame]; RCA_DAS_ffdas is [x, y, z, (angle,) frame]
        switch outputMode
            case 'opw'
                IQ = permute(IQ_CR + IQ_RC, [2, 1, 3, 4]);
                save(outFile, 'IQ', '-v7.3', '-nocompression')
            case 'rccr'
                IQ_CR = permute(IQ_CR, [2, 1, 3, 4]);
                IQ_RC = permute(IQ_RC, [2, 1, 3, 4]);
                save(outFile, 'IQ_CR', 'IQ_RC', '-v7.3', '-nocompression')
            case 'allangles'
                IQ = permute(cat(4, IQ_CR, IQ_RC), [2, 1, 3, 4, 5]);
                save(outFile, 'IQ', '-v7.3', '-nocompression')
        end
        clear IQ IQ_CR IQ_RC

        IQfiles{end+1} = outFile; %#ok<AGROW>
        fprintf('RF file %d reconstructed and saved as %s (%.1f s)\n', fileNum(k), outName, toc)
    end

    if exist('grid', 'var') % not set if every file was skipped
        ReconOpts = Opts; %#ok<NASGU>
        save(fullfile(savepath, 'ffdas_recon_params.mat'), 'grid', 'ReconOpts') % grid is the same for every file
    end

end


function S = setDefault(S, name, val)
    if ~isfield(S, name)
        S.(name) = val;
    end
end
