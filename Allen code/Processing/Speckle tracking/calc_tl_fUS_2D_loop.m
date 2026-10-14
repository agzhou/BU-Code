%% Description: tl-fUS processing for 2D fUS data
clearvars

%% Add the Processing folder (and subfolders) to path
codeDir = cd;
codeDir_split = split(string(codeDir), filesep);
% AllenVerasonicsCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "Allen code"))), '\') + "\Verasonics");
AllenProcessingCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "BU-Code"))), '\') + "\Allen Code\Processing\");
addpath(genpath(AllenProcessingCodePath))

%% Choose processed data directory, params, and RF timetags data file
if ~exist('PDpath', 'var')
    PDpath = uigetdir('D:\Allen\Data\', 'Select the Processed data (g1) path');
    PDpath = [PDpath, '\'];
end
load([PDpath, 'fUS_proc_params.mat'])

% Load acquisition parameters: params.mat
if ~exist('P', 'var')
    % Choose and load the params.mat file (from the acquisition)
    [params_filename, params_pathname, ~] = uigetfile('*.mat', 'Select the params file', [PDpath, '..\params.mat']);
    load([params_pathname, params_filename])
end

% Load beamforming parameters: PData.mat
% if ~exist('PData', 'var')
%     % Choose and load the RFTimeTagsPerChannel.mat file (from the acquisition and running readRFTimeTags.m)
%     [PData_filename, PData_pathname, ~] = uigetfile('*.mat', 'Select the PData file', [PDpath, '..\PData.mat']);
%     load([PData_pathname, PData_filename])
% end

%% Processing parameter prompt
procParamsPrompt = {'Start file number', 'End file number', 'Use RCA plane or not (1-yes, 0-no)'};
procParamsDefaults = {'1', '', '0'};
procParamsUserInput = inputdlg(procParamsPrompt, 'Input Parameters', 1, procParamsDefaults);

% define # of files manually for now
% str2double(parameterUserInput{});
startFile = str2double(procParamsUserInput{1});
endFile = str2double(procParamsUserInput{2});
numFiles = endFile - startFile + 1;
usingRCAPlane = str2double(procParamsUserInput{3}); % Using a plane from 3D data or not
tau1_index_CBF = 2;
tau2_index_CBF = 5;
tau1_index_CBV = 2;

% if usingRCAPlane
%     [x_mm, y_mm, z_mm] = getReconCoords3D(PData, P);
% 
% else
%     [x_mm, z_mm] = getReconCoords2D(PData, P);
% end
% z_mm = z_mm(voxelRange{1});
% x_mm = x_mm(voxelRange{2});

%% Process tl-fUS
save([PDpath, 'tlfUS_proc_params.mat'], 'usingRCAPlane', 'tau1_index_CBV', 'tau1_index_CBF', 'tau2_index_CBF')
for filenum = startFile:endFile
% for filenum = 1
    disp(filenum)
    tic
    load([PDpath, 'fUSdata-', num2str(filenum), '.mat'], 'g1')

    CBVi = cell(size(g1));
    CBFSi = cell(size(g1));

    for j = 1:3
        [CBFSi{j}, CBVi{j}] = g1_to_CBi(g1{j}, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV);
    end
    save([PDpath, 'tlfUSdata-', num2str(filenum), '.mat'], 'CBVi', 'CBFSi', '-v7.3')
    toc
end
