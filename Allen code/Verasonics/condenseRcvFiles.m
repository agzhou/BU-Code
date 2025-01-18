% Take a folder of RcvData files (each a group of frames) and put them into
% one matrix

datapath = 'G:\Allen\Data\01-17-2025 AZ001 ULM\L22-14v\run 1 allen code left eye\';

load([datapath, 'params.mat']) % load parameters
probe_name = P.Trans.name;
filename_structure = [probe_name, '-RcvData-', num2str(P.maxAngle), '-', num2str(P.na), '-', num2str(P.PRF), '-', num2str(P.frameRate), '-', num2str(P.numFramesPerBuffer), '-'];

r = []; % Create empty variable to append the data to

%% Stack RcvData files into one variable r
numFiles = 315; % regexp to get last file number?

for filenum = 1:numFiles
% for filenum = 1:2
    load([datapath, filename_structure, num2str(filenum), '.mat'])
    r = cat(1, r, RcvData); % concatenate along the time/z sample dimension

end

%% Recon
numFiles = 315; % regexp to get last file number?
P.numSubFrames = P.numFramesPerBuffer; % temporary fix
IQstack = []; % Create empty variable to append the data to

for filenum = 1:numFiles
% for filenum = 1:2
    load([datapath, filename_structure, num2str(filenum), '.mat'])
    IQ = LA_DAS(RcvData, P, P.wl/2);
    disp(strcat("recon on buffer ", num2str(filenum), " done"))
    IQstack = cat(1, IQstack, IQ); % concatenate along the time/z sample dimension
end