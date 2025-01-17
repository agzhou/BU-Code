% Take a folder of RcvData files (each a group of frames) and put them into
% one matrix

datapath = 'G:\Allen\Data\';

load('params.mat') % load parameters
filename_structure = [num2str(P.maxAngle), '-', num2str(P.na), '-', num2str(P.PRF), '-', num2str(P.frameRate), '-', num2str(P.numFramesPerBuffer), '-'];

r = []; % Create empty variable to append the data to


% regexp to get files?
for file = 1:numFiles
    load([filename_structure, num2str(P.bufferIndex), '.mat'])
    r = cat(1, r, RcvData); % concatenate along the time/z sample dimension

end