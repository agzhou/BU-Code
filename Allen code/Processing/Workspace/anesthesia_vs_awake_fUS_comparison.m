% Look at the artifacts for awake vs. anesthetized fUS measurements

%% load data (manually for now)
% datapath = 'E:\Allen BME-BOAS-27 Data Backup\AZ03 Stroke RC15gV\fUS\05-06-2025 pre-stroke\fUS processing results 25 to 150 SVs\';
% startFile = 1;
% endFile = 285;

datapath = 'D:\Allen\Data\AZ01 fUS RCA\06-05-2025 awake manual right whisker stim\run 1 5 trials wooden stick right whiskers 11 angles -5 to 5 deg 2500 Hz\fUS results SVs 30 to 180\';
startFile = 1;
endFile = 144;
%% RCA
load([datapath, 'tlfUSdata-', num2str(startFile)])
CBViallSF = zeros([size(CBVi), endFile - startFile + 1]);
CBViallSF(:, :, :, startFile) = CBVi;
for filenum = startFile + 1:endFile
% for filenum = 2
    load([datapath, 'tlfUSdata-', num2str(filenum)])
    CBViallSF(:, :, :, filenum) = CBVi;


end

% % Testing
figure; imagesc(squeeze(max(CBViallSF(:, :, :, 1), [], 1) .^ 0.5)'); colormap hot
figure; imagesc(squeeze(max(CBViallSF(:, :, :, 2), [], 1) .^ 0.5)'); colormap hot

% generateTiffStack_acrossframes(CBViallSF, [8.8, 8.8, 8], 'hot', 1:80)