% addpath('./subFunctions')
clearvars
close all

RFPath = uigetdir('C:\');
outPath = uigetdir('F:\');
%%
load([RFPath,'/params.mat'],'P')
RFDir = dir([RFPath '/*.rf']);
%%
for filenum = 1:length(RFDir) %RF is zero indexed lol
    filenum
    % tic
    RFData = loadBinaryRF([RFPath,'/', P.fsFileName, sprintf('%05d', filenum+256), '.rf'], P);
    RFData = unstackFrames(RFData, P);
    RFData = reshape(RFData, [], 2*P.na, P.Resource.Parameters.numRcvChannels, P.numFramesPerBuffer); % Reshape to dimensions (samples, angles, channels, frames)
    RcvData = RFData(:, :, P.Trans.Connector, :);
    filenameParams = {'RF',num2str(P.maxAngle),num2str(P.na),num2str(P.frameRate)...
        ,num2str(P.numFramesPerBuffer),'1',num2str(filenum+257)};
    filenameOut = [strjoin(filenameParams,'-') '.mat'];
    
    save(fullfile(outPath,filenameOut),'RcvData','-v7.3')
    % toc
end
