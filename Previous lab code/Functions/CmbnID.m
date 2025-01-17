%% combine ID data
function CmbnID (P0)

% disp(datestr(now,'dd-mm-yyyy HH:MM:SS FFF'));
addpath('D:\CODE\Functions');
load('D:\CODE\Mains\DAQParameters.mat');
Filename = ['ID',DAQInfo.FileNameFull(3:end-4)];
savepath=[DAQInfo.savepath,Filename,'\'];

for iFrame=1:P.numCCframes
    iFile=[Filename,'-',num2str(iFrame),'.mat'];
    load([savepath,iFile]);
    ID(:,:,iFrame)=iID;
end
clear P
P=P0;
% fileID = fopen([savepath,'DAQ-INFO.txt'],'at');
% TNow=datestr(now,'dd-mm-yyyy HH:MM:SS FFF');
% fmt = [IDFilename, ':' TNow,'\n\n'];
% fprintf(fileID,fmt);
% fclose(fileID);
savefast ([DAQInfo.savepath,Filename], 'ID','P')
% disp ('RF DATA SAVED!')
% pause(P.tIntPDI/1e6-P.numCCframes/P.CCFR-tSave);
rmdir(savepath,'s')