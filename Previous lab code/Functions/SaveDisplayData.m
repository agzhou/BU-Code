%% save RF data
function SaveDisplayData(iID)
tic
% disp(datestr(now,'dd-mm-yyyy HH:MM:SS FFF'));
addpath('D:\CODE\Functions');
load('D:\CODE\Mains\DAQParameters.mat');
Filename = ['ID',DAQInfo.FileNameFull(3:end-4)];
savepath=[DAQInfo.savepath,Filename,'\'];
if exist(savepath)~= 7
    mkdir (savepath)
end
Frame=1;
IDFilename = [Filename,'-',num2str(Frame)];
if exist([savepath,IDFilename,'.mat'])
    Frame=1;
    while (exist([savepath,IDFilename,'.mat'])==2)
        Frame=Frame+1;
        IDFilename = [Filename,'-',num2str(Frame)];
    end
end
% fileID = fopen([savepath,'DAQ-INFO.txt'],'at');
% TNow=datestr(now,'dd-mm-yyyy HH:MM:SS FFF');
% fmt = [IDFilename, ':' TNow,'\n\n'];
% fprintf(fileID,fmt);
% fclose(fileID);
savefast ([savepath,IDFilename], 'iID')
% disp ('RF DATA SAVED!')
% pause(P.tIntPDI/1e6-P.numCCframes/P.CCFR-tSave);
