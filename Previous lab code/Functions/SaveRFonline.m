%% save RF data
function SaveRFonline(RFRAW)
tic
disp(datestr(now,'dd-mm-yyyy HH:MM:SS FFF'));
addpath('D:\CODE\Functions');
load('D:\CODE\Mains\DAQParameters.mat');
load('D:\CODE\Mains\P.mat');
RFdataFilename = ['RF-',num2str(DAQInfo.CCangle),'-',num2str(DAQInfo.numAngles),'-',num2str(DAQInfo.CCFR),'-',num2str(DAQInfo.numCCframes),'-',num2str(DAQInfo.numSupFrames),'-',DAQInfo.Filename, '-','1'];
savepath=DAQInfo.savepath;
DAQInfo.Plane=1;
if exist([savepath,RFdataFilename,'.mat'])
    DAQInfo.Plane=1;
    while (exist([savepath,RFdataFilename,'.mat'])==2)
        DAQInfo.Plane=DAQInfo.Plane+1;
        RFdataFilename = ['RF-',num2str(DAQInfo.CCangle),'-',num2str(DAQInfo.numAngles),'-',num2str(DAQInfo.CCFR),'-',num2str(DAQInfo.numCCframes),'-',num2str(DAQInfo.numSupFrames),'-',DAQInfo.Filename,'-',num2str(DAQInfo.Plane)];
    end
end
fileID = fopen([savepath,'DAQ-INFO.txt'],'at');
TNow=datestr(now,'dd-mm-yyyy HH:MM:SS FFF');
fmt = [RFdataFilename, ':' TNow,'\n\n'];
fprintf(fileID,fmt);
fclose(fileID);
disp(['Saving SupFrame: ',num2str(DAQInfo.Plane)])
savefast ([savepath,RFdataFilename], 'RFRAW','P')
tSave=toc
disp ('RF DATA SAVED!')
% pause(P.tIntPDI/1e6-P.numCCframes/P.CCFR-tSave);
pause(P.tIntPDI/1e6-tSave);
