%% save PA RF data
function SavePARFonline(RFRAW)
tic
disp(datestr(now,'dd-mm-yyyy  HH:MM:SS FFF'));
addpath('D:\CODE\Functions');
load('D:\CODE\Mains\PAT\PATAQParameters_L22_14v.mat');
load('D:\CODE\Mains\PAT\P.mat');
RFdataFilename = ['RF-',PATAQInfo.Filename, '-','1'];
savepath=PATAQInfo.savepath;
PATAQInfo.Plane=1;
if exist([savepath,RFdataFilename,'.mat'])
    PATAQInfo.Plane=1;
    while (exist([savepath,RFdataFilename,'.mat'])==2)
        PATAQInfo.Plane=PATAQInfo.Plane+1;
        RFdataFilename = ['RF-',PATAQInfo.Filename,'-',num2str(PATAQInfo.Plane)];
    end
end
fileID = fopen([savepath,'PAT-DAQ-INFO.txt'],'at');
TNow=datestr(now,'dd-mm-yyyy HH:MM:SS FFF');
fmt = [RFdataFilename, ':' TNow,'\n\n'];
fprintf(fileID,fmt);
fclose(fileID);
disp(['Saving Frame: ',num2str(PATAQInfo.Plane)])
savefast([savepath,RFdataFilename],'RFRAW','P')
tSave=toc
disp ('RF DATA SAVED!')
% pause(P.tIntPDI/1e6-P.numCCframes/P.CCFR-tSave);
% pause(1/20-tSave); %20Hz pulse