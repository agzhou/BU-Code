%% for multiple files PDI process
clear all;
load('D:\CODE\Mains\DAQParameters.mat');
defaultpath=DAQInfo.savepath;
addpath('D:\CODE\Functions');
[FileName,FilePath]=uigetfile(defaultpath);  % read data of a small part of the brain cortex (IQR matrix)  
fileInfo=strsplit(FileName(1:end-4),'-');
startCP0=fileInfo{7}(3:end);
if isempty(startCP0)
    startCP0=0;
end
%% Load IQ data DAQ information
myFile=matfile([FilePath,FileName]);
P=myFile.P;
%% parameters
PRMT.fCC=P.CCFR;
nCC=P.numCCframes;
%% data processing parameters
prompt={'SVD Rank (low):', ['SVD Rank (High):(Max Rank: ',num2str(nCC),')'],'High pass cutoff frequency (Hz)',...
    ['nCC_process (nCC total: ',num2str(nCC),')'], 'Image refine scale','Transducer center frequency (Mhz)','Element Pitch (mm)','Cluter Rejection (SVD:0; HP: 1; None: 2)',...
    'Start file (CP)','Number of files (CPs)','Start Repeat', 'Number of Repeats'};
name='Power Doppler data processing';
defaultvalue={'25', num2str(nCC), '50',...
    num2str(nCC),'1', '18','0.1','0',...
    num2str(startCP0), '1',fileInfo{8},'1'};
numinput=inputdlg(prompt,name, 1, defaultvalue);
SVDRank=[str2num(numinput{1}), str2num(numinput{2})];
HPcut=str2num(numinput{3});
PRMT.nCC_proc=str2num(numinput{4});
PRMT.RefScale=str2num(numinput{5});     % image refine scale
PRMT.fCenter=str2num(numinput{6});     % transducer center frequency
PRMT.dx=str2num(numinput{7});          % transducer element pitch
PRMT.Method=str2num(numinput{8});     % transducer center frequency
startCP=str2num(numinput{9});
PRMT.nCP=str2num(numinput{10});          % number of coronal planes
startRpt=str2num(numinput{11});
PRMT.nRpt=str2num(numinput{12});          % number of repeat for each coronal plane

P.SVDRank=SVDRank;
P.HPcut=HPcut;
xCoor=P.xCoor;
zCoor=P.zCoor;
for iCP=startCP:startCP+PRMT.nCP-1
    for iRpt=startRpt:startRpt+PRMT.nRpt-1
        %% Load iFile
        iFileInfo=fileInfo;
        if PRMT.nCP>1
            iFileInfo{7}(3:end)=num2str(iCP);
        end            
        iFileInfo{8}=num2str(iRpt);
        iFileName=[strjoin(iFileInfo,'-'),'.mat'];
        iSaveName=['PDI-',strjoin(iFileInfo(2:end),'-'),'.mat'];
        disp(['Loading and processing: ', iFileName,' ...']);
        load([FilePath,iFileName]);
        [sIQ, sIQHP, Noise]=IQ2sIQ(IQ,SVDRank,HPcut,P.CCFR);
        [PDI]=sIQ2PDI(sIQHP);
        % Save data
        save([FilePath,iSaveName],'-V7.3','PDI','P');
        disp([iFileName,' is processed and saved!']);
    end
end