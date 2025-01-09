%% function for BB pair and track, ULM
% cluter rejection is based on singular value decomposition (SVD)
function SCC_IQ2ULM(datapath, filename)
%% Load data
load([datapath,'ULM-PRSinfo.mat'])
disp('Loading data ...');
load ([datapath, filename]);
disp('Data loaded!');
pathInfo=strsplit(datapath,'/');
fileInfo=strsplit(filename,'-');
%% BB identification
if PRSinfo.BBid==0 % frame-to-frame subtraction based BB identification
    cIQ=diff(IQ,1,3);
    cIQ(:,:,end+1)=IQ(:,:,end)-IQ(:,:,end-2);
    SavePath0=['/',strjoin(pathInfo(1:end-2),'/'),'/RESULTdIQ-',pathInfo{end-1}];
else
    disp('SVD Processing ...');
    [cIQ,Noise]=SVDfilter(IQ,PRSinfo.SignalRank);
    SavePath0=['/',strjoin(pathInfo(1:end-2),'/'),'/RESULTsIQ-',pathInfo{end-1}];
end
%% BB localization
[CoorBB]=dIQ2BB(abs(cIQ), PRSinfo);
[nz0,nx0,nt]=size(IQ);
nz=nz0*PRSinfo.rfn(1);
nx=nx0*PRSinfo.rfn(2);
PRSinfo.Dim=[nz,nx,nt];
Coor.z=[1:PRSinfo.Dim(1)]/PRSinfo.lPix;
Coor.x=[1:PRSinfo.Dim(2)]/PRSinfo.lPix;
PRSinfo.Coor=Coor;
if ~exist(SavePath0)
    mkdir(SavePath0);
end
PRSSinfo=PRSinfo;
SaveName=['CoorBB-',strjoin(fileInfo(2:end),'-')];
save([SavePath0,'/',SaveName],'-V7.3','CoorBB','PRSSinfo')
disp([SaveName, ' is saved']);
clear IQ
%% 2. BB pair and track
% disp(['Pair and track BB...']);
% [BBPD]=BBPT(CoorBB,PRSinfo);
vRange=[0 5; 0 10; 5 40]; % mm/s
NewDt=[10;10;4]*1e-3; % s
PRSinfo.IntPixSize=PRSinfo.lPix;
for iVRange=1:3
    PRSinfo.vRange=vRange(iVRange,:);
    PRSinfo.NewDt=NewDt(iVRange);
    PRSinfo.dCrit=PRSinfo.vRange*PRSinfo.NewDt*1e3; % um
    PRSinfo.nFitp=(PRSinfo.NewDt*PRSinfo.CCFR); % new frame intervel
    PRSinfo.rFrame=PRSinfo.CCFR/PRSinfo.nFitp; % new frame rate
    
    SavePath=[SavePath0,'-V[',num2str(PRSinfo.vRange),']-dt[',num2str(PRSinfo.NewDt*1e3),'ms]/'];
    disp(['Pair and track BB...']);
    CoorBBnF=reshape(CoorBB,[PRSinfo.nFitp,nt/PRSinfo.nFitp]);
    for iTP=1:PRSinfo.nFitp
        [iBBPD]=BBPT_BD(CoorBBnF(iTP,:),PRSinfo);
        BBPD(1+(iTP-1)*(nt/PRSinfo.nFitp-10):iTP*(nt/PRSinfo.nFitp-10),1)=iBBPD;
    end
    %% 3. vULM calculation
    disp(['vULM calculation']);
    [BB,BBV, BBVz]=BBPD2vBB(BBPD,PRSinfo);
    %%
    if ~exist(SavePath)
        mkdir(SavePath);
    end
    PRSSinfo=PRSinfo;
    SaveName=['ULM-',strjoin(fileInfo(2:end),'-')];
    save([SavePath,SaveName],'-V7.3','BB','BBV','BBVz','BBPD','PRSSinfo')
    disp([SaveName, ' processed and saved']);
end

