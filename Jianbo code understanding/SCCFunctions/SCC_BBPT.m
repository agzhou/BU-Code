%% function for BB pair and track, ULM
% cluter rejection is based on singular value decomposition (SVD)
function SCC_BBPT(datapath, filename)
load([datapath,filename]);
load([datapath,'BBPT-PRSinfo.mat'])
%%
[nz,nx,nt]=size(BB);
PRSinfo.Dim=[nz,nx,nt];
Coor.z=[1:PRSinfo.Dim(1)]/PRSinfo.lPix;
Coor.x=[1:PRSinfo.Dim(2)]/PRSinfo.lPix;
PRSinfo.thdTrk=5;
PRSinfo.Coor=Coor;
disp(['Pair and track BB...']);
% [BBPD]=BBPT_CL(BB,PRSinfo,PRSinfo.dCriteria/PRSinfo.lPix,PRSinfo.nTrack);
[BBPD]=BBPT(BB,PRSinfo);
disp(['PDBB-to- down and up flow...']);
BB0=BB;
clear BB;
% [BB,BBV, BBVz]=BBPD2DnUpV(BBPD,PRSinfo,PRSinfo.thdTrk);
[BB,BBV, BBVz]=BBPD2vBB(BBPD,PRSinfo);
%%
pathInfo=strsplit(datapath,'/');
fileInfo=strsplit(filename,'-');
SavePath=['/',strjoin(pathInfo(1:end-2),'/'),'/RESULT-',pathInfo{end-2},'/'];
if ~exist(SavePath)
    mkdir(SavePath);
end
SaveName=['vULM-',strjoin(fileInfo(2:end),'-')];
PRSSinfo=PRSinfo;
save([SavePath,SaveName],'-V7.3','BB','BBV','BBVz','BBPD','BB0','PRSSinfo','Coor')
disp([SaveName, ' processed and saved']);

