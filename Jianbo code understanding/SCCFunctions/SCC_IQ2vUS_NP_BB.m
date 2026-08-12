%% function for processing IQ data to V, fit negative and positive frequency signal separatly
% cluter rejection is based on singular value decomposition (SVD)
function SCC_IQ2vUS_NP_BB(datapath, filename)
% IQ: IQ data
% PRMT: data processing parameter
% vUS: Obtained vUS results
load([datapath,'vUS-PRSinfo.mat'])
ROI=PRSinfo.ROI;
%% SVD process 1 (direct SVD use MATLAB)
disp('Loading data ...');
load([datapath,'BB.mat']);
load ([datapath, filename]);
disp('Data loaded!');
%%
disp('SVD Processing ...');
cIQ=(IQ(min(ROI(:,1)):max(ROI(:,1)),min(ROI(:,2)):max(ROI(:,2)),1:PRSinfo.g1nT));
P.xCoor=interp(P.xCoor(min(ROI(:,2)):max(ROI(:,2))),PRSinfo.rfnScale);
P.zCoor=interp(P.zCoor(min(ROI(:,1)):max(ROI(:,1))),PRSinfo.rfnScale);
Coor.x=P.xCoor;
Coor.z=P.zCoor;
[sIQ, sIQHP, sIQHHP, eqNoise]=IQ2sIQ(cIQ,PRSinfo.SignalRank,25,P.CCFR,0);
%% PDI processing
disp('PDI Processing ...');
[PDI]=sIQ2PDI(sIQ);
[PDIHP]=sIQ2PDI(sIQHP);
[PDIHHP]=sIQ2PDI(sIQHHP);
clear IQ cIQ sIQ sIQHHP
%% vUS data processing
disp('vUS Processing ...');
g1Info.tStart=PRSinfo.g1startT;
g1Info.nt=PRSinfo.g1nT;
g1Info.nTau=PRSinfo.g1nTau;
for iBB=1:3
    FitMsk(:,:,iBB)=imresize(abs(BB(:,:,iBB)),PRSinfo.rfnScale/5);
end
[Ms, Mf, Vx, Vy, Vz,Vz0, Pvz,Vcz,R,GG]=sIQ2vUS_NP(sIQHP./eqNoise, FitMsk, g1Info, PRSinfo.Res, PRSinfo.rFrame,PRSinfo.rfnScale,PRSinfo.nItpT,1);
clear FitMsk
%% data saving
pathInfo=strsplit(datapath,'/');
SavePath=['/',strjoin(pathInfo(1:end-2),'/'),'/RESULT-',pathInfo{end-1},'-vUSBB/'];
if ~exist(SavePath)
    mkdir(SavePath);
end
save([SavePath,'vUS',filename(3:end)],'-v7.3','Ms','Mf','Vx','Vy','Vz','Vz0','R','Vcz','Coor','PDI','PDIHP','PDIHHP','eqNoise','BB','BBV','P','PBB');
disp('vUS data saved!')
