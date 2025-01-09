%% function for processing IQ data to V, fit negative and positive frequency signal separatly
% cluter rejection is based on singular value decomposition (SVD)
function SCC_IQ2vUS_NP(datapath, filename)
% IQ: IQ data
% PRMT: data processing parameter
% vUS: Obtained vUS results
load([datapath,'vUS-PRSinfo.mat'])
ROI=PRSinfo.ROI;
%% SVD process 1 (direct SVD use MATLAB)
disp('Loading data ...');
load ([datapath, filename]);
disp('Data loaded!');
%%
disp('SVD Processing ...');
cIQ=(IQ(min(ROI(:,1)):max(ROI(:,1)),min(ROI(:,2)):max(ROI(:,2)),1:PRSinfo.g1nT));
P.xCoor=interp(P.xCoor(min(ROI(:,2)):max(ROI(:,2))),PRSinfo.rfnScale);
P.zCoor=interp(P.zCoor(min(ROI(:,1)):max(ROI(:,1))),PRSinfo.rfnScale);
Coor.x=P.xCoor;
Coor.z=P.zCoor;
[sIQ, sIQHP, sIQHHP, eqNoise0]=IQ2sIQ(cIQ,PRSinfo); % 0: no noise equalization
[nz,nx,nt]=size(sIQ);
%% PDI processing
disp('PDI Processing ...');
[PDI0]=sIQ2PDI(sIQ);
[PDIHP0]=sIQ2PDI(sIQHP);
[PDIHHP0]=sIQ2PDI(sIQHHP);
for iD=1:3
    PDI(:,:,iD)=imresize(PDI0(:,:,iD),[nz,nx]*PRSinfo.rfnScale,'bilinear');
    PDIHP(:,:,iD)=imresize(PDIHP0(:,:,iD),[nz,nx]*PRSinfo.rfnScale,'bilinear');
    PDIHHP(:,:,iD)=imresize(PDIHHP0(:,:,iD),[nz,nx]*PRSinfo.rfnScale,'bilinear');
end
eqNoise=imresize(eqNoise0,[nz,nx]*PRSinfo.rfnScale,'bilinear');
clear IQ cIQ sIQHP sIQHHP PDI0 PDIHP0 PDIHHP0
%% vUS data processing
disp('vUS Processing ...');
g1Info.tStart=PRSinfo.g1StartT;
g1Info.nt=PRSinfo.g1nT;
g1Info.nTau=PRSinfo.g1nTau;
[Mf, Vz, V, pVz, Vcz, R, CR, Vx, Ms, pnR,GGf,sumGG0,sumGGV0,sumGGFV0,sumGG0_,GG]=sIQ2vUS_NPDV(sIQ, PRSinfo);
%% data saving
pathInfo=strsplit(datapath,'/');
if PRSinfo.useMsk==1
    SavePath=['/',strjoin(pathInfo(1:end-2),'/'),'/RESULT-',pathInfo{end-1},'-vUSBB/'];
else
    SavePath=['/',strjoin(pathInfo(1:end-2),'/'),'/RESULT-',pathInfo{end-1},'-vUS/'];
end
if ~exist(SavePath)
    mkdir(SavePath);
end

if PRSinfo.useMsk==1
    BB=PRSinfo.BB;
    BBV=PRSinfo.BBV;
    BBVz=PRSinfo.BBVz;
    save([SavePath,'vUSBB',filename(3:end)],'-v7.3','Ms','Mf','Vx','Vz','V','R','CR','pVz','Vcz','PDI','PDIHHP','eqNoise','BB', 'BBV','BBVz','PRSinfo','P');
else
    save([SavePath,'vUS',filename(3:end)],'-v7.3','Ms','Mf','Vx','Vz','V','R','CR','pVz','Vcz','PDI','PDIHP','PDIHHP','eqNoise','PRSinfo','P','sumGG0','sumGGV0','sumGGFV0','sumGG0_','GG');
end

% save([datapath,'GG',filename(3:end)],'-v7.3','gfit','GG','P');
disp('vUS data saved!')
