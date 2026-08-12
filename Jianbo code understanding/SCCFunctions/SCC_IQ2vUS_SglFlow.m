%% function for processing IQ data to V, fit negative and positive frequency signal separatly
% cluter rejection is based on singular value decomposition (SVD)
function SCC_IQ2vUS_SglFlow(datapath, filename)
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
cIQ=(IQ(min(ROI(:,1)):max(ROI(:,1)),min(ROI(:,2)):max(ROI(:,2)),:));
P.xCoor=interp(P.xCoor(min(ROI(:,2)):max(ROI(:,2))),PRSinfo.rfnScale);
P.zCoor=interp(P.zCoor(min(ROI(:,1)):max(ROI(:,1))),PRSinfo.rfnScale);
[sIQ, sIQHP, sIQHHP, Noise]=IQ2sIQ(cIQ,PRSinfo);
%% vUS data processing
disp('vUS Processing ...');
% g1Info.tStart=PRSinfo.g1startT;
% g1Info.nt=PRSinfo.g1nT;
% g1Info.nTau=PRSinfo.g1nTau;
% [Ms, Mf, Vx, Vy, Vz,Pvz,Vcz,R,GG]=sIQ2vUS_SglFD(sIQ, PRSinfo.FWHM, PRSinfo.rFrame, PRSinfo.rfnScale, PRSinfo.nItpT);
[Mf, Vx, Vz, V, pVz, R, Ms, CR, GGf,sumGG0,sumGGV0,sumGGFV0,sumGGV0_]=sIQ2vUS_SV(sIQ, PRSinfo);
[PDI]=sIQ2PDI(sIQ);
%[PDISVD]=sIQ2PDI(sIQ);
Vcz0=(ColorDoppler(sIQ,PRSinfo));
Vcz = Vcz0.*CR;
%% data saving
pathInfo=strsplit(datapath,'/');
SavePath=['/',strjoin(pathInfo(1:end-2),'/'),'/RESULT-',pathInfo{end-1},'/'];
if ~exist(SavePath)
    mkdir(SavePath);
end
save([SavePath,'vUS',filename(3:end)],'-v7.3','Ms','Mf','Vx','Vz','V','R','pVz','Vcz','PDI','PRSinfo','P','sumGG0','sumGGV0','sumGGFV0','sumGGV0_');
% save([datapath,'GG',filename(3:end)],'-v7.3','gfit','GG','P');
disp('vUS data saved!')
