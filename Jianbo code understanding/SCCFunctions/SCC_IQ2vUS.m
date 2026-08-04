%% function for processing IQ data to V, vUS
% cluter rejection is based on singular value decomposition (SVD)
function SCC_IQ2vUS(datapath, filename)
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
P.xCoor=P.xCoor(min(ROI(:,2)):max(ROI(:,2)));
P.zCoor=P.zCoor(min(ROI(:,1)):max(ROI(:,1)));
[sIQ, Noise]=SVDfilter(cIQ,PRSinfo.SignalRank); % sIQ: signal IQ
%% color Doppler data processing
[Vcz0]=ColorDoppler(sIQ,P.TWfrequency*1E6,P.CCFR);
Vcz=imresize(Vcz0,PRSinfo.rfnScale);
%% Power Doppler data processing
[PDI0,PDINeg0,PDIPos0]=PowerDoppler(sIQ,Noise);
PDI=imresize(PDI0,PRSinfo.rfnScale);
PDINeg=imresize(PDINeg0,PRSinfo.rfnScale);
PDIPos=imresize(PDIPos0,PRSinfo.rfnScale);
%% vUS data processing
disp('sIQ to GG ...');
GG=IQ2g1(sIQ,PRSinfo.g1startT,PRSinfo.g1nT,PRSinfo.g1nTau);
for iTau=1:PRSinfo.g1nTau
    GGrfn(:,:,iTau)=imresize(squeeze(GG(:,:,iTau)),PRSinfo.rfnScale);
end
disp('vUS Processing ...');
% [Ms, Mf, Vx, Vy, Vz,R,gfit]=US_g1fit_CPX(GGrfn, PRSinfo.Res, PRSinfo.rFrame);
[Ms, Mf, Vx, Vy, Vz,R,gfit]=US_g1fit_INVIVO(GGrfn, PRSinfo.Res, PRSinfo.rFrame);
%% data saving
save([datapath,'vUS',filename(3:end)],'-v7.3','Ms','Mf','Vx','Vy','Vz','R','Vcz','PDI','PDINeg','PDIPos','P');
% save([datapath,'GG',filename(3:end)],'-v7.3','gfit','GG','P');
disp('vUS data saved!')
