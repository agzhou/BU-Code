%% function for processing IQ data to V, vUS
% cluter rejection is based on singular value decomposition (SVD)
function SCC_IQ2GGAVG2vUS(datapath, filename)
% IQ: IQ data
% PRMT: data processing parameter
% vUS: Obtained vUS results
load([datapath,'vUS-PRSinfo.mat'])
ROI=PRSinfo.ROI;
%% SVD process 1 (direct SVD use MATLAB)
fileInfo=strsplit(filename,'-');
for iRpt=1:PRSinfo.nRPT
    ifilename=[strjoin(fileInfo(1:end-1),'-'),'-',num2str(PRSinfo.startRPT+iRpt-1)];
    disp('Loading data ...');
    load ([datapath, ifilename]);
    disp('Data loaded!');
    %%
    disp('SVD Processing ...');
    cIQ=(IQ(min(ROI(:,1)):max(ROI(:,1)),min(ROI(:,2)):max(ROI(:,2)),:));
    [sIQ, Noise]=SVDfilter(cIQ,PRSinfo.SignalRank); % sIQ: signal IQ
    %% color Doppler data processing
    [Vcz0(:,:,iRpt)]=ColorDoppler(sIQ,P.TWfrequency*1E6,P.CCFR);
    %% Power Doppler data processing
    [PDI0(:,:,iRpt),PDINeg0(:,:,iRpt),PDIPos0(:,:,iRpt)]=PowerDoppler(sIQ,Noise);
    %% vUS data processing
    disp('sIQ to GG ...');
    GG0(:,:,:,iRpt)=IQ2g1(sIQ,PRSinfo.g1startT,PRSinfo.g1nT,PRSinfo.g1nTau);
end
P.xCoor=P.xCoor(min(ROI(:,2)):max(ROI(:,2)));
P.zCoor=P.zCoor(min(ROI(:,1)):max(ROI(:,1)));
GG=mean(GG0,4);
Vcz=imresize(mean(Vcz0,3),PRSinfo.GGrfnSale);
PDI=imresize(mean(PDI0,3),PRSinfo.GGrfnSale);
PDINeg=imresize(mean(PDINeg0,3),PRSinfo.GGrfnSale);
PDIPos=imresize(mean(PDIPos0,3),PRSinfo.GGrfnSale);
clear GG0 Vcz0 PDI0 PDINeg0 PDIPos0
for iTau=1:PRSinfo.g1nTau
    GGrfn(:,:,iTau)=imresize(squeeze(GG(:,:,iTau)),PRSinfo.GGrfnSale);
end
disp('vUS Processing ...');
[Ms, Mf, Vx, Vy, Vz,R,gfit]=US_g1fit_CPX(GGrfn, PRSinfo.Res, PRSinfo.rFrame);
% [Ms, Mf, Vx, Vy, Vz,R,gfit]=US_g1fit_INVIVO(GGrfn, PRSinfo.Res, PRSinfo.rFrame);
%% data saving
save([datapath,'vUS',filename(3:end)],'-v7.3','Ms','Mf','Vx','Vy','Vz','R','Vcz','PDI','PDINeg','PDIPos','P');
% save([datapath,'GG',filename(3:end)],'-v7.3','gfit','GG','P');
disp('vUS data saved!')
