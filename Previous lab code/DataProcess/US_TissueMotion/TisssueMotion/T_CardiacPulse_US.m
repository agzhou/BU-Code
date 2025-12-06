% cardiac signal and phase
addpath('D:\OneDrive\Work\PROJ - FUS\CODE\Functions') % Path on JTOPTICS
load('H:\PJ - USI\PROJ - fUS Velocimetry\0611InvivoTest\Cardiac&Respiratory\10-5-200-1000-1-CP3-1-IQ.mat');
% cIQ=hilbert(squeeze(sum(IQ,3))); % Coherence Compounding IQ data, then Hilbert transform to get complext Coherence compounded IQ data (envelop+phase)
cIQ=IQ;
[nz,nx,nt]=size(IQ);
fRate=P.CCFR; % CC frame rate
P.zCoor=linspace(P.zCoor(1),P.zCoor(1)+P.zCoor(end),length(P.zCoor)); 
% P.xCoor=P.xCoor/2;
%% 1. calculate PDI
prompt={'SVD Rank (low):', ['SVD Rank (High):(Max Rank: ',num2str(P.numCCframes),')'],...
    ['nCC_process (nCC total: ',num2str(P.numCCframes),')']};
name='Power Doppler data processing';
defaultvalue={'10', num2str(P.numCCframes),num2str(P.numCCframes)};
numinput=inputdlg(prompt,name, 1, defaultvalue);
handles.RankLow=str2num(numinput{1});
handles.RankHigh=str2num(numinput{2});
nCC_proc=str2num(numinput{3});
handles.DefnCC_proc=nCC_proc;
IQR=IQ(:,:,1:nCC_proc);
% 1.1 SVD process (eigen-to-SVD use MATLAB)
rank=[handles.RankLow:handles.RankHigh];
[nz,nx,nt]=size(IQR);
S=reshape(IQR,[nz*nx,nt]);
S_COVt=(S'*S);
[V,D]=eig(S_COVt); % V is the right singular Vector of S/eigenvector; D is the eigenvalue/square of Singular value
for it=1:nt 
    Ddiag(it)=abs(sqrt(D(it,it)));
end
Ddiag=20*log10(Ddiag/max(Ddiag)); % singular value in db
[Ddesc, Idesc]=sort(Ddiag,'descend');
% figure,plot(Ddesc);
for it=1:nt
    Vdesc(:,it)=V(:,Idesc(it));
end
Vrank=zeros(size(Vdesc));
Vrank(:,rank)=Vdesc(:,rank);
Vnoise=zeros(size(Vdesc));
Vnoise(:,end)=Vdesc(:,min(300,nt));
UDelta=S*Vdesc;
sBlood0=reshape(UDelta*Vrank',[nz,nx,nt]);
% sBlood=sBlood0./repmat(std(abs(sBlood0),1,2),[1,nx]);
%%%% Noise equalization 
sNoise=reshape(UDelta*Vnoise',[nz,nx,nt]);
sNoiseMed=medfilt2(abs(squeeze(mean(sNoise,3))),[30 30],'symmetric');
sNoiseMedNorm=sNoiseMed/min(sNoiseMed(:));
sBlood=sBlood0./repmat(sNoiseMedNorm,[1,1,nt]);
% 1.2 power doppler
PDI=mean(abs(sBlood).^2,3); 
PDIdb=10*log10(PDI./max(PDI(:)));
% 1.3 Plot PDI
figure,
imagesc(P.xCoor,P.zCoor,PDIdb);
axis tight equal
title(['SVD-based PDI, Rank=[',num2str([handles.RankLow handles.RankHigh]),']'])
xlabel('x [mm]');
ylabel('z [mm]');
colormap hot
%% 2. Cardiac and respiratory map from IQ
[nz,nx,nt]=size(IQ);
tCoor=linspace(1/P.CCFR,nt/P.CCFR,nt)*1e3;
% 2.1 Fourier transfom of g1(iTau)

Lfft=2^nextpow2(nt);
FphaseIQ=(fft(squeeze(angle(IQ(:,:,:))-mean(angle(IQ(:,:,:)),3)),Lfft,3));
FrealIQ=(fft(squeeze(real(IQ(:,:,:))-mean(real(IQ(:,:,:)),3)),Lfft,3));
FabsIQ=(fft(squeeze(abs(IQ(:,:,:))-mean(abs(IQ(:,:,:)),3)),Lfft,3));
fCoor=linspace(0, fRate/2,Lfft/2);
fCoorStep=mean(diff(fCoor));
% 2.2 locate the cardiac rate and respiratory rate
figure,
subplot(3,1,1);
plot(fCoor,abs(squeeze(mean(mean(FphaseIQ(:,:,1:Lfft/2),1),2))))
xlabel('Frequency [HZ]')
ylabel('Power (phase)')
subplot(3,1,2);
plot(fCoor,abs(squeeze(mean(mean(FrealIQ(:,:,1:Lfft/2),1),2))))
xlabel('Frequency [HZ]')
ylabel('Power (real)')
subplot(3,1,3);
plot(fCoor,abs(squeeze(mean(mean(FabsIQ(:,:,1:Lfft/2),1),2))))
xlabel('Frequency [HZ]')
ylabel('Power (magnitude)')

prompt={'Cardiac Frequency', 'Cardiac Frequency Range (CF+/- CFR)',...
    'Respiratory Frequency', 'Respiratory Frequency Range (RF+/- RFR)'};
name='CF and RF processing';
defaultvalue={'8','0','1.2','0'};
numinput=inputdlg(prompt,name, 1, defaultvalue);
CF=str2num(numinput{1});
CFR=str2num(numinput{2});
RF=str2num(numinput{3});
RFR=str2num(numinput{4});
% 2.3 get cardiac and respiratory maps 
CFrange=round((CF-CFR)/fCoorStep):round((CF+CFR)/fCoorStep);
Cardiac=squeeze(max(abs(FabsIQ(:,:,CFrange)),[],3));
RFrange=round((RF-RFR)/fCoorStep):round((RF+RFR)/fCoorStep);
Resp=squeeze(max(abs(FrealIQ(:,:,RFrange)),[],3));
% 2.4 plot cardiac and respiratory maps 
fig=figure;
set(fig,'Position',[200 400 1200 300])
subplot(1,2,1);imagesc(P.xCoor,P.zCoor,Cardiac);
axis tight equal
title(['Cardiac map - IQ'])
xlabel('x [mm]');
ylabel('z [mm]');
subplot(1,2,2);imagesc(P.xCoor,P.zCoor,Resp);
axis tight equal
title(['Respiratory map - IQ'])
xlabel('x [mm]');
ylabel('z [mm]');

%% 2. Cardiac and respiratory map from G1
% GG=IQ2g1(IQ,1,nt,800);
[nz,nx,nTau]=size(GG);
tCoor=linspace(1/P.CCFR,nTau/P.CCFR,nTau)*1e3;
% 2.1 Fourier transfom of g1(iTau)

Lfft=2^nextpow2(nTau);
FphaseGG=(fft(squeeze(angle(GG(:,:,:))-mean(angle(GG(:,:,:)),3)),Lfft,3));
FrealGG=(fft(squeeze(real(GG(:,:,:))-mean(real(GG(:,:,:)),3)),Lfft,3));
FabsGG=(fft(squeeze(abs(GG(:,:,:))-mean(abs(GG(:,:,:)),3)),Lfft,3));
fCoor=linspace(0, fRate/2,Lfft/2);
fCoorStep=mean(diff(fCoor));
% 2.2 locate the cardiac rate and respiratory rate
figure,
subplot(3,1,1);
plot(fCoor,abs(squeeze(mean(mean(FphaseGG(:,:,1:Lfft/2),1),2))))
xlabel('Frequency [HZ]')
ylabel('Power (phase)')
subplot(3,1,2);
plot(fCoor,abs(squeeze(mean(mean(FrealGG(:,:,1:Lfft/2),1),2))))
xlabel('Frequency [HZ]')
ylabel('Power (real)')
subplot(3,1,3);
plot(fCoor,abs(squeeze(mean(mean(FabsGG(:,:,1:Lfft/2),1),2))))
xlabel('Frequency [HZ]')
ylabel('Power (magnitude)')

prompt={'Cardiac Frequency', 'Cardiac Frequency Range (CF+/- CFR)',...
    'Respiratory Frequency', 'Respiratory Frequency Range (RF+/- RFR)'};
name='CF and RF processing';
defaultvalue={'8','0','1.2','0'};
numinput=inputdlg(prompt,name, 1, defaultvalue);
CF=str2num(numinput{1});
CFR=str2num(numinput{2});
RF=str2num(numinput{3});
RFR=str2num(numinput{4});
% 2.3 get cardiac and respiratory maps 
CFrange=ceil((CF-CFR)/fCoorStep):ceil((CF+CFR)/fCoorStep);
Cardiac=squeeze(max(abs(FabsGG(:,:,CFrange)),[],3));
RFrange=ceil((RF-RFR)/fCoorStep):ceil((RF+RFR)/fCoorStep);
Resp=squeeze(max(abs(FphaseGG(:,:,RFrange)),[],3));
% 2.4 plot cardiac and respiratory maps 
fig=figure;
set(fig,'Position',[200 400 1200 300])
subplot(1,2,1);imagesc(P.xCoor,P.zCoor,Cardiac);
axis tight equal
title(['Cardiac map - GG'])
xlabel('x [mm]');
ylabel('z [mm]');
subplot(1,2,2);imagesc(P.xCoor,P.zCoor,Resp);
axis tight equal
title(['Respiratory map - GG'])
xlabel('x [mm]');
ylabel('z [mm]');
%% 3. plot selected point IQ and GG
figure,imagesc((squeeze(PDIdb)));
colormap hot
nSlt=1;
[xSlt,zSlt]=ginput(nSlt);
xSlt=round(xSlt);
zSlt=round(zSlt);
clear IQslt GGslt FphaseGGslt FabsGGslt FphaseIQslt FabsIQslt
for iSlt=1:nSlt
    IQslt (:,iSlt)= squeeze(IQ(zSlt(iSlt),xSlt(iSlt),:));
    GGslt(:,iSlt)= squeeze(GG(zSlt(iSlt),xSlt(iSlt),:));
    FphaseGGslt(:,iSlt)= squeeze(FphaseGG(zSlt(iSlt),xSlt(iSlt),:));
    FabsGGslt(:,iSlt)= squeeze(FabsGG(zSlt(iSlt),xSlt(iSlt),:));
    FphaseIQslt(:,iSlt)= squeeze(FphaseIQ(zSlt(iSlt),xSlt(iSlt),:));
    FabsIQslt(:,iSlt)= squeeze(FabsIQ(zSlt(iSlt),xSlt(iSlt),:));
end
[nz,nx,nt]=size(IQ);
tCoorIQ=linspace(1/fRate,nt/fRate,nt);
tCoorGG=linspace(1/fRate,nTau/fRate,nTau);
fig=figure;
set(fig,'Position',[400 400 800 500])
subplot(2,2,1); plot(tCoorIQ,squeeze((angle(IQslt(:,1:nSlt)))))
xlabel('t [s]');
ylabel(['phase(IQ)'])
subplot(2,2,2); plot(tCoorIQ,squeeze(abs(IQslt(:,1:nSlt))))
xlabel('t [s]');
ylabel(['abs(IQ)'])

subplot(2,2,3); plot(tCoorGG,squeeze((angle(GGslt(:,1:nSlt)))))
xlabel('t [s]');
ylabel(['phase(GG)'])
subplot(2,2,4); plot(tCoorGG,squeeze(abs(GGslt(:,1:nSlt))))
xlabel('t [s]');
ylabel(['abs(GG)'])

fig=figure;
set(fig,'Position',[400 400 800 500])
subplot(2,2,1);
plot(fCoor,abs(squeeze(FphaseIQslt(1:Lfft/2,1:nSlt))))
xlabel('Frequency [HZ]')
ylabel('Power (phase IQ)')
xlim([0 20])
subplot(2,2,2);
plot(fCoor,abs(squeeze(FabsIQslt(1:Lfft/2,1:nSlt))))
xlabel('Frequency [HZ]')
ylabel('Power (magnitude IQ)')
xlim([0 20])
subplot(2,2,3);
plot(fCoor,abs(squeeze(FphaseGGslt(1:Lfft/2,1:nSlt))))
xlabel('Frequency [HZ]')
ylabel('Power (phase GG)')
xlim([0 20])
subplot(2,2,4);
plot(fCoor,abs(squeeze(FabsGGslt(1:Lfft/2,1:nSlt))))
xlabel('Frequency [HZ]')
ylabel('Power (magnitude GG)')
xlim([0 20])
%% 4. pulse wave
% nt_g1=50; % time interval: ~6 ms
% nt=4;   % calculate fourth time lag g1
% iTau=3;
% GG=zeros(nz,nx,nt,nt-nt_g1-nt);
% for it=1:nt-nt_g1-nt
%     for iz=1:nz
%         for ix=1:nx
%             GG (iz,ix,:,it)= squeeze(IQ2g1(IQ(iz,ix,:), it, nt_g1, nt));
%         end
%     end
%     it
% end
% iTauGG=squeeze(GG(:,:,iTau,:));
% 
% %% plot time couse of selected time lag
% tCoor=linspace(1/fRate,(nt-nt_g1-nt)/fRate,nt-nt_g1-nt);
% figure,
% subplot(2,1,1); plot(tCoor,squeeze(unwrap(angle(GG(iTau,1:nSlt,:)))))
% xlabel('t [s]');
% ylabel(['phase(g_1(',num2str(iTau/BRate*1e3),' ms))'])
% title(['P1x=',num2str(xSlt(1)*L_xStep),'um, P2x=',num2str(xSlt(2)*L_xStep), 'um'])
% subplot(2,1,2); plot(tCoor,squeeze(abs(GG(iTau,1:nSlt,:))))
% xlabel('t [s]');
% ylabel(['abs(g_1(',num2str(iTau/BRate*1e3),' ms))'])
% title(['P1x=',num2str(xSlt(1)*L_xStep),'um, P2x=',num2str(xSlt(2)*L_xStep), 'um'])
% %% Fourier transfom of g1(iTau)
% Lfft=2^nextpow2(nt-nt_g1-nt);
% FphaseGG=fftshift(fft(squeeze(angle(GG(iTau,1:nSlt,:))-mean(angle(GG(iTau,1:nSlt,:)),3)),Lfft,2));
% FrealGG=fftshift(fft(squeeze(abs(GG(iTau,1:nSlt,:))-mean(abs(GG(iTau,1:nSlt,:)),3)),Lfft,2));
% fCoor=linspace(-BRate/2, BRate/2,Lfft);
% figure,
% subplot(2,1,1); plot(fCoor,squeeze(abs(FphaseGG)))
% xlabel('f [HZ]');
% ylabel(['Power (angle)'])
% title(['P1x=',num2str(xSlt(1)*L_xStep),'um, P2x=',num2str(xSlt(2)*L_xStep), 'um'])
% xlim([0 30])
% subplot(2,1,2); plot(fCoor,squeeze(abs(FrealGG)))
% xlabel('t [s]');
% ylabel(['Power (magnitude)'])
% title(['P1x=',num2str(xSlt(1)*L_xStep),'um, P2x=',num2str(xSlt(2)*L_xStep), 'um'])
% xlim([0 30])
% %% cross correlation of g1(iTau)
% % iSlt=2;
% % [r, lags,~]=crosscorr(squeeze(angle(GG(iTau,1,:))),squeeze(angle(GG(iTau,iSlt,:))));
% % % r=r/(mean(angle(GG(iTau,1,:)).^2)*mean(angle(GG(iTau,2,:)).^2));
% % figure,plot(lags*1e3/BRate, r)
% % xlabel('Time lag [ms]');
% % ylabel('Cross correlation coefficient')
% % xlim([-200 200])
% % title(['P1x=',num2str(xSlt(1)*L_xStep),'um, P2x=',num2str(xSlt(iSlt)*L_xStep), 'um'])
% % plot multiple
% figure
% Pcolor=jet(nSlt);
% 
% for ix=1:nSlt
%     [r, lags,~]=crosscorr(squeeze(angle(GG(iTau,cSlt,:))),squeeze(angle(GG(iTau,ix,:))));
%     nLags=numel(lags);
%     [~, Imax(ix)]=max(r(round(nLags/2)-round(BRate*30e-3)+1:round(nLags/2)+round(BRate*30e-3))); % find the peak within [-30 30] ms time lag range
%     Imax(ix)=Imax(ix)+round(nLags/2)-round(BRate*30e-3);
%     Tmax(ix)=lags(Imax(ix))*1e3/BRate;
%     hold on;
%     % r=r/(mean(angle(GG(iTau,1,:)).^2)*mean(angle(GG(iTau,2,:)).^2));
%     plot(lags*1e3/BRate, r,'color',Pcolor(ix,:))
%     xlabel('Time lag [ms]');
%     ylabel('Cross correlation coefficient')
%     xlim([-200 200])
%     IndText=round(numel(lags/2));
% %     text(double(lags(IndText)),double(r(IndText)),num2str(xSlt(iSlt)*L_xStep));
% 
% end
% % title(['P1x=',num2str(xSlt(1)*L_xStep),'um, P2x=',num2str(xSlt*L_xStep), 'um'])
% title(['PxC=',num2str(xSlt(cSlt)*L_xStep),'um, PzC=',num2str(zSlt(1)*1.5),'um'])
% legend({num2str((xSlt-xSlt(cSlt))'*L_xStep)})
% grid on
% axis tight
% figure,plot((xSlt(1:nSlt)-xSlt(cSlt))*L_xStep,Tmax(1:nSlt))
% hold on, plot(0,Tmax(cSlt),'ro')
% ylim([-30 30])
% xlabel('x distance from selected source point [um]');
% ylabel('tDelay [ms]')
% title(['z=',num2str(zSlt(1)*1.5),' um'])