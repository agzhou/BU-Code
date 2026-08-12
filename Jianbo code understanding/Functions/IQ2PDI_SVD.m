%% function for processing IQ data to obtained Power Doppler Image (PDI)
% cluter rejection is based on singular value decomposition (SVD)
function PDI=IQ2PDI_SVD(IQ,PRSinfo)
% IQ: IQ data
% PRMT: data processing parameter
% PDI: Obtained Power Doppler Image
% DefCutFreq=PRMT.DefCutFreq;
% fCC=PRMT.fCC;
%% SVD process 
IQR=IQ(:,:,1:PRSinfo.nCC_proc); % IQ data used for data processing
[nz,nx,nt]=size(IQR);
rank=[PRSinfo.SignalRank(1):PRSinfo.SignalRank(2)];
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
UDelta=S*Vdesc;
%%%% Noise equalization 
Vnoise=zeros(size(Vdesc));
Vnoise(:,end)=Vdesc(:,end);
sNoise=reshape(UDelta*Vnoise',[nz,nx,nt]);
sNoiseMed=medfilt2(abs(squeeze(mean(sNoise,3))),[30 30],'symmetric');
sNoiseMedNorm=sNoiseMed/min(sNoiseMed(:));
%% SVD-based Power Doppler after noise equalization
Vrank=zeros(size(Vdesc));
Vrank(:,rank)=Vdesc(:,rank);
sBlood0=reshape(UDelta*Vrank',[nz,nx,nt]);
sBlood=sBlood0./repmat(sNoiseMedNorm,[1,1,nt]);
PDISVD=mean(abs(sBlood).^2,3); 
% PDISVDdb=10*log10(PDISVD./max(PDISVD(:))); % SVD-based PD image in dB
% PDISVDdb=log10(PDISVD); % SVD-based PD image log
% PDI=imresize(PDISVDdb,RefScale);
PDI=PDISVD;
%% high pass filter - based Power Doppler
% [B,A]=butter(4,DefCutFreq/fCC*2,'high');    %coefficients for the high pass filter
% IQR1(:,:,21:20+nCC_proc)=IQR-repmat(IQR(:,:,1),[1,1,nCC_proc]);
% for iCC=1:20
%     IQR1(:,:,iCC)=IQR1(:,:,41-iCC);
% end
% sBloodHP=filter(B,A,IQR1,[],3);    % blood signal (filtering in the time dimension)
% sBloodHP=sBloodHP(:,:,21:end)./repmat(sNoiseMedNorm,[1,1,nt]);           % the first 4 temporal samples are eliminates (filter oscilations)
% PDIHP=mean(abs(sBloodHP).^2,3); 
% PDIHPdb=10*log10(PDIHP./max(PDIHP(:)));% High pass filter-based PD image in dB
% PDI=imresize(PDIHPdb,RefScale);
% %% figure plot
% fig=figure;
% % set(fig,'Position',[300 400 900 300]);
% % subplot(1,2,1);
% imagesc(xCoor, zCoor, PDI);
% caxis([-25 0]);
% colormap(hot);
% colorbar;
% title ('SVD-based PDI')
% xlabel('X [mm]'); ylabel('Z [mm]');
% axis equal tight
% % subplot(1,2,2);
% % imagesc(xCoor, zCoor, PDI);
% % caxis([-45 0]);
% % colormap(hot);
% % colorbar;
% % title ('HP-based PDI')
% % xlabel('X [mm]'); ylabel('Z [mm]');
% % axis equal tight