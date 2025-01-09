%% function for processing IQ data to obtained Power Doppler Image (PDI)
% cluter rejection is based on singular value decomposition (SVD)
function SCC_IQ2PDI(datapath, filename)
% IQ: IQ data
% PRSinfo: data processing parameter
% PDI: Obtained Power Doppler Image

%% Load data
disp('Loading data...')
load([datapath,'PDI-PRSinfo.mat'])
load([datapath,filename]);
pathInfo=strsplit(datapath,'/');
P.PRSinfo=PRSinfo;
Coor.x=P.xCoor;
Coor.z=P.zCoor;
cIQ=IQ(:,:,1:PRSinfo.nCC_proc); % IQ data used for data processing
clear IQ
%% clutter rejection
disp('SVD & HP...')
[sIQ, sIQHP, sIQHHP, eqNoise]=IQ2sIQ(cIQ,PRSinfo); % 0: no noise equalization
% clear cIQ
disp('sIQ to PDI...')
[PDI]=sIQ2PDI(sIQ);
[PDIHP]=sIQ2PDI(sIQHP);
[PDIHHP]=sIQ2PDI(sIQHHP);


%% g1-fUS
[nz0,nx0,nt]=size(sIQ);
% I.0 all frequency signal and SNR
PRSSinfo.g1StartT=1;
PRSSinfo.g1nTau=100;
PRSSinfo.g1nT=nt;
[GG0, Numer, Denom] = sIQ2GG(sIQHP, PRSSinfo); % g1 of whole frequency signal
temp_deno=(conj(cIQ(:,:,PRSSinfo.g1StartT:PRSSinfo.g1StartT-1+PRSSinfo.g1nT))).*(cIQ(:,:,PRSSinfo.g1StartT:PRSSinfo.g1StartT-1+PRSSinfo.g1nT));
IQDenom=repmat(mean(temp_deno,3),[1,1,PRSSinfo.g1nTau]); 
GIQ = Numer./IQDenom;

% fIQ = fft(sIQHP,nt,3);
% for iNP=1:2
%     iFIQ=zeros(size(fIQ));
%     iFIQ(:,:,(floor(nt/2))*(iNP-1)+1:floor(nt/2)*iNP)=fIQ(:,:,(floor(nt/2))*(iNP-1)+1:floor(nt/2)*iNP);
%     % III. 3 iGG and ggCR
%     iIQ=(ifft(iFIQ,nt,3));
%     GGnp(:,:,:,iNP)=sIQ2GG(iIQ, PRSSinfo);
% end
% GGnp = single(GGnp);

Ag1=abs(GG0(:,:,1))-min(abs(GG0(:,:,10)),[],3);
SUMg1=sum(abs(GG0(:,:,1:10)),3);
SumIg1=sum(abs(imag(GG0(:,:,1:10))),3);

%% SAVE PDI dat %%%%%%%%%%%%%%%%
fileInfo=strsplit(filename,'-');
NamePDI=['PDI-',strjoin(fileInfo(2:end),'-')];
% save([datapath,NamePDI],'PDISVD','PDINeg','PDIPos','PDIHP','sNoiseMedNorm','P');
SavePath=['/',strjoin(pathInfo(1:end-2),'/'),'/RESULT-',pathInfo{end-1},'-PDI/'];
% if ~exist(SavePath)
    mkdir(SavePath);
% end
save([SavePath,NamePDI],'-v7.3','PDI','PDIHP','PDIHHP','eqNoise','Ag1','SUMg1','SumIg1','Coor','P','GG0','Numer','Denom','GIQ','IQDenom');
disp('PDI data saved!')
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [nz,nx,nt]=size(IQR);
% rank=[RankLow:RankHigh];
% S=reshape(IQR,[nz*nx,nt]);
% S_COVt=(S'*S);
% [V,D]=eig(S_COVt); % V is the right singular Vector of S/eigenvector; D is the eigenvalue/square of Singular value
% for it=1:nt 
%     Ddiag(it)=abs(sqrt(D(it,it)));
% end
% Ddiag=20*log10(Ddiag/max(Ddiag)); % singular value in db
% [Ddesc, Idesc]=sort(Ddiag,'descend');
% % figure,plot(Ddesc);
% for it=1:nt
%     Vdesc(:,it)=V(:,Idesc(it));
% end
% UDelta=S*Vdesc;
% %%%% Noise equalization 
% Vnoise=zeros(size(Vdesc));
% Vnoise(:,end)=Vdesc(:,end);
% sNoise=reshape(UDelta*Vnoise',[nz,nx,nt]);
% sNoiseMed=medfilt2(abs(squeeze(mean(sNoise,3))),[30 30],'symmetric');
% sNoiseMedNorm=sNoiseMed/min(sNoiseMed(:));
% %% SVD-based Power Doppler after noise equalization
% Vrank=zeros(size(Vdesc));
% Vrank(:,rank)=Vdesc(:,rank);
% sBlood0=reshape(UDelta*Vrank',[nz,nx,nt]);
% sBlood=sBlood0./repmat(sNoiseMedNorm,[1,1,nt]);
% % sBloodavg=MoveAvg(abs(sBlood),10,3);
% % PDISVD=mean((abs(sBlood)./sBloodavg+median(sBloodavg,3)).^2,3); % new 
% PDISVD=mean(abs(sBlood).^2,3);  % Original method
% % PDISVDdb=10*log10(PDISVD./max(PDISVD(:))); % SVD-based PD image in dB
% % PDISVDdb=log10(PDISVD); % SVD-based PD image log
% % PDI=imresize(PDISVDdb,RefScale);
% %% high pass filter - based Power Doppler
% [B,A]=butter(4,DefCutFreq/fCC*2,'high');    %coefficients for the high pass filter
% IQR1(:,:,21:20+nCC_proc)=IQR-repmat(IQR(:,:,1),[1,1,nCC_proc]);
% for iCC=1:20
%     IQR1(:,:,iCC)=IQR1(:,:,41-iCC);
% end
% sBloodHP=filter(B,A,IQR1,[],3);    % blood signal (filtering in the time dimension)
% sBloodHP=sBloodHP(:,:,21:end)./repmat(sNoiseMedNorm,[1,1,nt]);           % the first 4 temporal samples are eliminates (filter oscilations)
% % % sBloodHPavg=MoveAvg(abs(sBloodHP),10,3);
% % % PDIHP=mean((abs(sBloodHP)./sBloodHPavg+median(sBloodHPavg,3)).^2,3); % new
% PDIHP=mean(abs(sBloodHP).^2,3); % original
% % % PDIHPdb=10*log10(PDIHP./max(PDIHP(:)));% High pass filter-based PD image in dB
% % % PDI=imresize(PDIHPdb,RefScale);
% %% directional Power Doppler data processing %%%%%%%%%%%%%%%%%%%%%%%%%%
% Fs=P.CCFR;
% [nz,nx,nt]=size(sBlood);
% nf= 2^nextpow2(2*nt);           % Fourier transform points
% P.fCoor=linspace(-Fs/2,Fs/2,nf);  % frequency coordinate
% 
% SpecBlood=fftshift(fft(sBlood,nf,3),3);
% PDINeg=squeeze(sum(abs(SpecBlood(:,:,1:floor(nf/2)-1)).^2,3));
% % PDINegdb=10*log10(PDINeg/max(PDINeg(:)));
% % PDINegRef=imresize(PDINegdb,RefScale);
% PDIPos=squeeze(sum(abs(SpecBlood(:,:,floor(nf/2)+1:nf)).^2,3));
% % PDIPosdb=10*log10(PDIPos/max(PDIPos(:)));
% % PDIPosRef=imresize(PDIPosdb,RefScale);
% %% SAVE PDI dat %%%%%%%%%%%%%%%%
% fileInfo=strsplit(filename,'-');
% NamePDI=['PDI-',strjoin(fileInfo(2:end),'-')];
% % save([datapath,NamePDI],'PDISVD','PDINeg','PDIPos','PDIHP','sNoiseMedNorm','P');
% save([datapath,NamePDI],'PDISVD','PDIHP','PDINeg','PDIPos','sNoiseMedNorm','P');
% disp('PDI data saved!')
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