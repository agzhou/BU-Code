clear all;
defaultpath='E:\0301FlowRBCPhantom';
addpath('E:\0301FlowRBCPhantom');
[FileName,FilePath]=uigetfile(defaultpath);  % read data of a small part of the brain cortex (IQR matrix)  
disp('Loading data...');
load ([FilePath, FileName]);
disp('Data loaded!');
dataRAW=IQ;
clear IQ0 IQData;
[nz, nx, nt]=size(dataRAW);
fileinfo=strsplit(FileName(1:end-4),'-');
Agl=str2num(fileinfo{1});
nAgl=str2num(fileinfo{2});
fCC=str2num(fileinfo{3});
nCC=str2num(fileinfo{4});
iSupFrame=str2num(fileinfo{5});
%% data processing parameters
prompt={'SVD Rank (low):', ['SVD Rank (High):(Max Rank: ',num2str(nCC),')'],'High pass cutoff frequency (Hz)',...
    ['nCC_process (nCC total: ',num2str(nCC),')'], 'Image refine scale','Transducer center frequency (Mhz)','Element Pitch (mm)'};
name='Power Doppler data processing';
% defaultvalue={'15', '110', '40',num2str(nCC),'5', num2str(P.frequency)};
defaultvalue={'15', '110', '40',num2str(nCC),'5', '15.625','0.1'};
numinput=inputdlg(prompt,name, 1, defaultvalue);
RankLow=str2num(numinput{1});
RankHigh=str2num(numinput{2});
DefCutFreq=str2num(numinput{3});
nCC_proc=str2num(numinput{4});
RefScale=str2num(numinput{5});     % image refine scale
fCenter=str2num(numinput{6});     % transducer center frequency
dx=str2num(numinput{7});          % transducer element pitch
dz=1.540/fCenter/2;             % axial samppling step, in mm, the default is 4 samppling points for one wavelength  
%% x,z coordinates after image refined
xCoor=linspace(0,nx*dx,nx*RefScale);
zCoor=linspace(0,nz*dz,nz*RefScale);

%% SVD process 1 (direct SVD use MATLAB)
% [nz,nx,nt]=size(IQR);
% S=reshape(IQR,[nz*nx,nt]);
% [UU,DD,VV]=svd(S);
% DD0=zeros(size(DD));
% DD0(:,21:end)=DD(:,21:end);
% sBlood=reshape(UU*DD0*VV',[nz,nx,nt]);
%% SVD process 2 (eigen-to-SVD use MATLAB)
IQR=dataRAW(:,:,1:nCC_proc); % IQ data used for data processing
[nz,nx,nt]=size(IQR);
rank=[RankLow:RankHigh];
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
Vnoise(:,end)=Vdesc(:,end);
UDelta=S*Vdesc;
sBlood0=reshape(UDelta*Vrank',[nz,nx,nt]);
% sBlood=sBlood0./repmat(std(abs(sBlood0),1,2),[1,nx]);
%%%% Noise equalization 
sNoise=reshape(UDelta*Vnoise',[nz,nx,nt]);
B=ones(30, 30);
sNoiseMed=convn(abs(squeeze(mean(sNoise,3))),B,'same');
sNoiseMedNorm=sNoiseMed/min(sNoiseMed(:));
sBlood=sBlood0./repmat(sNoiseMedNorm,[1,1,nt]);
%% high pass filter
[B,A]=butter(4,DefCutFreq/fCC*2,'high');    %coefficients for the high pass filter
IQR1(:,:,21:20+nCC_proc)=IQR-repmat(IQR(:,:,1),[1,1,nCC_proc]);
for iCC=1:20
    IQR1(:,:,iCC)=IQR1(:,:,41-iCC);
end
sBloodHP=filter(B,A,IQR1,[],3);    % blood signal (filtering in the time dimension)
sBloodHP=sBloodHP(:,:,21:end)./repmat(sNoiseMedNorm,[1,1,nt]);           % the first 4 temporal samples are eliminates (filter oscilations)
%% power doppler
PDISVD=mean(abs(sBlood).^2,3); 
PDISVDdb=10*log10(PDISVD./max(PDISVD(:))); % SVD-based PD image in dB
PDIHP=mean(abs(sBloodHP).^2,3); 
PDIHPdb=10*log10(PDIHP./max(PDIHP(:)));% High pass filter-based PD image in dB
% PDISVDRef=imresize(PDISVDdb,RefScale);
% PDIHPRef=imresize(PDIHPdb,RefScale);
%% figure plot
% figure,
% imagesc(PDIdb);
% caxis([-22 0]);
% colormap(hot);
% colorbar;
fig=figure;
set(fig,'Position',[300 400 900 300]);
subplot(1,2,1);imagesc(xCoor, zCoor, PDISVDdb);
caxis([-25 0]);
colormap(hot);
colorbar;
title ('SVD-based PDI')
xlabel('X [mm]'); ylabel('Z [mm]');
axis equal tight
subplot(1,2,2);imagesc(xCoor, zCoor, PDIHPdb);
caxis([-45 0]);
colormap(hot);
colorbar;
title ('HP-based PDI')
xlabel('X [mm]'); ylabel('Z [mm]');
axis equal tight