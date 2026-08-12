%% function for processing IQ data to obtained Power Doppler Image (PDI)
% cluter rejection is based on high pass filtering
function IQHP=IQFilt_HP(IQ,PRMT)
% IQ: IQ data
% PRMT: data processing parameter
% PDI: Obtained Power Doppler Image
RankLow=PRMT.RankLow;
RankHigh=PRMT.RankHigh;
DefCutFreq=PRMT.DefCutFreq;
fCC=PRMT.fCC;
nCC_proc=PRMT.nCC_proc;
RefScale=PRMT.RefScale;     % image refine scale
fCenter=PRMT.fCenter;     % transducer center frequency
dx=PRMT.dx;          % transducer element pitch
dz=1.540/fCenter/2;             % axial samppling step, in mm, the default is 4 samppling points for one wavelength  
%% x,z coordinates after image refined
[nz,nx,nt]=size(IQ);
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
IQR=IQ(:,:,1:nCC_proc); % IQ data used for data processing
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
UDelta=S*Vdesc;
%%%% Noise equalization 
Vnoise=zeros(size(Vdesc));
Vnoise(:,end)=Vdesc(:,end);
sNoise=reshape(UDelta*Vnoise',[nz,nx,nt]);
sNoiseMed=medfilt2(abs(squeeze(mean(sNoise,3))),[30 30],'symmetric');
sNoiseMedNorm=sNoiseMed/min(sNoiseMed(:));
%% high pass filter - based Power Doppler
[B,A]=butter(4,DefCutFreq/fCC*2,'high');    %coefficients for the high pass filter
IQR1(:,:,21:20+nCC_proc)=IQR-repmat(IQR(:,:,1),[1,1,nCC_proc]);
for iCC=1:20
    IQR1(:,:,iCC)=IQR1(:,:,41-iCC);
end
sBloodHP=filter(B,A,IQR1,[],3);    % blood signal (filtering in the time dimension)
IQHP=sBloodHP(:,:,21:end)./repmat(sNoiseMedNorm,[1,1,nt]);           % the first 20 temporal samples are eliminates (filter oscilations)
