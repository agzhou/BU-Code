%% singular value decomposition filter function
function [Signal, Noise, Bulk]=SVDfilter(Data,SignalRank,cBLK)
if nargin<3
    cBLK=0;
else
    cBLK=1;
end
[nz,nx,nSVD]=size(Data);
S=reshape(Data,[nz*nx,nSVD]);
S_COVt=(S'*S);
[V,D]=eig(S_COVt); % V is the right singular Vector of S/eigenvector; D is the eigenvalue/square of Singular value
for it=1:nSVD 
    Ddiag(it)=abs(sqrt(D(it,it)));
end
Ddiag=20*log10(Ddiag/max(Ddiag)); % singular value in db
[Ddesc, Idesc]=sort(Ddiag,'descend');
% figure,plot(Ddesc);
for it=1:nSVD
    Vdesc(:,it)=V(:,Idesc(it));
end
UDelta=S*Vdesc;
%% SVD filtered
Vrank=zeros(size(Vdesc));
rank=SignalRank(1):SignalRank(2);
Vrank(:,rank)=Vdesc(:,rank);
Signal=reshape(UDelta*Vrank',[nz,nx,nSVD]);
%% Noise 
Vnoise=zeros(size(Vdesc));
Vnoise(:,min(400,nSVD):end)=Vdesc(:,min(400,nSVD):end);
sNoise=reshape(UDelta*Vnoise',[nz,nx,nSVD]);
sNoiseMed=medfilt2(abs(squeeze(mean(sNoise,3))),[50 50],'symmetric');
Noise=sNoiseMed/min(sNoiseMed(:));
if cBLK==1
    %% Bulk motion, tissue
    Vblk=zeros(size(Vdesc));
    Vblk(:,1:3)=Vdesc(:,1:3);
    Bulk=reshape(UDelta*Vblk',[nz,nx,nSVD]);
else 
    Bulk=0;
end


