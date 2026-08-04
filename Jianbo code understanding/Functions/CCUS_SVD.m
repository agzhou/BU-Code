%% Singular Value Decomposition of CCUS data
function [Signal_SVD]=CCUS_SVD(DATA,Rej_SVD)
[nz,nx,nt]=size(DATA);
S=reshape(DATA, [nz*nx,nt]); % reshape 3D stack data to spatical-temporal 2D data (Casorati matrix)
rank=[20:180];
%% direct SVD using matlab function
[UU,DD,VV]=svd(S);
for it=1:nt
    DDdiag(it)=abs((DD(it,it)));
end
DDdiag=20*log10(DDdiag/max(DDdiag)); % singular value in db
[DDdesc0, IIdesc]=sort(DDdiag,'descend');
figure,plot(DDdesc0);
for it=1:nt
    VVdesc(:,it)=VV(:,IIdesc(it));
    DDdesc(:,it)=DD(:,IIdesc(it));
end
DDrank=zeros(size(DDdesc));
DDrank(:,rank)=DDdesc(:,rank);
ss=reshape(UU*DDrank*VVdesc',[nz,nx,nt]);
fig=figure;set(fig,'Position',[500 300 1100 300]);
subplot(1,3,1);imagesc(log(squeeze(abs(ss(:,:,3))))); colormap(jet); 
subplot(1,3,2);plot((squeeze(ss(11,14,:))));
subplot(1,3,3);plot(squeeze(abs(ss(11,14,:))));
%% eigen-based calculation for SVD
S_COVt=(S'*S);
[V,D]=eig(S_COVt); % V is the right singular Vector of S/eigenvector; D is the eigenvalue/square of Singular value
for it=1:nt 
    Ddiag(it)=abs(sqrt(D(it,it)));
end
Ddiag=20*log10(Ddiag/max(Ddiag)); % singular value in db
[Ddesc, Idesc]=sort(Ddiag,'descend');
figure,plot(Ddesc);
for it=1:nt
    Vdesc(:,it)=V(:,Idesc(it));
end
Vrank=zeros(size(Vdesc));
Vrank(:,rank)=Vdesc(:,rank);
UDelta=S*Vdesc;
s=reshape(UDelta*Vrank',[nz,nx,nt]);
fig=figure;set(fig,'Position',[500 300 1100 300]);
subplot(1,3,1);imagesc(log(squeeze(abs(s(:,:,3))))); colormap(jet); 
subplot(1,3,2);plot((squeeze(s(11,14,:))));
subplot(1,3,3);plot(squeeze(abs(s(11,14,:))));
%% depth depended SVD
for iz=1:nz
    Siz=squeeze(DATA(iz,:,:)); 
    [UUiz,DDiz,VViz]=svd(Siz);
    DDiz0=zeros(size(DDiz));
    DDiz0(:,21:end)=DDiz(:,21:end);
    sBlood(iz,:,:)=(UUiz*DDiz0*VViz'); 
    DDzdiag(iz,:)=(max(DDiz,[],1));
end
figure,plot(20*log10((DDzdiag([100 150 200 250],:)/max(DDzdiag(:))))')

