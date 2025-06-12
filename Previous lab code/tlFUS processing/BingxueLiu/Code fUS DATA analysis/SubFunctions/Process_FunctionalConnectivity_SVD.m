%% FC analysis using SVD
[nz,nx,nSVD]=size(cData);
rData=reshape(cData,[nz*nx,nSVD]);
Sn=rData./mean(rData,2); % temporal normalization
Sm=Sn-mean(Sn,2);
S=Sm;
% S=rData./max(rData,[],2);
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
U=reshape(UDelta,[nz,nx,nt]);
%%
[VzCmap, VzCmapDn,VzCmapUp]=Colormaps_fUS;
for r=1:12
    figure,imagesc(U(:,:,r));colormap(VzCmap); 
    caxis([-1 1]); colorbar; axis equal tight;
    title(['r=',num2str(r)])
end
