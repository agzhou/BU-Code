%% bulk image shift compensaton function
function [xShift, zShift, Xcoef]=ImgSftCmp(IQ, P, nOrgPix,nRfnImg, nRfnCoef)
switch nargin
    case 2
        nOrgPix=0.3;
        nRfnImg=10;
        nRfnCoef=20;
    case 3
        nRfnImg=10;
        nRfnCoef=20;
    case 4
        nRfnCoef=20;
end
[nz,nx,nt]=size(IQ);
IQRfn=imresize(squeeze(IQ(:,:,round(nt/2))),nRfnImg);
maxNpix=max(nOrgPix,1);
RefIQ=IQRfn(nRfnImg*maxNpix+1:end-nRfnImg*maxNpix,nRfnImg*maxNpix+1:end-nRfnImg*maxNpix);
%% calcuate cross correlation of the image/ROI
for it=1:nt
    iIQ0=squeeze(IQ(:,:,it));
    iIQrfn=imresize(iIQ0,nRfnImg);
    for iShiftDx=-nOrgPix*nRfnImg:nOrgPix*nRfnImg
        for iShiftDz=-nOrgPix*nRfnImg:nOrgPix*nRfnImg
            iIQ=iIQrfn(nRfnImg*maxNpix+1+iShiftDz:nRfnImg*maxNpix+iShiftDz+(nz-2*maxNpix)*nRfnImg,nRfnImg*maxNpix+1+iShiftDx:nRfnImg*maxNpix+iShiftDx+(nx-2*maxNpix)*nRfnImg);
            Xcoef0(iShiftDz+nOrgPix*nRfnImg+1,iShiftDx+nOrgPix*nRfnImg+1,it)=sum(iIQ(:).*conj(RefIQ(:)))/(sqrt(sum(abs(iIQ(:)).^2))*sqrt(sum(abs(RefIQ(:)).^2)));
        end
    end
end
% XcoefOrg=squeeze(Xcoef(nOrgPix*nRfnImg+1,nOrgPix*nRfnImg+1,:));
Xcoef=squeeze(max(max(Xcoef0,[],1),[],2));
% plot xCoef map
% iXcoef=abs(squeeze(Xcoef(:,:,1)));
% zCoorCoef=[-nOrgPix*nRfnImg:nOrgPix*nRfnImg]*P.dzImg*1e3/nRfnImg;
% xCoorCoef=[-nOrgPix*nRfnImg:nOrgPix*nRfnImg]*P.dxImg*1e3/nRfnImg;
%         figure,imagesc(zCoorCoef,xCoorCoef,iXcoef)
%         colormap(jet)
%         colorbar
%         xlabel('X shift [um]')
%         ylabel('Z shift [um]')
%         axis equal tight
%% refine cross correlation map and find the image shift
zCoefRfn=linspace(-nOrgPix*nRfnImg,nOrgPix*nRfnImg,(2*nOrgPix*nRfnImg+1)*nRfnCoef);
xCoefRfn=linspace(-nOrgPix*nRfnImg,nOrgPix*nRfnImg,(2*nOrgPix*nRfnImg+1)*nRfnCoef);
for it=1:nt
    iXcoef=abs(squeeze(Xcoef0(:,:,it)));
    iXcoefFilt=imgaussfilt(iXcoef,2);
    iXcoefRfn=imresize(iXcoefFilt,nRfnCoef);
    [zMax(it),xMax(it)]=find(iXcoefRfn==max(iXcoefRfn(:)));
end
zShiftPix=zMax-(nOrgPix*nRfnImg+1/2)*nRfnCoef;
xShiftPix=xMax-(nOrgPix*nRfnImg+1/2)*nRfnCoef;
zShift=zCoefRfn(zMax)*P.dzImg*1e3/nRfnImg;
xShift=xCoefRfn(xMax)*P.dxImg*1e3/nRfnImg;