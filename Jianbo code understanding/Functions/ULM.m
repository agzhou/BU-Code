%% ultrasound localized microscopy data processing function, IQ-based
function [BB]=ULM(dIQ,sysPSF, rfn, PSFtd)
if nargin<4
    PSFtd=0.5;
end
% addpath('D:\OneDrive\Work\PROJ - FUS\CODE\Functions') % Path on JTOPTICS
% sImg(sImg<0)=0;
[nz,nx,nt]=size(dIQ);
%% Non-local means (nlm) spatiotemporal filter
nFperSeg=100;
nSeg=floor(nt/nFperSeg);
dIQ=(dIQ-min(dIQ(:)))/(max((dIQ(:)))-min((dIQ(:))));
for itSeg=1:nSeg
    tic
    tIQ=squeeze(dIQ(:,:,(itSeg-1)*nFperSeg+1:itSeg*nFperSeg));
    Options.kernelratio=2;
    Options.windowratio=2;
    Options.filterstrength=0.008;
    nlmIQ(:,:,(itSeg-1)*nFperSeg+1:itSeg*nFperSeg)=NLMF(tIQ,Options); % NLM filtered Img
    toc
end
%% BB localization
for it=1:nt
    %% 1. it-th frame
    disp(['Processing it=',num2str(it),'...']);
%     fImg=double((squeeze(sImg(:,:,it))));
    fIQ=double((squeeze(nlmIQ(:,:,it))));
    %     figure,imagesc(fImg);colormap hot
%     thd=mean(fImg(:))-1*std(fImg(:));
%     thdImg=fImg;
%     thdImg(thdImg<thd)=thd;%min(thdIQ(:));
    thdIQ=fIQ;
    %     figure,imagesc(thdImg);colormap hot
    %% 2. refine image
    refIQ=imresize(thdIQ,[nz*rfn(1),nx*rfn(2)],'bilinear');
    %     figure,imagesc(refImg);colormap jet
    %% 3. normalized crosscorrelation
    PSFxIQ=normxcorr2(sysPSF,refIQ);
    PSFxIQ(abs(PSFxIQ)<PSFtd)=0;
    %     figure,imagesc(PSFxImg(31:end-30,31:end-30));colormap jet
    %% 4. reject regions with small pixel number
    [rL,rN]=bwlabel(PSFxIQ,4);
    CC=bwconncomp(PSFxIQ,4);
    s=regionprops(CC,'Area');
    nPixRegion=cat(1,s.Area);
    nPixThd=max(median(nPixRegion),rfn(1)*rfn(2)*1);
    for iR=1:rN
        irN=sum(rL(rL==iR))/iR;
        if irN<nPixThd     %refScale^2*3
            rL(rL==iR)=0;
        end
    end
    rL(rL>0)=1;
    PSFxIQ=PSFxIQ.*rL;
%         figure,imagesc(PSFxImg(31:end-30,31:end-30));colormap jet
    %% 5. get the centroid of remainning regions
    BBcent=(imregionalmax(PSFxIQ));
    BBcent(BBcent.*PSFxIQ<PSFtd)=0;
    %     figure,imagesc(BBcent);colormap(hot);caxis([0 2])
    BB(:,:,it)=uint8(BBcent(31:end-30,31:end-30));
end
    
    