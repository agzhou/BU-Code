%% ultrasound localized microscopy data processing function, IQ-based
% PRSSinfo.sysPSF
% PRSSinfo.PSFtd
% PRSSinfo.rfn
function [CoorBB]=dIQ2BB(dIQ,PRSSinfo)
% addpath('D:\OneDrive\Work\PROJ - FUS\CODE\Functions') % Path on JTOPTICS
% sImg(sImg<0)=0;
[nz,nx,nt]=size(dIQ);
%% Non-local means (nlm) spatiotemporal filter
% nFperSeg=100;
% nSeg=floor(nt/nFperSeg);
% dIQ=(dIQ-min(dIQ(:)))/(max((dIQ(:)))-min((dIQ(:))));
% for itSeg=1:nSeg
%     tic
%     tIQ=squeeze(dIQ(:,:,(itSeg-1)*nFperSeg+1:itSeg*nFperSeg));
%     Options.kernelratio=2;
%     Options.windowratio=2;
%     Options.filterstrength=0.008;
%     nlmIQ(:,:,(itSeg-1)*nFperSeg+1:itSeg*nFperSeg)=NLMF(tIQ,Options); % NLM filtered Img
%     toc
% end
%% BB localization
nTperChk=50; 
nChk=nt/nTperChk;
for iChk=1:nChk
    %% 1. it-th frame
    disp(['Processing iChk=',num2str(iChk),'...']);
%     fIQ=double((squeeze(nlmIQ(:,:,it))));
    fIQ=double(reshape(squeeze(dIQ(:,:,(iChk-1)*nTperChk+1:iChk*nTperChk)),[nz,nx*nTperChk]));
    %     figure,imagesc(thdImg);colormap hot
    %% 2. refine image
    refIQ=int16(imresize(fIQ/max(fIQ(:))*256*5,[nz*PRSSinfo.rfn(1),nx*nTperChk*PRSSinfo.rfn(2)],'bilinear'));
    clear fIQ;
    %     figure,imagesc(refImg);colormap jet
    %% 3. normalized crosscorrelation
    PSFxIQ=gather((normxcorr2(PRSSinfo.sysPSF,gpuArray(refIQ))));
    PSFxIQ(abs(PSFxIQ)<PRSSinfo.PSFtd)=0;
    clear refIQ;
    %     figure,imagesc(PSFxIQ(31:end-30,31:end-30));colormap jet
    %% 4. reject regions with small pixel number
    [rL,rN]=bwlabel(PSFxIQ,4);
    nPixThd=PRSSinfo.rfn(1)*PRSSinfo.rfn(2)*3;
    sBW=bwareaopen(rL,nPixThd);
    PSFxIQ=PSFxIQ(31:end-30,31:end-30).*sBW(31:end-30,31:end-30);
%     PSFxIQ=PSFxIQ(31:end-30,31:end-30);
    %% 5. get the centroid coordinates of identified BBs
    CC=bwconncomp(PSFxIQ,4);
    propROI=regionprops(CC,PSFxIQ,'Area','WeightedCentroid');
    peakCoor=cat(1,propROI.WeightedCentroid);
    for iframe=1:nTperChk
        index=logical((peakCoor(:,1)>(iframe-1)*nx*PRSSinfo.rfn(2)).*(peakCoor(:,1)<(iframe*nx*PRSSinfo.rfn(2))));
        CoorBB{(iChk-1)*nTperChk+iframe}=[peakCoor(index,2),peakCoor(index,1)-(iframe-1)*nx*PRSSinfo.rfn(1)]; % [z,x]
    end
    %% plot identified BB
%     startShow=10;
%     nShow=10;
%     ColorSet=num2cell(jet(nShow), 2);
%     figure;
%     for itShow=startShow:startShow+nShow-1
%         hold on,
%         plot(CoorBB{itShow}(:,2),nz*PRSSinfo.rfn(2)+1-CoorBB{itShow}(:,1),'.','Color',ColorSet{itShow-startShow+1});
%     end
    %% Form a 2D image
%     BBcnt=zeros(size(PSFxIQ));
%     Ind=sub2ind(size(PSFxIQ),round(peakCoor(:,2)),round(peakCoor(:,1)));
%     BBcnt(Ind)=1;
%     %     figure,imagesc(conv2(BBcnt(:,1:51200),PRSSinfo.sysPSF,'same'));colormap(hot);caxis([0 1.5])
%     BB(:,:,(iChk-1)*nTperChk+1:iChk*nTperChk)=uint8(reshape(BBcnt,[nz*PRSSinfo.rfn(1),nx*PRSSinfo.rfn(2),nTperChk]));
end
    
    