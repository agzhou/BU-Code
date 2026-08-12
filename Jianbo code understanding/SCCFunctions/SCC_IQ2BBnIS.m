%% function for processing IQ data to BB and IS, ULM
% cluter rejection is based on singular value decomposition (SVD)
function SCC_IQ2BBnIS(datapath, filename)
% IQ: IQ data
% PRMT: data processing parameter
% vUS: Obtained vUS results
load([datapath,'BBnIS-PRSinfo.mat'])
ROI=PRSinfo.ROI;

%% bulk image shift calculation
% nPxROI=floor(abs(diff(ROI(:,2)))/PRSinfo.ISnSubROI(2)); % number of Xpixel for each ROI
% nPzROI=floor(abs(diff(ROI(:,1)))/PRSinfo.ISnSubROI(1)); % number of Zpixel for each ROI
% RefPlane=squeeze(Bulk(:,:,round(nt/2))); % in block reference plane
% for ixROI=1:PRSinfo.ISnSubROI(2)
%     tic
%     for izROI=1:PRSinfo.ISnSubROI(1)
%         iIQ=Bulk((izROI-1)*nPzROI+1:izROI*nPzROI,(ixROI-1)*nPxROI+1:ixROI*nPxROI,:);
%         [xShift(:,izROI,ixROI), zShift(:,izROI,ixROI), Xcoef(:,izROI,ixROI)]=ImgSftCmp(iIQ, PRSinfo,PRSinfo.ISmaxOrgPix,PRSinfo.ISnImgRfn,PRSinfo.ISnxCoefRfn);
%     end
%     toc
% end
% %% inter block image shift calculation
% BlockRefPlane(:,:,1)=PRSinfo.RefBulkBlock(min(ROI(:,1)):max(ROI(:,1)),min(ROI(:,2)):max(ROI(:,2)));
% BlockRefPlane(:,:,2)=RefPlane;
% for ixROI=1:PRSinfo.ISnSubROI(2)
%     tic
%     for izROI=1:PRSinfo.ISnSubROI(1)
%         iIQBlock=BlockRefPlane((izROI-1)*nPzROI+1:izROI*nPzROI,(ixROI-1)*nPxROI+1:ixROI*nPxROI,:);
%         [xShiftBLK(izROI,ixROI), zShiftBLK(izROI,ixROI), XcoefBLK(izROI,ixROI)]=ImgSftCmp(iIQBlock, 2,5,20);
%     end
%     toc
% end
% ShiftBLK(izROI,ixROI,1): xShiftBLK; ShiftBLK(izROI,ixROI,2): zShiftBLK
%%
P.PRSinfo_IQ2BB=PRSinfo;
P.lPix=PRSinfo.lPix;
SavePath=[datapath, '/A-BBIS/'];
if ~exist(SavePath)
    mkdir(SavePath);
end
SaveName=['BBIS',filename(3:end)];
% save([SavePath,SaveName],'-v7.3','BB','xShift','zShift','Xcoef','ROI','P','RefPlane','PRSinfo','xShiftBLK','zShiftBLK','XcoefBLK');
save([SavePath,SaveName],'-v7.3','BB','P');
disp([SaveName]);

