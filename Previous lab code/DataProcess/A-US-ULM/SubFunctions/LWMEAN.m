%% 2D local weighted averaging
function DATAout=LWMEAN(DATAin, lW)
% lW, length of local window, odd
X=[1:lW]-(ceil(lW/2)); Y=X;
for iX=1:lW
    for iY=1:lW
        W(iX,iY)=sqrt(X(iX)^2+Y(iY)^2)+1;
    end
end
W=1./W;
% DATAout=convn(DATAin,W,'same')/sum(W(:));
[nx,ny]=size(DATAin);
for ix=1:nx
    for iy=1:ny
        
        lcWD=DATAin(max(ix-floor(lW/2),1):min(ix+floor(lW/2),nx),max(iy-floor(lW/2),1):min(iy+floor(lW/2),ny)).*...
            W(max(ix-floor(lW/2),1)-ix+ceil(lW/2):min(ix+floor(lW/2),nx)-ix+ceil(lW/2),max(iy-floor(lW/2),1)-iy+ceil(lW/2):min(iy+floor(lW/2),ny)-iy+ceil(lW/2));
        DATAout(ix,iy)=sum(lcWD(:));
    end
end