%% microbubble pair and track - based on closest criteria only
function [BBPD]=BBPT_CL(BB,P,Thd,nTrkF)
% BB, bubble position obtained from Img2ULM
% Thd, pair distance criteria, in pixel. Thd(1): maximum searching distance
% between two frames; Thd(2): maximum footprint pairing distance
% BBPT, Paired and Tracked bubble 
% Revised from Xiaojun's code by JT, 08162018
if nargin<3
    Thd=[8 2]; % 
    nTrkF=5; % number of time frames to track
elseif nargin<4
    nTrkF=5;
end

BBp=squeeze(BB(:,:,:)); % Positive bubble corresponding to bubble identified at current frame
[nz,nx,nt]=size(BBp); % nubmer of x,y,t points
%% obtain all bubbles' coordinates
for it=1:nt
    [zBBp,xBBp]=find(BBp(:,:,it)==1);
    CoorBBp{it}=[zBBp,xBBp];
end
BBPD=cell(nt-nTrkF,1);
%% 
for it=1:nt-nTrkF
    itCoorBBp=CoorBBp{it};  % coordinates (z, x) of BB in frame it
    nBBp=size(itCoorBBp,1); % number of total positive bubbles
    ipdBB=1; % number of paired BB in frame it
    %% Find out BB satisfying the criteria
    for iBB=1:nBBp
        trackable=0; % check if iBB is trackable
        for iTrkF=1:nTrkF % try to pair and track the same bubble in the next nTrkF frames
            if iTrkF==1
                CoorCBBp=itCoorBBp(iBB,:); % Coordinate of current frame positive bubble of iBB
            end
            %% pairing frame
            CoorNBBp=CoorBBp{it+iTrkF}; % Coordinates of next frame positive bubble
            
            dBBp=sqrt(sum(((CoorCBBp-CoorNBBp).^2),2)); % distance between current positive bubble (iBB) and next frame positive bubble
            sPositionInd=find(dBBp==0);
            dBBp(sPositionInd)=[];
            CoorNBBp(sPositionInd,:)=[];
            if ~isempty(dBBp)
                if min(dBBp)<Thd(1)  % criteria for bubble pairing
                    indPD=find(dBBp==min(dBBp));
                    CoorNPD=CoorNBBp(indPD(1),:); % new coordinate of the paired BB in the new frame
                    dircZ=CoorNPD(1)-CoorCBBp(1);
                    dircX=CoorNPD(2)-CoorCBBp(2);
                    if dircZ>=0 && dircX>=0
                        dirc=1;
                    elseif dircZ>=0 && dircX<0
                        dirc=2;
                    elseif dircZ<0 && dircX<0
                        dirc=3;
                    elseif dircZ<0 && dircX>=0
                        dirc=4;
                    end
                    v=min(dBBp)*P.lPix/1e3*P.CCFR;
                    vz=(CoorNPD(1)-CoorCBBp(1))*P.lPix/1e3*P.CCFR;
                    iBBPD(:,iTrkF,ipdBB)=[CoorCBBp(1),CoorCBBp(2),CoorNPD(1),CoorNPD(2),iTrkF,v,vz,dirc]'; % paired BB, [zCoor,xCoor,zCoorPair,xCoorPair,nTrackable,speed,direction]
                    CoorCBBp=CoorNPD; % update the coordiante of iBB
                    trackable=1;
                else
                    break;
                end
            else
                break;
            end
        end
        if trackable==1
            ipdBB=ipdBB+1;
        end
    end
    if ipdBB>1
        BBPD{it}=iBBPD;
        clear iBBPD
    end
    
end

% %% plot result
% [VzCmap1, VzCmap2, ActCmap]=Colormaps_fUS;
% BB0=zeros(nz,nx);
% BBV=zeros(nz,nx);
% for it=1:nt-nTrkF
%     inPDBB=size(PDBB{it},3);
%     for iBB=1:inPDBB
%         nTrk=max(PDBB{it}(5,:,iBB));
%         v=sqrt((PDBB{it}(1,1,iBB)-PDBB{it}(3,nTrk,iBB))^2+(PDBB{it}(2,1,iBB)-PDBB{it}(4,nTrk,iBB))^2)*P.lPix/1e3*P.CCFR; % Average speed
%         vDir=sign(PDBB{it}(1,1,iBB)-PDBB{it}(3,nTrk,iBB));
%         BB0(PDBB{it}(1,1,iBB),PDBB{it}(2,1,iBB))=BB0(PDBB{it}(1,1,iBB),PDBB{it}(2,1,iBB))+1;
%         BB0(PDBB{it}(3,nTrk,iBB),PDBB{it}(4,nTrk,iBB))=BB0(PDBB{it}(3,nTrk,iBB),PDBB{it}(4,nTrk,iBB))+1;
%         for iTrk=1:nTrk
%             BBV(PDBB{it}(3,iTrk,iBB),PDBB{it}(4,iTrk,iBB))=v*vDir;
%         end
%     end
% end
%     
% fig=figure;
% set(fig,'Position',[200 300 700 300])
% h1=subplot(1,2,1);
% imagesc(BB0);
% colormap(h1,hot);
% caxis([0 3]);
% axis equal tight
% h2=subplot(1,2,2);
% imagesc(BBV);
% 
% colormap(h2,VzCmap1);
% caxis([-40 40]);
% axis equal tight
    
    