%% microbubble pair and track - based on closest criteria only
% function [BBPD]=BBPT(BB,PRSSinfo) % use 2D image of BB 
% % BB, bubble position obtained from Img2ULM
% % PRSSinfo.dCrit, pair distance criteria, in pixel. PRSSinfo.dCrit(1): maximum searching distance
% % between two frames; PRSSinfo.dCrit(2): maximum footprint pairing distance
% % BBPD, Paired and Tracked bubble 
% [nz,nx,nt]=size(BB); % nubmer of x,y,t points
% %% obtain all bubbles' coordinates
% for it=1:nt
%     [zBB,xBB]=find(BB(:,:,it)==1);
%     CoorBB{it}=[zBB,xBB];
% end
function [BBPD]=BBPT1(CoorBB,PRSSinfo) % use the coordinates of BB
nt=length(CoorBB);
BBPD=cell(nt-PRSSinfo.nTrack,1);
%% 
for it=1:nt-PRSSinfo.nTrack
    itCoorBBp=CoorBB{it};  % coordinates (z, x) of BB in frame it
    nBB=size(itCoorBBp,1); % number of total bubbles in the it frame
    ipdBB=1; % initialize number of paired BB in frame it
    %% Find out BB satisfying the criteria
    for iBB=1:nBB
        trackable=0; % check if iBB is trackable
        nMis=0;   % account for missed BB during BB identification
        for iTrkF=1:PRSSinfo.nTrack % try to pair and track the same bubble in the next PRSSinfo.nTrack frames
            if iTrkF==1
                iBBCoorCBB=itCoorBBp(iBB,:); % Coordinate of current frame bubble of iBB
                estCoorNPD=iBBCoorCBB;
                dCrit=12; % maximum distance between BB(t+1) and BB(t) 
                dEstCrit=12; % maximum distance between BBest(t+1) and BB(t+1) 
                mV=0;
            else
                dCrit=max(min((mV)*2,10),3);% maximum distance between BB(t+1) and BB(t) 
                dEstCrit=max(min(mV*2,6),2); % maximum distance between BBest(t+1) and BB(t+1)
            end
            %% pairing frames
            CoorCBB=CoorBB{it+iTrkF-1}; % Coordinates of all identified bubble in the it+iTrkF-1 frame 
            CoorNBB=CoorBB{it+iTrkF}; % Coordinates of all identified bubble in the it+iTrkF frame 
            % show the iBB
%             figure,plot(CoorCBB(:,2),CoorCBB(:,1),'r.'); hold on,plot(CoorNBB(:,2),CoorNBB(:,1),'b.')
%             hold on,plot(iBBCoorCBB(2),iBBCoorCBB(1),'r<'); hold on,plot(estCoorNPD(2),estCoorNPD(1),'bo');
            
            dBBEst=sqrt(sum(((estCoorNPD-CoorNBB).^2),2)); % distance between the estimated traking bubble (iBB) with all bubbles in the it+iTrkF frame 
            sPositionInd=find(dBBEst==0);
            dBBEst(sPositionInd)=[];
            CoorNBB(sPositionInd,:)=[];
            if ~isempty(dBBEst)
                dnbBBnext=min(dBBEst); % shortdistance between the closest BB in it+iTrkF frame
                indPD=find(dBBEst==dnbBBnext);
                iBBCoorNPD=CoorNBB(indPD(1),:); % the forward paired BB coordinate in the it+iTrkF frame
                dBB=sqrt(sum(((iBBCoorCBB-iBBCoorNPD).^2),2)); % distance of iBB(it+iTrkF-1) and iBBpd(it+iTrkF)
                if (dBB<dCrit || dnbBBnext<dEstCrit) && dnbBBnext>0
                    %% backward pairing
                    dBBpre=sqrt(sum(((iBBCoorNPD-CoorCBB).^2),2)); % distance between the forward paried bubble (iBB) in the it+iTrkF frame with all bubbles in the previous frame
                    sPositionInd=find(dBBpre==0);
                    dBBpre(sPositionInd)=[];
                    CoorCBB(sPositionInd,:)=[];
                    dnbBBpre=min(dBBpre); % shortest distance between BBs in previous frame
                    indPD=find(dBBpre==dnbBBpre);
                    iBBCoorCPD=CoorCBB(indPD(1),:); % the backward paired BB coordinate in the it+iTrkF-1 frame of the forward paried BB
                    MutualMin=sum(iBBCoorCBB-iBBCoorCPD); % check if the BB in the it+iTrkF-1 frame is also the closest BB of the paried BB in the it+iTrkF frame
                    if (dnbBBpre<dCrit || dBB<dEstCrit) && MutualMin==0  % criteria of maximum distance for bubble pairing
                        dircZ=iBBCoorNPD(1)-iBBCoorCBB(1);
                        dircX=iBBCoorNPD(2)-iBBCoorCBB(2);
                        if dircZ>=0 && dircX>=0
                            dirc=1;
                        elseif dircZ>=0 && dircX<0
                            dirc=2;
                        elseif dircZ<0 && dircX<0
                            dirc=3;
                        elseif dircZ<0 && dircX>=0
                            dirc=4;
                        end
                        v=min(dBBEst)*PRSSinfo.lPix/1e3*PRSSinfo.CCFR/(nMis+1);
                        vz=(iBBCoorNPD(1)-iBBCoorCBB(1))*PRSSinfo.lPix/1e3*PRSSinfo.CCFR/(nMis+1);
%                         iBBPD(:,iTrkF,ipdBB)=[iBBCoorCBB(1),iBBCoorCBB(2),estCoorNPD(1),estCoorNPD(2),iTrkF,v,vz,dirc]'; % paired BB, [zCoor,xCoor,zCoorPair,xCoorPair,nTrackable,V,Vz,direction]
                        iBBPD(:,iTrkF,ipdBB)=[iBBCoorCBB(1),iBBCoorCBB(2),iBBCoorNPD(1),iBBCoorNPD(2),iTrkF,v,vz,dirc]'; % paired BB, [zCoor,xCoor,zCoorPair,xCoorPair,nTrackable,V,Vz,direction]
                        mV=mean(iBBPD(6,max(iTrkF-2,1):end,ipdBB),2);
                        estCoorNPD=iBBCoorNPD+(iBBCoorNPD-iBBCoorCBB); % estimate the next possible position of the BB
                        % remove the paired BB from the it+iTrkF-1 frame
%                         temp=(CoorCBB-iBBCoorCBB);
%                         iCBBPDind=find(temp(:,1)==0);
%                         CoorBB{it+iTrkF-1}(iCBBPDind,:)=[];
                        iBBCoorCBB=iBBCoorNPD; % update the coordiante of iBB
                        
                        trackable=1;
                        nMis=0;
                    end
                else
                    nMis=nMis+1;
                    if nMis>1 || (iTrkF==1 && nMis==1) % account for missed BB during BB identification
                        break;
                    end
                end
            else
                break;
            end
        end
%         if iTrkF>2
%             iBBPD(:,1,ipdBB)=iBBPD(:,2,ipdBB);
%         end
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
% for it=1:nt-PRSSinfo.nTrack
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
    
    