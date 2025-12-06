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
function [BBPD]=BBPT2(CoorBB,PRSSinfo) % use the coordinates of BB
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
            else
%                 dCrit=(PRSSinfo.dCrit(1)/PRSSinfo.lPix);% maximum distance between BB(t+1) and BB(t) 
%                 dEstCrit=6; % maximum distance between BBest(t+1) and BB(t+1) 
                dCrit=min((mV)*2,6);% maximum distance between BB(t+1) and BB(t) 
                dEstCrit=min(mV*1,6); % maximum distance between BBest(t+1) and BB(t+1)
            end
            %% pairing frames ?current and next frames?
            if nMis==0
                CoorCBB=CoorBB{it+iTrkF-1}; % Coordinates of all identified bubble in the it+iTrkF-1 frame 
                CoorNBB=CoorBB{it+iTrkF}; % Coordinates of all identified bubble in the it+iTrkF frame
                dBBnext=sqrt(sum(((iBBCoorCBB-CoorNBB).^2),2)); % distance between the traking bubble (iBB) in the it+iTrkF-1 frame
                                                           % with all bubbles in the next frame
            else
                CoorCBB=CoorBB{it+iTrkF-1}; % Coordinates of all identified bubble in the it+iTrkF-1 frame 
                CoorCBB=[estCoorNPD;CoorCBB]; % put the estimated BB location in current (it+iTrkF-1) frame
                iBBCoorCBB=estCoorNPD;
                CoorNBB=CoorBB{it+iTrkF}; % Coordinates of all identified bubble in the it+iTrkF frame
                dBBnext=sqrt(sum(((estCoorNPD-CoorNBB).^2),2)); % distance between the estimated traking bubble (iBB) in the it+iTrkF-1 frame
                                                           % with all bubbles in the next frame
            end
            sPositionInd=find(dBBnext==0);
            dBBnext(sPositionInd)=[];
            CoorNBB(sPositionInd,:)=[];
            if ~isempty(dBBnext)
                indPD=find(dBBnext==min(dBBnext));
                iBBCoorNPD=CoorNBB(indPD(1),:); % the forward paired BB coordinate in the it+iTrkF frame
                dEST=sqrt(sum(((estCoorNPD-iBBCoorNPD).^2),2)); % identified position vs estimated position
                if (min(dBBnext)<dCrit || dEST<dEstCrit) && min(dBBnext)>0
                    %% backward pairing
                    dBBpre=sqrt(sum(((iBBCoorNPD-CoorCBB).^2),2)); % distance between the forward paried bubble (iBB) in the it+iTrkF frame with all bubbles in the previous frame
                    sPositionInd=find(dBBpre==0);
                    dBBpre(sPositionInd)=[];
                    CoorCBB(sPositionInd,:)=[];
                    indPD=find(dBBpre==min(dBBpre));
                    iBBCoorCPD=CoorCBB(indPD(1),:); % the backward paired BB coordinate in the it+iTrkF-1 frame of the forward paried BB
                    MutualMin=sum(iBBCoorCBB-iBBCoorCPD); % check if the BB in the it+iTrkF-1 frame is also the closest BB of the paried BB in the it+iTrkF frame
                    if (min(dBBnext)<dCrit || dEST<6) && MutualMin==0 % criteria of maximum distance for bubble pairing
%                         dircZ=iBBCoorNPD(1)-iBBCoorCBB(1);
%                         dircX=iBBCoorNPD(2)-iBBCoorCBB(2);
%                         if dircZ>=0 && dircX>=0
%                             dirc=1;
%                         elseif dircZ>=0 && dircX<0
%                             dirc=2;
%                         elseif dircZ<0 && dircX<0
%                             dirc=3;
%                         elseif dircZ<0 && dircX>=0
%                             dirc=4;
%                         end
                        v=min(dBBnext)*PRSSinfo.lPix/1e3*PRSSinfo.CCFR/(nMis+1);
                        vz=(iBBCoorNPD(1)-iBBCoorCBB(1))*PRSSinfo.lPix/1e3*PRSSinfo.CCFR/(nMis+1);
                        preMove=(iBBCoorNPD-iBBCoorCBB);
                        iBBPD(:,iTrkF,ipdBB)=[iBBCoorCBB(1),iBBCoorCBB(2),estCoorNPD(1),estCoorNPD(2),iTrkF,v,vz]'; % paired BB, [zCoor,xCoor,zCoorPair,xCoorPair,nTrackable,V,Vz,direction]
                        if nMis>0
                            v=abs(sqrt(sum(iBBCoorCBBpre-iBBCoorCBB).^2))*PRSSinfo.lPix/1e3*PRSSinfo.CCFR/(nMis+1);
                            vz=(iBBCoorCBBpre(1)-iBBCoorCBB(1))*PRSSinfo.lPix/1e3*PRSSinfo.CCFR/(nMis+1);
                            iBBPD(:,iTrkF-1,ipdBB)=[iBBCoorCBBpre(1),iBBCoorCBBpre(2),iBBCoorCBB(1),iBBCoorCBB(2),iTrkF-1,v,vz]'; % paired BB, [zCoor,xCoor,zCoorPair,xCoorPair,nTrackable,V,Vz,direction]
                        end
                        mV=mean(iBBPD(6,:,ipdBB),2);
                        % remove the paired BB from the it+iTrkF-1 frame
%                         temp=(CoorCBB-iBBCoorCBB);
%                         iCBBPDind=find(temp(:,1)==0);
%                         CoorBB{it+iTrkF-1}(iCBBPDind,:)=[];
                        estCoorNPD=iBBCoorNPD+(iBBCoorNPD-iBBCoorCBB);  % estimate the next possible position of iBB                        
                        iBBCoorCBB=iBBCoorNPD; % update the coordiante of iBB as the current frame (it+iTrkF)
                        
                        trackable=1;
                        nMis=0;
                    end
                else
                    nMis=nMis+1;
                    iBBCoorCBBpre=iBBCoorCBB; % save the previous coordiante of iBB (frame before the missing frame)
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
    
    