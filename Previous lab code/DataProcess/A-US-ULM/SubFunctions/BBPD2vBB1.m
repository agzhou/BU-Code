%% calculate velocty and ascending and descending vasculature, no interpolation
function [BB,BBV, BBVz]=BBPD2vBB1(BBPD,PRSSinfo)
% input: track and paired BB info
BB=zeros(PRSSinfo.Dim(1),PRSSinfo.Dim(2),2);
BBV=zeros(PRSSinfo.Dim(1),PRSSinfo.Dim(2),2);
BBVz=zeros(PRSSinfo.Dim(1),PRSSinfo.Dim(2),2);
nt=size(BBPD,1);
for it=1:nt
    inPDBB=size(BBPD{it},3);
    for iBB=1:inPDBB
        if ~isempty (BBPD{it})
            Trks=(BBPD{it}(5,:,iBB));
            goodTrk=Trks(Trks>0);
            nTrk=numel(goodTrk);
            iBBPD=BBPD{it}(:,goodTrk,iBB);
            if nTrk>=PRSSinfo.thdTrk % Trackable threshold
                for iTrk=1:nTrk
                    if iBBPD(1,iTrk)~=0
%                         v=median((BBPD{it}(6,:,iBB)));
%                         vz=median((BBPD{it}(7,:,iBB)));
                        if iTrk+2>nTrk
                            startAdjust=iTrk+2-nTrk;
                        else
                            startAdjust=0;
                        end
                        if iTrk-2<1
                            endAdjust=1-(iTrk-2);
                        else
                            endAdjust=0;
                        end
%                         v=median(iBBPD(6,max(iTrk-1-startAdjust,1):min(iTrk+1+endAdjust,nTrk)));
%                         vz=median((iBBPD(7,max(iTrk-1-startAdjust,1):min(iTrk+1+endAdjust,nTrk))));
                        v=mean(iBBPD(6,max(iTrk-1-startAdjust,1):min(iTrk+1+endAdjust,nTrk)));
                        vz=mean((iBBPD(7,max(iTrk-1-startAdjust,1):min(iTrk+1+endAdjust,nTrk))));
%                         vz=((BBPD{it}(3,iTrk,iBB))-(BBPD{it}(1,iTrk,iBB)))*P.lPix/1e3*P.CCFR;
                        vDir=sign(mean(iBBPD(3,:)-iBBPD(1,:)));
                        % Path interpolation
                        if iBBPD(1,iTrk)==iBBPD(3,iTrk)
                            intZ=round(linspace(iBBPD(1,iTrk),iBBPD(3,iTrk),abs(iBBPD(4,iTrk)-iBBPD(2,iTrk))+1));
                            intX=round(linspace(iBBPD(2,iTrk),iBBPD(4,iTrk),abs(iBBPD(4,iTrk)-iBBPD(2,iTrk))+1));
                        else
                            intZ=(linspace(iBBPD(1,iTrk),iBBPD(3,iTrk),abs(iBBPD(3,iTrk)-iBBPD(1,iTrk))+1));
                            intX = round(interp1([iBBPD(1,iTrk),iBBPD(3,iTrk)],[iBBPD(2,iTrk),iBBPD(4,iTrk)],intZ));
                            intZ=round(intZ);
                        end
                        if vDir<=0 % upwards flow
%                             BBup(BBPD{it}(3,iTrk,iBB),BBPD{it}(4,iTrk,iBB))=BBup(BBPD{it}(3,iTrk,iBB),BBPD{it}(4,iTrk,iBB))+1;
                            for iC=1:numel(intZ)
                                BB(intZ(iC),intX(iC),1)=BB(intZ(iC),intX(iC),1)+1;
                                if BBV(intZ(iC),intX(iC),1)==0
                                    BBVz(intZ(iC),intX(iC),1)=abs(vz)*vDir;
                                    BBV(intZ(iC),intX(iC),1)=v*vDir;
                                else
                                    BBVz(intZ(iC),intX(iC),1)=(BBVz(intZ(iC),intX(iC),1)+abs(vz)*vDir)/2;
                                    BBV(intZ(iC),intX(iC),1)=(BBV(intZ(iC),intX(iC),1)+v*vDir)/2;
                                end
                            end
                        else % downwards flow
%                             BBdn(BBPD{it}(3,iTrk,iBB),BBPD{it}(4,iTrk,iBB))=BBdn(BBPD{it}(3,iTrk,iBB),BBPD{it}(4,iTrk,iBB))+1;
                            for iC=1:numel(intZ)
                                BB(intZ(iC),intX(iC),2)=BB(intZ(iC),intX(iC),2)+1;
                                if BBV(intZ(iC),intX(iC),2)==0
                                    BBVz(intZ(iC),intX(iC),2)=abs(vz)*vDir;
                                    BBV(intZ(iC),intX(iC),2)=v*vDir;
                                else
                                    BBVz(intZ(iC),intX(iC),2)=(BBVz(intZ(iC),intX(iC),2)+abs(vz)*vDir)/2;
                                    BBV(intZ(iC),intX(iC),2)=(BBV(intZ(iC),intX(iC),2)+v*vDir)/2;
                                end
                            end
                        end
                    end
                end
            end
            
        end
%         hold on,imagesc(BBV(:,:,2));caxis([0 20]); colormap(hot)
    end
end