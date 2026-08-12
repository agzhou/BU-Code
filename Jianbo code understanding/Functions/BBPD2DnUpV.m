%% calculate velocty and ascending and descending vasculature
function [BB,BBV, BBVz]=BBPD2DnUpV(BBPD,P,thdTrk)
% input: track and paired BB info
if nargin <3
    thdTrk=3;
end
BB=zeros(P.Dim(1),P.Dim(2),2);
BBV=zeros(P.Dim(1),P.Dim(2),2);
BBVz=zeros(P.Dim(1),P.Dim(2),2);
nt=size(BBPD,1);
for it=1:nt
    inPDBB=size(BBPD{it},3);
    for iBB=1:inPDBB
        if ~isempty (BBPD{it})
            nTrk=max(BBPD{it}(5,:,iBB));
            if nTrk>=thdTrk % Trackable threshold
                for iTrk=1:nTrk
                    if BBPD{it}(1,iTrk,iBB)~=0
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
                        v=median((BBPD{it}(6,max(iTrk-1-startAdjust,1):min(iTrk+1+endAdjust,nTrk),iBB)));
                        vz=median((BBPD{it}(7,max(iTrk-1-startAdjust,1):min(iTrk+1+endAdjust,nTrk),iBB)));
%                         vz=((BBPD{it}(3,iTrk,iBB))-(BBPD{it}(1,iTrk,iBB)))*P.lPix/1e3*P.CCFR;
                        vDir=sign(mean(BBPD{it}(3,:,iBB)-BBPD{it}(1,:,iBB)));
                        % Path interpolation
                        if BBPD{it}(1,iTrk,iBB)==BBPD{it}(3,iTrk,iBB)
                            intZ=linspace(BBPD{it}(1,iTrk,iBB),BBPD{it}(3,iTrk,iBB),abs(BBPD{it}(4,iTrk,iBB)-BBPD{it}(2,iTrk,iBB))+1);
                            intX=linspace(BBPD{it}(2,iTrk,iBB),BBPD{it}(4,iTrk,iBB),abs(BBPD{it}(4,iTrk,iBB)-BBPD{it}(2,iTrk,iBB))+1);
                        else
                            intZ=linspace(BBPD{it}(1,iTrk,iBB),BBPD{it}(3,iTrk,iBB),abs(BBPD{it}(3,iTrk,iBB)-BBPD{it}(1,iTrk,iBB))+1);
                            intX = round(interp1([BBPD{it}(1,iTrk,iBB),BBPD{it}(3,iTrk,iBB)],[BBPD{it}(2,iTrk,iBB),BBPD{it}(4,iTrk,iBB)],intZ));
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
    end
end