%% calculate velocty and ascending and descending vasculature
function [BB,BBV,BBVz]=BBPD2vBB2(BBPD,PRSSinfo)
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
                mV =movmean(squeeze(iBBPD(6,1:nTrk)),3);
                mVz=movmean(squeeze(iBBPD(7,1:nTrk)),3);
                vDir=sign(mean(iBBPD(3,:)-iBBPD(1,:)));
                vDir(vDir==0)=1;
                vDir(vDir==-1)=0;
                vDir=vDir+1; % 1: upwards flow; 2: downwards flow
                for iTrk=1:nTrk
                    if iBBPD(1,iTrk)~=0
                        v=mV(iTrk);
                        vz=mVz(iTrk);
                        % Path interpolation
                        lZ=ceil(abs(iBBPD(1,iTrk)-iBBPD(3,iTrk)));
                        lX=ceil(abs(iBBPD(2,iTrk)-iBBPD(4,iTrk)));
                        nPixInt=max(lZ,lX);
                        intZ=round(linspace(iBBPD(1,iTrk),iBBPD(3,iTrk),nPixInt));
                        intX=round(linspace(iBBPD(2,iTrk),iBBPD(4,iTrk),nPixInt));
                        intV=(linspace(mV(iTrk),mV(min(iTrk+1,nTrk)),nPixInt));      
                        intVz=(linspace(mVz(iTrk),mVz(min(iTrk+1,nTrk)),nPixInt));   
                        for iC=1:numel(intZ)
                            BB(intZ(iC),intX(iC),vDir)=BB(intZ(iC),intX(iC),vDir)+1;
                            BBVz(intZ(iC),intX(iC),vDir)=(BBVz(intZ(iC),intX(iC),vDir)+abs(intVz(iC))*sign(vDir-1.5));
                            BBV(intZ(iC),intX(iC),vDir)=(BBV(intZ(iC),intX(iC),vDir)+intV(iC)*sign(vDir-1.5));
                        end
%                         for iC=1:numel(intZ)
%                             BB(intZ(iC),intX(iC),vDir)=BB(intZ(iC),intX(iC),vDir)+1;
%                             BBVz(intZ(iC),intX(iC),vDir)=(BBVz(intZ(iC),intX(iC),vDir)+abs(vz)*sign(vDir-1.5));
%                             BBV(intZ(iC),intX(iC),vDir)=(BBV(intZ(iC),intX(iC),vDir)+v*sign(vDir-1.5));
%                         end
                    end
                end
            end
            
        end
%         hold on,imagesc(BBV(:,:,2));caxis([0 20]); colormap(hot)
    end
end