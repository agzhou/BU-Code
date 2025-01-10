%% SCC function for Coherence compounding Beamform RF data

function SCC_RF2IQ(datapath, filename, PRCSinfo, iSeg)
nSeg=PRCSinfo(1);    % number of splited chunks for each super frame, 
nCCpSeg=PRCSinfo(2); % each chunk/segment contains nCCpSeg CCframes
RlcThrld=PRCSinfo(3);
% note: nCC per SupFrame=nSeg*nCCpSeg
fileInfo=strsplit(filename,'-');

fileRef=matfile([datapath,'/',filename]);
P=fileRef.P;
%% generate subfolder for processed IQ data
subfolder=strjoin(fileInfo(1:7),'-');
% savepath=[datapath,'/',subfolder,'/'];
savepath=[datapath,'/'];
if exist(savepath)
else
    mkdir(datapath,subfolder);
end
%% load the beamforming matrix file, which includes IndCtriMatrix, ApodChn, NA, xCoor, and zCoor
xElemCoor=(0.5:P.nCh)*P.pitch;   
if exist([datapath,'/', strjoin(fileInfo(1:6),'-'),'-BFMatrix.mat'])
    load([datapath,'/', strjoin(fileInfo(1:6),'-'),'-BFMatrix.mat']);
    xCoor=P.xCoor;
    zCoor=P.zCoor;
    NA=P.NA;
    nRFref=P.nRFref;
    nz=length(P.zCoor);
    nx=length(P.xCoor);
    %% RF2IQ for ith seg
    disp('Loading Data');
    RFRAW=fileRef.RFRAW((P.actZsamples*P.numAngles)*(iSeg-1)*nCCpSeg+1:(P.actZsamples*P.numAngles)*(iSeg*nCCpSeg),:);
    RF=reshape(RFRAW,[P.actZsamples, P.numAngles, nCCpSeg,P.nCh]);
    %% FOR in vivo data processing only, to remove strong tissue/bone reflection
    if RlcThrld==1
        for iAgl=1:P.numAngles
            for iCh=1:P.nCh
                RF(abs(squeeze(mean(RF(:,iAgl,:,iCh),3)))>mean(abs(squeeze(mean(RF(:,iAgl,:,iCh),3))))*2+std(abs(squeeze(mean(RF(:,iAgl,:,iCh),3))))*3,iAgl,:,iCh)=0;
            end
        end
        disp('Data loaded');
    end
    %% RF2IQ
    tic
    IQ=double(zeros(nz,nx,P.numAngles,P.numCCframes));
    for iAgl=1:P.numAngles
        %% Resample RF data
        iRF=zeros(P.actZsamples*P.nRFref,P.numCCframes,P.nCh);
        if P.nRFref==1
            iRF=squeeze(RF(:,iAgl,:,:));
        else
            for iCh=1:P.nCh
                iRF0=squeeze(RF(:,iAgl,:,iCh));
                iRF(:,:,iCh)=imresize(iRF0,[P.actZsamples*P.nRFref,P.numCCframes]);
            end
        end
        %% Beamforming
        for ix=1:nx
            for iz=1:nz
                [~, ctrChn]=find(ApodChn(iz,ix,:,iAgl)>0);
                nCtrChn=numel(ctrChn);
                for iCh=ctrChn(1):ctrChn(end)
                    WT=exp(-(xCoor(ix)-xElemCoor(iCh))^2/(2*(zCoor(iz)*NA)^2));
                    IQ(iz,ix,iAgl,:)=squeeze(IQ(iz,ix,iAgl,:))'+double(squeeze(iRF(IndCtriMatrix(iz,ix,iCh,iAgl),:,iCh)))*WT;
                end
                IQ(iz,ix,iAgl,:)=IQ(iz,ix,iAgl,:)/nCtrChn;
            end
        end
        disp(['iAngle ', num2str(iAgl),' is processed.'])
    end
%     for iCC=1:nCCpSeg
%         for iAgl=1:P.numAngles
%             iRF=squeeze(RF(:,iAgl,iCC,:));
%             iBF=zeros(nz,nx);
%             for ix=1:nx
%                 for iz=1:nz
%                     [~, ctrChn]=find(ApodChn(iz,ix,:,iAgl)>0);
%                     for iCh=ctrChn(1):ctrChn(end)
%                         iBF(iz,ix)=iBF(iz,ix)+iRF(IndCtriMatrix(iz,ix,iCh,iAgl),iCh);
%                     end
%                 end
%             end
%             BFP(:,:,iAgl,iCC)=iBF;%hilbert(iBF);
%         end
%         
%         disp(['iCCFrame ', num2str(iCC),' is processed.'])
%     end
    % end
    toc
    %% data saving
    P.Prcs_nSeg=nSeg;
    P.Prcs_iSeg=iSeg;
    P.Prcs_nCCpSeg=nCCpSeg;
    NameIQ=[filename(1:end-7),'-',num2str(nSeg),'-',num2str(iSeg),'-',num2str(nCCpSeg),'-IQ'];   % save name for coherence compounded IQ data
    disp('Data saving...');
    save([savepath,NameIQ,'.mat'],'-v7.3','IQ','P');
    disp('Data saved!');
else
    disp('Error: No Delay&Sum reference matrix existed!')
end
