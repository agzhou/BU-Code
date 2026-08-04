%% Generate SCC job file for RF2IQ data processing
function JOBS_RF2IQ(datapath,filename0, startFile,nFile)
jobfilename='RF2IQ.sh';
jobsavepath= '/projectnb/npbfus/s/Jianbo/CODE/BU-SCC/JOBS/';  % job.txt file save path
SubFunctionPath='/projectnb/npbfus/s/Jianbo/CODE/Functions'; % sub functions Path on SCC server
SCCFunctionPath='/projectnb/npbfus/s/Jianbo/CODE/BU-SCC/SCCFunctions'; % SCC submit function Path on SCC server
%% RF2IQ data process
% DeftPath='/projectnb/npboctiv/ns/Jianbo/USI/EXPERIMENT/0309BubblePhantom/';
% [filename0,datapath]=uigetfile(DeftPath);
%% Load RF data acquisition information
disp('Loading data and calculating the first IQ image...')
load([datapath, filename0]);
fileInfo=strsplit(filename0,'-');
RF0=reshape(RFRAW,[P.actZsamples,P.numAngles,P.numCCframes,P.nCh]);
%% Resample RF data
P.nRFref=1;
if P.nRFref==1
    RF=RF0;
else
    for iCC=1%:P.numCCframes
        for iAgl=1:P.numAngles
            iRF0=squeeze(RF0(:,iAgl,iCC,:));
            RF(:,iAgl,iCC,:)=imresize(iRF0,[P.actZsamples*P.nRFref,P.nCh]);
        end
    end
end
clear RF0
%% for in vivo data processing, to remove strong tissue/bone reflection signal
% for iAgl=1:P.numAngles
%     for iCh=1:P.nCh
%         RF(abs(squeeze(mean(RF(:,iAgl,:,iCh),3)))>mean(abs(squeeze(mean(RF(:,iAgl,:,iCh),3))))*2+std(abs(squeeze(mean(RF(:,iAgl,:,iCh),3))))*4,iAgl,:,iCh)=0;
%     end
% end
startDepthDAQ=P.startDepth;
% wavelength=P.vSound/(P.TWfrequency*1e3); % mm
wavelength=P.vSound/(18*1e3); % mm
zDelay=startDepthDAQ*wavelength/2;
%% calculate the reference matrix
if isprop(P,'nSmplPerWvlnth')==0
    P.nSmplPerWvlnth=4;
end
xCoor=[0:P.pitch:P.pitch*P.nCh];
zCoor=[zDelay:wavelength:(zDelay+15)];
nx=length(xCoor);
nz=length(zCoor);
IndCtriMatrix=zeros(nz,nx,P.nCh,P.numAngles);
ApodChn=zeros(nz,nx,P.nCh,P.numAngles);
P.startDepthRec=2*P.startDepth-15; % calibrated startDepth value
[IndCtriMatrix,ApodChn]=RefIndMatrix(P,xCoor,zCoor, 0.28);
%% calculate one compounded image for reference
for iCC=1 %P.numCCframes
    for iAgl=1:P.numAngles
        iRF=squeeze(RF(:,iAgl,iCC,:));
        iBF=zeros(nz,nx);
        for ix=1:nx
            for iz=1:nz
                [~, ctrChn]=find(ApodChn(iz,ix,:,iAgl)>0);
                for iCh=ctrChn(1):ctrChn(end)
                    iBF(iz,ix)=iBF(iz,ix)+iRF(IndCtriMatrix(iz,ix,iCh,iAgl),iCh);
                end
            end
        end
        IQ(:,:,iAgl,iCC)=iBF;%hilbert(iBF);
    end
    disp(['iCCFrame ', num2str(iCC),' is processed.'])
end
fig=figure;
set(fig,'Position',[200 500 400 300]);
imagesc(xCoor,zCoor,abs(hilbert(squeeze(sum(IQ(:,:,:,1),3)))));
colormap(gray)
xlabel('x [mm]');
ylabel('z [mm]');
axis equal tight
%% Beamform image parameter
prompt={['zImgStart, zDelay=', num2str(zDelay),' (mm)'],...
    ['zImg, zDelay=', num2str(zDelay),' (mm)'],...
    ['dzImg, wavelength=',num2str(wavelength),' (mm)'],...
    ['xImg, probeLength=',num2str(P.nCh*P.pitch),' (mm)'],...
    ['dxImg, pitch=',num2str(P.pitch),' (mm)'],'nRF refine','NA',...
    'start processing file','number of files to process',...
    ['number of split Segs for one SupFrame, nCC=',num2str(P.numCCframes),', nAgl=',num2str(P.numAngles)],...
    'Strong reflection signal thresholding? (Y: 1, N: 0)'};
name='Beamforming';
defaultvalue={'4','8', '0.05', num2str(P.nCh*P.pitch),num2str(P.pitch/2),'4','0.28',num2str(startFile),num2str(nFile),'1','0'};
numinput=inputdlg(prompt,name, 1, defaultvalue);
zDelay=str2num(numinput{1});
zImg=str2num(numinput{2});
dzImg=str2num(numinput{3});
xImg=str2num(numinput{4});
dxImg=str2num(numinput{5});
P.nRFref=str2num(numinput{6});
NA=str2num(numinput{7});
startFile=str2num(numinput{8});
nFile=str2num(numinput{9});
nSegPerFile=str2num(numinput{10});
RlcThrld=str2num(numinput{11});
close (fig)
%% beamform coordinates
clear IndCtriMatrix ApodChn xCoor zCoor
fileInfo{end-1}=num2str(startFile);
filename=[strjoin(fileInfo,'-')];
xCoor=[0:dxImg:xImg];
zCoor=[zDelay:dzImg:(zDelay+zImg)];
P.xCoor=xCoor;
P.dxImg=dxImg;
P.zCoor=zCoor;
P.dzImg=dzImg;
P.NA=NA;
%% calculate the Delay&Sum reference matrix
clear IndCtriMatrix ApodChn
[IndCtriMatrix,ApodChn]=RefIndMatrix(P,xCoor,zCoor, NA);
save([datapath,'/', strjoin(fileInfo(1:6),'-'),'-BFMatrix.mat'],'IndCtriMatrix','ApodChn','P');
% for ifile=startFile:startFile-1+nFile
%     fileinfo{end-1}=num2str(ifile);
%     ifilename=[strjoin(fileinfo,'-')];
%     %% generate subfolder for processed IQ data
%     subfolder=strjoin(fileinfo(1:7),'-');
%     savepath=[datapath,'/',subfolder,'/'];
%     if exist(savepath)
%     else
%         mkdir(datapath,subfolder);
%     end
%     save([savepath,'/', ifilename(1:end-4),'-BFMatrix.mat'],'IndCtriMatrix','ApodChn','NA','xCoor','zCoor');
% end
%% SCC file info
PRCSinfo(1)=nSegPerFile;    % number of splited chunks for each super frame, 
nCCpSeg=P.numCCframes/nSegPerFile;
PRCSinfo(2)=nCCpSeg; % each chunk/segment contains nCCpSeg CCframes
PRCSinfo(3)=RlcThrld; % threshold out strong tissue/bone reflection
%% SCC node/core request 
pathInfo=strsplit(datapath,'/');
jobname=[pathInfo{end-1},'-',filename];
prompt2={'# jobs','# cores per job','Memory per core','wall time','job name','email notify? (Y/N)'};
inputSCCinfo=inputdlg(prompt2,'SCC request parameter', 1,{num2str(nSegPerFile*nFile),'2','8','6',jobname,'N'});
Njobs=(inputSCCinfo{1});       % number of jobs to be submitted for each file
Ncores=(inputSCCinfo{2});      % number of cores requested for each job
MemperCore=(inputSCCinfo{3});  % memory per core
WallTime=(inputSCCinfo{4});    % specify wall time
JobName=['RF2IQ-',inputSCCinfo{5}];              % Job name
EmailNote=inputSCCinfo{6};             % Email notify ?
if EmailNote=='Y' || EmailNote=='y'
    prompt4={'email address'};
    inputEmail=inputdlg(prompt4,'User email', 1,{'jianbo@bu.edu'});
    Emadd=(inputEmail{1});       % user's email address
    Emailinfo=['#$ -m ea','\n',...
        ['#$ -M ' Emadd],'\n'];
else
    Emailinfo=['\n'];
end

%%%%
fileInfo=strsplit(filename,'-');
%% save RF2IQ process information
fid=fopen([jobsavepath, jobfilename],'wt');
fileNameBase=[strjoin(fileInfo(1:6),'-'),'-'];
fileNameTail=['-',fileInfo{end}];
JobName=['RF2IQ-',pathInfo{end-1},'-',filename];

job_cmd=['#! /bin/bash -l','\n'...
    ['#$ -pe omp ', Ncores],'\n',...
    ['#$ -l mem_per_core=',MemperCore,'G'], '\n',...
    ['#$ -l h_rt=',WallTime,':00:00'], '\n',...
    ['#$ -N ',JobName],'\n',...
    ['#$ -t 1-',Njobs],'\n',...
    Emailinfo,...
    'module load matlab/2017a','\n',...
    ['matlab -nodisplay -r ',...
    '"diary on; '...
    'addpath(''',SCCFunctionPath,'''); '...
    'addpath(''',SubFunctionPath,'''); ',...
    'datapath=''', datapath,'''; ',...
    'iJob=$SGE_TASK_ID; ',...
    'nSegPerFile=',num2str(nSegPerFile),';',...
    'ifile=',num2str(startFile),'+floor((iJob-1)/nSegPerFile);',...
    'iSeg=iJob-(ifile-',num2str(startFile),')*nSegPerFile;',...
    'fileNameBase=''',fileNameBase,''';',...
    'fileNameTail=''',fileNameTail,''';',...
    'filename=[fileNameBase, num2str(ifile),fileNameTail]; ',...
    'PRCSinfo=[', num2str(PRCSinfo), ']; ',...
    'SCC_RF2IQ(datapath, filename,PRCSinfo, iSeg); ',...
    'diary off; exit"'],'\n\n'];
fprintf(fid,job_cmd);
fclose(fid);
%%

%% save RF2IQ process information in data folder

% for ifile=startFile:startFile-1+nFile
%     fileInfo{end-1}=num2str(ifile);
%     ifilename=[strjoin(fileInfo,'-')];
%     %% generate subfolder for processed IQ data
%     subfolder=strjoin(fileInfo(1:7),'-');
%     savepath=[datapath,'/',subfolder,'/'];
%     
%     fid=fopen([savepath, jobfilename],'wt');
%     filename00=['',ifilename,''];
%     job_cmd=['#! /bin/bash -l','\n'...
%         ['#$ -pe omp ', Ncores],'\n',...
%         ['#$ -l mem_per_core=',MemperCore,'G'], '\n',...
%         ['#$ -l h_rt=',WallTime,':00:00'], '\n',...
%         ['#$ -N ',JobName],'\n',...
%         ['#$ -t 1-',Njobs],'\n',...
%         Emailinfo,...
%         'module load matlab/2017a','\n',...
%         ['matlab -nodisplay -r ',...
%         '"diary on; '...
%         'addpath(''',SCCFunctionPath,'''); '...
%         'addpath(''',SubFunctionPath,'''); ',...
%         'datapath=''', datapath,'''; ',...
%         'filename=''',filename00,'''; ',...
%         'PRCSinfo=[', num2str(PRCSinfo), ']; ',...
%         'iSeg=$SGE_TASK_ID; ',...
%         'SCC_RF2IQ(datapath, filename,PRCSinfo, iSeg); ',...
%         'diary off; exit"'],'\n'];
%     fprintf(fid,job_cmd);
%     fclose(fid);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Job RF2IQ.txt saved,', datestr(now,'DD:HH:MM')]);
disp([datapath,filename]);