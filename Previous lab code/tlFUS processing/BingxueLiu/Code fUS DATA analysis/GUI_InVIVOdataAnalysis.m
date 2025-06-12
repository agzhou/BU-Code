function varargout = GUI_InVIVOdataAnalysis(varargin)
% GUI_INVIVODATAANALYSIS MATLAB code for GUI_InVIVOdataAnalysis.fig
%      GUI_INVIVODATAANALYSIS, by itself, creates a new GUI_INVIVODATAANALYSIS or raises the existing
%      singleton*.
%
%      H = GUI_INVIVODATAANALYSIS returns the handle to a new GUI_INVIVODATAANALYSIS or the handle to
%      the existing singleton*.
%
%      GUI_INVIVODATAANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_INVIVODATAANALYSIS.M with the given input arguments.
%
%      GUI_INVIVODATAANALYSIS('Property','Value',...) creates a new GUI_INVIVODATAANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_InVIVOdataAnalysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_InVIVOdataAnalysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_InVIVOdataAnalysis

% Last Modified by GUIDE v2.5 20-Aug-2020 11:24:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_InVIVOdataAnalysis_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_InVIVOdataAnalysis_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before GUI_InVIVOdataAnalysis is made visible.
function GUI_InVIVOdataAnalysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_InVIVOdataAnalysis (see VARARGIN)

% Choose default command line output for GUI_InVIVOdataAnalysis
handles.output = hObject;
handles.output = hObject;
handles.DefPath='G:\PROJ-R-Stroke';
addpath('D:\BingxueLiu\Code fUS DATA analysis')
addpath('D:\BingxueLiu\Code fUS DATA analysis\SubFunctions')
handles.DnUpSlt=1;
handles.DataSlt=1;
handles.nROI=1;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_InVIVOdataAnalysis wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_InVIVOdataAnalysis_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in BTN_LoadData.
function BTN_LoadData_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_LoadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;
addpath('D:\CODE\DataProcess\BingxueLiu\Code fUS DATA analysis')
addpath('D:\CODE\DataProcess\BingxueLiu\Code fUS DATA analysis\SubFunctions')
handles.vUSorPDI = questdlg('Data Type?', ...
                         'Select', ...
                         'PDI', 'vUS', 'Cancel','Cancel');
handles.V0=[];
handles.Vz0=[];
handles.Vcz0=[];
handles.Vx0=[];
handles.VxR0=[];
handles.Msk0=[];
handles.PDI0=[];
handles.PDIHP0=[];
handles.PDIHHP0=[];
handles.PDISVD0=[];
handles.eqNoise0=[];
handles.eqNoise=[];
handles.BB=[];
switch handles.vUSorPDI
    case 'PDI'
        %% data loading
        [FileName,FilePath]=uigetfile(handles.DefPath);
        handles.DefPath=FilePath;
        guidata(hObject, handles);
        
        fileInfo=strsplit(FileName(1:end-4),'-');
        startCP0=fileInfo{7}(3:end);
        if isempty(startCP0)
            startCP0=0;
        end
        myFile=matfile([FilePath,FileName]);
        P=myFile.P;
        % data processing parameters
        prompt={'Start Repeat', 'Number of Repeats','nRfn','vUS Time Interval [s]'};
        name='File info';
        defaultvalue={num2str(fileInfo{8}),'270','1',num2str(P.tIntPDI*1e-6)};
        numinput=inputdlg(prompt,name, 1, defaultvalue);
        handles.startRpt=str2num(numinput{1});
        handles.nRpt=str2num(numinput{2});          % number of repeat for each coronal plane
        nRfn=str2num(numinput{3});          % image refine scale
        handles.tIntV=str2num(numinput{4});
        guidata(hObject, handles);
        % load all velocity data
        indSkipped=1;
        for iRpt=handles.startRpt:handles.startRpt+handles.nRpt-1
            iFileInfo=fileInfo;
            iFileInfo{8}=num2str(iRpt);
            iFileName=[strjoin(iFileInfo,'-'),'.mat'];
            %     load([FilePath,iFileName]);
            if exist([FilePath,iFileName],'file')
                myFile=matfile([FilePath,iFileName]);
                PDI=myFile.PDI;
                PDIHP=myFile.PDIHP; %PDIHP
                PDIHHP=myFile.PDIHHP;%PDIHHP
                eqNoise=myFile.eqNoise;
                disp([iFileName,' was loaded!'])
            else
                disp([iFileName, ' skipped!'])
                handles.SkipFile(indSkipped)=iRpt;
                indSkipped=indSkipped+1;
            end
            handles.PDI0(:,:,:,iRpt)=PDI;
            handles.PDIHP0(:,:,:,iRpt)=PDIHP;
            handles.PDIHHP0(:,:,:,iRpt)=PDIHHP;
            handles.eqNoise(:,:,iRpt)=eqNoise;
        end
    case 'vUS'
        [FileName,FilePath]=uigetfile(handles.DefPath);
        handles.DefPath=FilePath;
        guidata(hObject, handles);
        
        fileInfo=strsplit(FileName(1:end-4),'-');
        startCP0=fileInfo{7}(3:end);
        if isempty(startCP0)
            startCP0=0;
        end
        myFile=matfile([FilePath,FileName]);
        P=myFile.P;
        % data processing parameters
        prompt={'Start Repeat', 'Number of Repeats','nRfn','vUS Time Interval [s]'};
        name='File info';
        defaultvalue={num2str(fileInfo{8}),'270','1',num2str(P.tIntPDI*1e-6)};
        numinput=inputdlg(prompt,name, 1, defaultvalue);
        handles.startRpt=str2num(numinput{1});
        handles.nRpt=str2num(numinput{2});          % number of repeat for each coronal plane
        nRfn=str2num(numinput{3});          % image refine scale
        handles.tIntV=str2num(numinput{4});
        guidata(hObject, handles);
        % load all velocity data
        indSkipped=1;
        for iRpt=handles.startRpt:handles.startRpt+handles.nRpt-1
            iFileInfo=fileInfo;
            iFileInfo{8}=num2str(iRpt);
            iFileName=[strjoin(iFileInfo,'-'),'.mat'];
            %     load([FilePath,iFileName]);
            if exist([FilePath,iFileName],'file')
                myFile=matfile([FilePath,iFileName]);
                R=myFile.R;
                Mf=myFile.Mf;
                Vx=myFile.Vx;
                Vz=myFile.Vz;
                Vcz=myFile.Vcz;
                PDI=myFile.PDI;
%                 PDIHP=myFile.PDIHP;
                PDIHHP=myFile.PDIHHP;
                eqNoise=myFile.eqNoise;
                disp([iFileName,' was loaded!'])
            else
                disp([iFileName, ' skipped!'])
                handles.SkipFile(indSkipped)=iRpt;
                indSkipped=indSkipped+1;
            end
            %% Vx weighted with R
%             Rmsk=0.01*ones(size(R));
%             Rmsk(R>0.1)=R(R>0.1);
%             Rmsk(Rmsk>0.3)=1;
%             %         VxR=Vx(:,:,1:2).*Rmsk(:,:,1:2);
%             
%             MfR=Mf(:,:,1:2);
%             MfRmsk=0.01*ones(size(MfR));
%             MfRmsk(MfR>0.1)=MfR(MfR>0.1);
%             MfRmsk(MfRmsk>0.4)=1;
%             VxR=Vx.*MfRmsk(:,:,1:2).*Rmsk(:,:,1:2);
%             VzR=Vz.*Rmsk(:,:,1:2);
            
            VzR=Vz;
            VxR=Vx;
            %% all velocity
            handles.V0(:,:,:,iRpt)=sign(VzR(:,:,1:2)).*sqrt(VxR(:,:,1:2).^2+VzR(:,:,1:2).^2);
            handles.Vz0(:,:,:,iRpt)=VzR(:,:,1:2);
            handles.Vcz0(:,:,:,iRpt)=Vcz(:,:,1:2);
            handles.Vx0(:,:,:,iRpt)=Vx(:,:,1:2);
            handles.VxR0(:,:,:,iRpt)=VxR(:,:,1:2);
%             handles.Msk(:,:,:,iRpt)=MfRmsk;
            handles.PDI0(:,:,:,iRpt)=PDI;
%             handles.PDIHP0(:,:,:,iRpt)=PDIHP;
            handles.PDIHHP0(:,:,:,iRpt)=PDIHHP;
            handles.eqNoise(:,:,iRpt)=eqNoise;
        end
end
eqNoise=mean(handles.eqNoise,3);
handles.PDI0=handles.PDI0./eqNoise.^1.8;
handles.PDIHP0=handles.PDIHP0./eqNoise.^1.8;
handles.PDIHHP0=handles.PDIHHP0./eqNoise.^1.8;

ButtonName = questdlg('Stimulation Pattern?', ...
                         'Select', ...
                         'NEW', 'Load Exist', 'UseExist','Cancel');
switch ButtonName
    case 'NEW'
        prompt={'nBaseline','nStim','nRecover', 'nStimTrial', 'nSample'};
        name='Stimulation info';
        defaultvalue={'30','2','30','10',num2str(handles.nRpt)};
        numinput=inputdlg(prompt,name, 1, defaultvalue);
        handles.nBase=str2num(numinput{1});
        handles.nStim=str2num(numinput{2});
        handles.nRecover=str2num(numinput{3});
        handles.nStimTrial=str2num(numinput{4});
        handles.nSample=str2num(numinput{5});
        for iStim=1:handles.nStimTrial
            prompt{iStim}=['Stim on ', num2str(iStim)];
            defaultvalue{iStim}=[num2str(handles.nBase+(iStim-1)*(handles.nStim+handles.nRecover)+1),':',...
                num2str(handles.nBase+(iStim-1)*(handles.nStim+handles.nRecover)+handles.nStim)];
        end
        numinput=inputdlg(prompt,name, 1, defaultvalue);
        handles.Stim0=zeros(1,handles.nSample);
        for iStim=1:handles.nStimTrial
            handles.Stim0(str2num(numinput{iStim}))=1;
        end
        Stim=handles.Stim0;
        save([handles.DefPath,'StimPattern.mat'],'Stim');
    case 'Load Exist'
        disp('Select the Stimulation Pattern file')
        [FileName,FilePath]=uigetfile(handles.DefPath);  % read data of a small part of the brain cortex (IQR matrix)
        load([FilePath,FileName]);
        handles.Stim0=Stim;
end

handles.P=P;
guidata(hObject, handles);
%% using literature data for code test
% [FileName,FilePath]=uigetfile(handles.DefPath);
% handles.DefPath=FilePath;
% guidata(hObject, handles);
% load([FilePath,FileName]);
% handles.PDI(:,:,1,:)=PDI;
% handles.PDI(:,:,2,:)=PDI;
% handles.PDI(:,:,3,:)=PDI;
% handles.tIntV=1.5;
% [nz,nx,handles.nRpt]=size(PDI);
% handles.startRpt=1;
% guidata(hObject, handles);


% --- Executes on button press in BTN_MotionCorrect.
function BTN_MotionCorrect_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_MotionCorrect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;
addpath('D:\CODE\DataProcess\BingxueLiu\Code fUS DATA analysis')
addpath('D:\CODE\DataProcess\BingxueLiu\Code fUS DATA analysis\SubFunctions')
%% motion artifacts correction
handles.nRest = 5;
nRest = handles.nRest;
handles.PDI0 = handles.PDIHP0; % use SVD + highpass(30Hz) filtered data
nRpt = handles.nRpt; % save the nRpt information; modified by Bingxue
disp('Motion artifacts correction...')
bulkPDI=squeeze(mean(mean(handles.PDI0(:,:,3,:),1),2));
bulkThd=median(bulkPDI)+3*std(medfilt1(bulkPDI,10,'truncate'));
bulkIndex=(bulkPDI>bulkThd);
figure;
plot(squeeze(mean(mean(handles.PDI0(:,:,3,:),1),2)));
hold on, plot(repmat(bulkThd,handles.nRpt))
hold on, plot(handles.Stim0*(bulkThd-min((bulkPDI)))+min((bulkPDI)),'r')
ylim([min((bulkPDI)), bulkThd*1.4])
% for iStim=1:handles.nStimTrial
%     text(handles.nBase+(iStim-1)*(handles.nStim+handles.nRecover)+1, double(bulkThd),num2str(iStim),'FontSize',12);
% end
title('mean-original PDI')
MCORRslt = questdlg('Do Motion Correction?', ...
                         'Select', ...
                         'YES', 'NO', 'Cancel');
switch MCORRslt
    case 'YES'
        RMVslt = questdlg('Remove strong motion trial?', ...
            'Select', ...
            'NO', 'YES','Cancel');
        switch RMVslt
            case 'YES'
                prompt={'Index of trial removal', };
                name='Remove strong motion trial (0: no remove)';
                defaultvalue={'0'};
                numinput=inputdlg(prompt,name, 1, defaultvalue);
                RmvIndex=str2num(numinput{1});
                if RmvIndex(1)==0
                    PDI=handles.PDI0;
                    PDIHHP=handles.PDIHHP0;
                    V=handles.V0;
                    Vx=handles.Vx0;
                    Vz=handles.Vz0;
                    Vcz=handles.Vcz0;
                    handles.Stim=handles.Stim0;
                    handles.nTrial=handles.nStimTrial;
                    PDIrm = handles.PDI0;
                else
                    handles.nRpt= size(handles.PDI0,4); %modified by Bingxue
                    handles.nTrial=handles.nStimTrial-numel(RmvIndex);
                    Samples=ones(handles.nRpt,1);
                    stimSamples = ones(handles.nRpt,1);
                    for iTrial=1:numel(RmvIndex)
                        Samples((handles.nBase-nRest+(RmvIndex(iTrial)-1)*(handles.nStim+handles.nRecover)+1):...
                            (handles.nBase-nRest+(RmvIndex(iTrial))*(handles.nStim+handles.nRecover)))=0;  % modified by Bingxue, keep 5frames before stim trials
                        stimSamples((handles.nBase-nRest+(RmvIndex(iTrial)-1)*(handles.nStim+handles.nRecover)+1):...
                            (handles.nBase-nRest+(RmvIndex(iTrial))*(handles.nStim+handles.nRecover)))=0;
                    end
                    Samples=(Samples>0);
                    stimSamples=(stimSamples>0);
                    PDIrm=handles.PDI0(:,:,:,Samples); % PDIrm: remove motion trial
                    PDIHHP=handles.PDIHHP0(:,:,:,Samples);
%                     V=handles.V0(:,:,:,Samples);
%                     Vx=handles.Vx0(:,:,:,Samples);
%                     Vz=handles.Vz0(:,:,:,Samples);
%                     Vcz=handles.Vcz0(:,:,:,Samples); %comment for PDI processing
                    handles.Stim=handles.Stim0(stimSamples);
                    handles.Samples = Samples;
                end
            case 'NO'
                PDIrm=handles.PDI0;
                PDIHHP=handles.PDIHHP0;
                V=handles.V0;
                Vx=handles.Vx0;
                Vz=handles.Vz0;
                Vcz=handles.Vcz0;
                handles.Stim=handles.Stim0;
                handles.nTrial=handles.nStimTrial;
        end
        
          
     Accelslt = questdlg('Use accelerometer data?', ...
            'Select', ...
            'YES', 'NO','Cancel'); % modified by Bingxue
   switch Accelslt
       case 'YES'
     handles.Acc = [];
     handles.Acc = getAcceldata(handles.Stim0);
   end
        PDI = PDIrm;        
        handles.nRpt= size(PDI,4); % size(Vx,4); modified by Bingxue
        % PDI intensity correction
        [nx, ny, nt] = size(squeeze(PDI(:,:,3,:))); % nt is data points the after removing motion trials; 
        handles.PDI=PDI;
        bulkPDI=squeeze(mean(mean(handles.PDI(end-5:end,:,3,:),1),2)); % bottom for clutter area
        bulkThd=median(bulkPDI)+5*1.48*mad(bulkPDI,1); % modified by Bingxue Liu, use truncate to smooth edges
        bulkIndex=(bulkPDI<bulkThd);% index for kepted data
        Fig=figure;set(Fig, 'Position',[400 300 1000 300]);
        plot(bulkPDI);hold on; plot(bulkThd*ones(size(bulkPDI))); hold on; area(1:handles.nRpt, handles.Stim*max(bulkPDI).*ones(1,handles.nRpt),'FaceColor','r','EdgeColor','none'); alpha(0.1);
        %%% replace %%%
%         %%%% trial
        vrmIndex = (bulkPDI>bulkThd); % index for removed data
        keptIndex = reshape(bulkIndex(handles.nBase+1:end), [handles.nStim+handles.nRecover, handles.nTrial]);
        rmIndex = reshape(vrmIndex(handles.nBase+1:end), [handles.nStim+handles.nRecover, handles.nTrial]);
         figure; imagesc(keptIndex); axis image;
        PDI(:,:,3,vrmIndex) = zeros([nx, ny, sum(vrmIndex(:))]); 
        trialPDI = reshape(PDI(:,:,3,handles.nBase+1:end), [nx, ny, handles.nStim+handles.nRecover, handles.nTrial]);
        vtrialPDI = squeeze(PDI(:,:,3,handles.nBase+1:end));
        rplcPDI = sum(trialPDI, 4)./reshape(repmat(sum(keptIndex,2)',[nx*ny,1]),[nx, ny,handles.nStim+handles.nRecover]);
        rplcIndex = rmIndex.*repmat([1:(handles.nStim+handles.nRecover)]',[1,handles.nTrial]);
        for it = 1: handles.nStim+handles.nRecover
            vtrialPDI(:,:,find(rplcIndex == it)) = repmat(rplcPDI(:,:,it),[1,1,size(find(rplcIndex == it),1)]) ;
        end
        
        baseIndex = bulkIndex(1:handles.nBase);
        basePDI = squeeze(PDI(:,:,3,1:handles.nBase));
        rplcBasePDI = sum(squeeze(PDI(:,:,3,baseIndex)),3)./(sum(baseIndex)*ones(nx, ny));
        basePDI(:,:,find(baseIndex==0)) = repmat(rplcBasePDI, [1,1, size(find(baseIndex==0),1)]);
        PDI(:,:,3,:) = cat(3, basePDI, vtrialPDI);
%         %%% interp
%         vt0 = 1:nt;
%         vt = vt0(bulkIndex);
%         [x0,y0,z0] = meshgrid(1:ny,1:nx,vt);
%         [x1,y1,z1] = meshgrid(1:ny,1:nx,1:nt);
%         PDIclr = PDI(:,:,:,bulkIndex); % clutter rejection
%         interpPDI = interp3(x0,y0,z0,squeeze(PDIclr(:,:,3,:)),x1,y1,z1);
% %         rplc=repmat(median(PDI,4)+1.2*std(medfilt1(PDI,15,[],4,'truncate'),1,4),[1,1,1,handles.nRpt]); % modified by Bingxue Liu, replace movemean by medfilt
% %         PDI(:,:,:,bulkIndex)=rplc(:,:,:,bulkIndex);
% %         PDI=medfilt1(PDI,10,[],4,'truncate');% modified by Bingxue Liu, replace movmean by medfilt
% %         handles.PDI(:,:,:,bulkIndex)=PDI(:,:,:,bulkIndex);
%         PDI(:,:,3,4:end) = interpPDI(:,:,4:end);% first 3 are 0

        %%% generate msk based on raw data
%         m1 = squeeze(squeeze(PDIrm(150,100,3,:))); % 150 120 (0528 data)
%         m2 = squeeze(squeeze(PDIrm(94, 12, 3, :))); % 100 155 (0528 data)
%         m3 = squeeze(squeeze(PDIrm(44, 154, 3, :)));
%         corrmap1 = CorrMap(squeeze(PDIrm(:,:,3,:)), m1');
%         corrmap2 = CorrMap(squeeze(PDIrm(:,:,3,:)), m2');
%         corrmap3 = CorrMap(squeeze(PDIrm(:,:,3,:)), m3');
%         handles.BKmsk = (corrmap1<0.3).*(corrmap2<0.3).*(corrmap3<0.3);
%         figure; imagesc(handles.BKmsk); axis image;
        handles.PDIrm = PDIrm;


        PDI = medfilt1(PDI,3,[],4,'truncate');
        

        handles.PDI = PDI;
        guidata(hObject, handles);
%         handles.PDI = medfilt1(PDI,3,[],4,'truncate');% modified by Bingxue Liu, medfilt
        
        % clear PDI
        % PDI=handles.PDIHP0;
        % rplc=repmat(median(PDI,4)+std(movmean(PDI,15,4),1,4),[1,1,1,handles.nRpt]);
        % PDI(:,:,:,bulkIndex)=rplc(:,:,:,bulkIndex);
        % PDI=movmean(PDI,15,4);
        % handles.PDIHP=handles.PDIHP0;
        % handles.PDIHP(:,:,:,bulkIndex)=PDI(:,:,:,bulkIndex);
        
    
        
        switch handles.vUSorPDI
            case 'PDI'
                figure;
                subplot(2,1,1);
                sumPDIrm = squeeze(mean(mean(PDIrm(:,:,3,:),1),2));
                plot(squeeze(mean(mean(PDIrm(:,:,3,:),1),2)));
%                 hold on, plot(repmat(bulkThd,handles.nRpt))
                hold on; plot(handles.Stim*(max(sumPDIrm)-min(sumPDIrm))+min(sumPDIrm))
                title('mean-original PDI')
                subplot(2,1,2);
                sumPDI = squeeze(mean(mean(handles.PDI(:,:,3,:),1),2));
                plot(squeeze(mean(mean(handles.PDI(:,:,3,:),1),2)));
                hold on; plot(handles.Stim*(max(sumPDIrm)-min(sumPDIrm))+min(sumPDIrm))
                title('mean-corrected PDI')
            case 'vUS'
                % v correction
                handles.V=V;
                rplc=repmat(median(V,4)+1.2*std(movmean(V,15,4),1,4),[1,1,1,handles.nRpt]);
                V(:,:,:,bulkIndex)=rplc(:,:,:,bulkIndex);
                V=movmean(V,15,4);
                handles.V(:,:,:,bulkIndex)=V(:,:,:,bulkIndex);
                
                handles.Vz=Vz;
                rplc=repmat(median(Vz,4)+1.2*std(movmean(Vz,15,4),1,4),[1,1,1,handles.nRpt]);
                Vz(:,:,:,bulkIndex)=rplc(:,:,:,bulkIndex);
                Vz=movmean(Vz,15,4);
                handles.Vz(:,:,:,bulkIndex)=Vz(:,:,:,bulkIndex);
                
                handles.Vx=Vx;
                rplc=repmat(median(Vx,4)+1.2*std(movmean(Vx,15,4),1,4),[1,1,1,handles.nRpt]);
                Vx(:,:,:,bulkIndex)=rplc(:,:,:,bulkIndex);
                Vx=movmean(Vx,15,4);
                handles.Vx(:,:,:,bulkIndex)=Vx(:,:,:,bulkIndex);
                
                handles.Vcz=Vcz;
                rplc=repmat(median(Vcz,4)+1.2*std(movmean(Vcz,15,4),1,4),[1,1,1,handles.nRpt]);
                Vcz(:,:,:,bulkIndex)=rplc(:,:,:,bulkIndex);
                Vcz=movmean(Vcz,15,4);
                handles.Vcz(:,:,:,bulkIndex)=Vcz(:,:,:,bulkIndex);
                
                Fig=figure;
                set(Fig, 'Position',[400 300 1400 800]);
                subplot(2,3,1);
                plot(squeeze(mean(mean(handles.PDI0(:,:,3,:),1),2)));
                hold on, plot(repmat(bulkThd,handles.nRpt))
                title('mean-original PDI')
                
                subplot(2,3,2);
                
                plot(squeeze(mean(mean(handles.V0(:,:,2,:),1),2)));
                title('mean-original V')
                
                subplot(2,3,3);
                plot(squeeze(mean(mean(handles.Vz0(:,:,2,:),1),2)));
                title('mean-original Vz')
                
                subplot(2,3,4);
                plot(squeeze(mean(mean(handles.PDI(:,:,3,:),1),2)));
                title('mean-corrected PDI')
                
                subplot(2,3,5);
                plot(squeeze(mean(mean(handles.V(:,:,2,:),1),2)));
                title('mean-corrected V')
                
                subplot(2,3,6);
                plot(squeeze(mean(mean(handles.Vz(:,:,2,:),1),2)));
                title('mean-corrected Vz')
        end
    case 'NO'
         RMVslt = questdlg('Remove strong motion trial?', ...
            'Select', ...
            'NO', 'YES','Cancel');
        switch RMVslt
            case 'YES'
                prompt={'Index of trial removal', };
                name='Remove strong motion trial (0: no remove)';
                defaultvalue={'0'};
                numinput=inputdlg(prompt,name, 1, defaultvalue);
                RmvIndex=str2num(numinput{1});
                if RmvIndex(1)==0
                    PDI=handles.PDI0;
                    PDIHHP=handles.PDIHHP0;
                    V=handles.V0;
                    Vx=handles.Vx0;
                    Vz=handles.Vz0;
                    Vcz=handles.Vcz0;
                    handles.Stim=handles.Stim0;
                    handles.nTrial=handles.nStimTrial;
                    PDIrm = handles.PDI0;
                else
                    handles.nRpt= size(handles.PDI0,4); %modified by Bingxue
                    handles.nTrial=handles.nStimTrial-numel(RmvIndex);
                    Samples=ones(handles.nRpt,1);
                    stimSamples = ones(handles.nRpt,1);
                    for iTrial=1:numel(RmvIndex)
                        Samples((handles.nBase-nRest+(RmvIndex(iTrial)-1)*(handles.nStim+handles.nRecover)+1):...
                            (handles.nBase-nRest+(RmvIndex(iTrial))*(handles.nStim+handles.nRecover)))=0;  % modified by Bingxue, keep 5frames before stim trials
                        stimSamples((handles.nBase-nRest+(RmvIndex(iTrial)-1)*(handles.nStim+handles.nRecover)+1):...
                            (handles.nBase-nRest+(RmvIndex(iTrial))*(handles.nStim+handles.nRecover)))=0;
                    end
                    Samples=(Samples>0);
                    stimSamples=(stimSamples>0);
                    PDIrm=handles.PDI0(:,:,:,Samples); % PDIrm: remove motion trial
                    PDIHHP=handles.PDIHHP0(:,:,:,Samples);
%                     V=handles.V0(:,:,:,Samples);
%                     Vx=handles.Vx0(:,:,:,Samples);
%                     Vz=handles.Vz0(:,:,:,Samples);
%                     Vcz=handles.Vcz0(:,:,:,Samples); %comment for PDI processing
                    handles.Stim=handles.Stim0(stimSamples);
                    handles.Samples = Samples;
                end
            case 'NO'
                PDIrm=handles.PDI0;
                PDIHHP=handles.PDIHHP0;
                V=handles.V0;
                Vx=handles.Vx0;
                Vz=handles.Vz0;
                Vcz=handles.Vcz0;
                handles.Stim=handles.Stim0;
                handles.nTrial=handles.nStimTrial;
                handles.PDIrm = PDIrm;
        end
        
        PDI = PDIrm;
        handles.nRpt= size(PDI,4);
        handles.PDI = medfilt1(PDI,3,[],4,'truncate');% modified by Bingxue Liu, medfilt
        guidata(hObject, handles);
                figure;
                sumPDI = squeeze(mean(mean(handles.PDIHP0(:,:,3,:),1),2));
                subplot(2,1,1);
                plot(squeeze(mean(mean(handles.PDI0(:,:,3,:),1),2)));
                hold on, plot(handles.Stim*(max(sumPDI)-min(sumPDI))+min(sumPDI));
                title('mean-original PDI')
                subplot(2,1,2);
                plot(squeeze(mean(mean(handles.PDI(:,:,3,:),1),2)));
                hold on, plot(handles.Stim*(max(sumPDI)-min(sumPDI))+min(sumPDI));
                title('mean-corrected PDI')
                
             Accelslt = questdlg('Use accelerometer data?', ...
            'Select', ...
            'YES', 'NO','Cancel'); % modified by Bingxue
   switch Accelslt
       case 'YES'
     handles.Acc = [];
     handles.Acc = getAcceldata(handles.Stim0);
   end
%         handles.PDI=handles.PDI0;
%         handles.PDIHHP=handles.PDIHHP0;
%         handles.V=handles.V0;
%         handles.Vx=handles.Vx0;
%         handles.Vz=handles.Vz0;
%         handles.Vcz=handles.Vcz0;
%         handles.Stim=handles.Stim0;
%         handles.nTrial = handles.nStimTrial;
end
       
guidata(hObject, handles);


% --- Executes on button press in BTN_PLOT.
function BTN_PLOT_Callback(hObject, ~, handles)
% hObject    handle to BTN_PLOT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;
handles.nRest = 5;
nRest = handles.nRest;
guidata(hObject, handles);

addpath('D:\CODE\DataProcess\BingxueLiu\Code fUS DATA analysis')
addpath('D:\CODE\DataProcess\BingxueLiu\Code fUS DATA analysis\SubFunctions')
[VzCmap, VzCmapDn,VzCmapUp]=Colormaps_fUS;
prompt={'PDI(0),V(1),Vz(2),Vcz(3)','Fused (3), Up Flow (1), Down Flow (2)'};
name='File info';
defaultvalue={num2str(handles.DataSlt),num2str(handles.DnUpSlt)};
numinput=inputdlg(prompt,name, 1, defaultvalue);
handles.DataSlt=str2num(numinput{1}); % 
handles.DnUpSlt=str2num(numinput{2}); % 1: up, 2: down
handles.vUSroi=[];
handles.vzUSroi=[];
handles.vBBroi=[];
handles.PDIN=[];
handles.ImgShow=[];
handles.data=[];
if handles.DataSlt==0
    [nz,nx,~]=size(handles.PDI);
    %% PDI data and show data
    PDIslt = questdlg('Data Type?', ...
                         'Select', ...
                         'PDI', 'PDIHP', 'PDIHHP','Cancel');
    switch PDIslt
        case 'PDI'
            PDI=handles.PDI;
        case 'PDIHP'
            PDI=handles.PDIHP;
        case 'PDIHHP'
            PDI=handles.PDIHHP;
    end
    handles.PDIN(:,:,1,:)=1*(PDI(:,:,1,:)-repmat(min(min(PDI(:,:,1,:),[],1),[],2),[nz,nx,1,1]))./... % normalize to 0~1
        repmat(max(max(PDI(:,:,1,:),[],1),[],2)-min(min(PDI(:,:,1,:),[],1),[],2),[nz,nx,1,1]);
    handles.PDIN(:,:,2,:)=1*(PDI(:,:,2,:)-repmat(min(min(PDI(:,:,2,:),[],1),[],2),[nz,nx,1,1]))./...
        repmat(max(max(PDI(:,:,2,:),[],1),[],2)-min(min(PDI(:,:,2,:),[],1),[],2),[nz,nx,1,1]);
    handles.PDIN(:,:,3,:)=1*(PDI(:,:,3,:)-repmat(min(min(PDI(:,:,3,:),[],1),[],2),[nz,nx,1,1]))./...
        repmat(max(max(PDI(:,:,3,:),[],1),[],2)-min(min(PDI(:,:,3,:),[],1),[],2),[nz,nx,1,1]);
    
    handles.ImgShow=(mean(handles.PDIN(:,:,:,1:10),4)).^0.25;
    handles.data=PDI;   % process PDI data based on this
    handles.cAxis=[-3 0.8];
elseif handles.DataSlt==1
    handles.ImgShow=mean(handles.V0(:,:,:,1:10),4);
    handles.ImgShow(:,:,3)=max(abs(mean(handles.V0(:,:,:,1:10),4)),[],3);
    handles.cAxis=[-25 25];
    handles.data=handles.V;
elseif handles.DataSlt==2
    handles.ImgShow=mean(handles.Vz0(:,:,:,1:10),4);
    handles.ImgShow(:,:,3)=max(abs(mean(handles.Vz0(:,:,:,1:10),4)),[],3);
    handles.cAxis=[-20 20];
    handles.data=handles.Vz;
elseif handles.DataSlt==3
    handles.ImgShow=mean(handles.Vcz0(:,:,:,1:10),4);
    handles.ImgShow(:,:,3)=max(abs(mean(handles.Vcz0(:,:,:,1:10),4)),[],3);
    handles.cAxis=[-20 20];
    handles.data=handles.Vcz;
end
handles.cData=[];
%% figure plot
if handles.DnUpSlt==3 % all flow direction
    if handles.DataSlt==0 % PDI, all frequency
        ImgShowV=handles.ImgShow(:,:,3);
        axes(handles.axes1)
        imagesc(ImgShowV);axis image;
        colormap(gray);
        colorbar
        caxis([min(min(ImgShowV(:))*1.04, max(ImgShowV(:))*0.96) max(min(ImgShowV(:))*1.04, max(ImgShowV(:))*0.96)]);
        handles.cData=squeeze(handles.data(:,:,3,:));
    else 
        axes(handles.axes1)
        imagesc(handles.ImgShow(:,:,3));axis image;
        colormap(VzCmap);
        caxis(handles.cAxis);
        colorbar
        
%         axes(handles.axes1)
%         h1=imagesc(handles.ImgShow(:,:,1)); % up flow
%         colormap(VzCmap);
%         caxis(handles.cAxis);
%         colorbar
%         
%         hold on;
%         h2=imagesc(handles.ImgShow(:,:,2)); % down flow
%         set(h2,'AlphaData',1*double(abs(handles.V(:,:,2))>0.2))
%         colormap(VzCmap);
%         caxis(handles.cAxis);
%         colorbar
%         hold off
        handles.cData=squeeze(max(abs(handles.data(:,:,1:2,:)),[],3));
    end
else
    if handles.DataSlt==0 % PDI
        ImgShowV=handles.ImgShow(:,:,handles.DnUpSlt);
        axes(handles.axes1)
        imagesc(ImgShowV);axis image;
        colormap(gray);
        colorbar
        caxis([min(min(ImgShowV(:))*1.04, max(ImgShowV(:))*0.96) max(min(ImgShowV(:))*1.04, max(ImgShowV(:))*0.96)]);
        handles.cData=squeeze(handles.data(:,:,handles.DnUpSlt,:));
        title(['Up(1)/Down(2): ', num2str(handles.DnUpSlt)])
    else
        [VzCmap, VzCmapDn,VzCmapUp]=Colormaps_fUS;
        axes(handles.axes1)
        imagesc(handles.ImgShow(:,:,handles.DnUpSlt));axis image;
        axis tight
        caxis(handles.cAxis);
        colormap(VzCmap)
        colorbar;
        handles.cData=squeeze(handles.data(:,:,handles.DnUpSlt,:));
    end
end
%% remove large motion artifacts
% mAll=squeeze(mean(mean(handles.cData,1),2));
% MotionIndex=find(mAll>(mean(mAll)+1*std(mAll)));
% for iMI=1:length(MotionIndex)
%     handles.cData(:,:,MotionIndex(iMI))=handles.cData(:,:,MotionIndex(iMI)-1);
% end
%% plot negative-positive fused image
PLOTslt = questdlg('Plot Fused Image?', ...
                         'Select', ...
                         'YES', 'NO', 'Cancel','Cancel');
switch PLOTslt
    case 'YES'
        ImgShowV=mean(handles.V(:,:,:,1:1),4);
        ImgShowVz=mean(handles.Vz(:,:,:,1:1),4);
        Fig=figure;
        set(Fig,'Position',[400 400 1300 450])
        subplot(1,2,1)
        h1=imagesc(ImgShowV(:,:,1)); % up flow
        colormap(VzCmap);
        caxis([-30 30]);
        colorbar
        axis equal tight;
        hold on;
        h2=imagesc(ImgShowV(:,:,2)); % down flow
        set(h2,'AlphaData',1*double(abs(ImgShowV(:,:,2))>0.2))
        colormap(VzCmap);
        caxis([-30 30]);
        colorbar
        hold off
        axis equal tight;
        title('V')
        
        subplot(1,2,2)
        h1=imagesc(ImgShowVz(:,:,1)); % up flow
        colormap(VzCmap);
        caxis([-30 30]);
        colorbar
        axis equal tight;
        hold on;
        h2=imagesc(ImgShowVz(:,:,2)); % down flow
        set(h2,'AlphaData',1*double(abs(ImgShowVz(:,:,2))>0.2))
        colormap(VzCmap);
        caxis([-30 30]);
        colorbar
        hold off
        axis equal tight;
        title('Vz')
        
end
guidata(hObject, handles);
% --- Executes on button press in BTN_TrailAverageMap.
function BTN_TrailAverageMap_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_TrailAverageMap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;
nRest = handles.nRest;
addpath('D:\CODE\DataProcess\BingxueLiu\Code fUS DATA analysis')
addpath('D:\CODE\DataProcess\BingxueLiu\Code fUS DATA analysis\SubFunctions')
[VzCmap,VzCmapDn,VzCmapUp,pdiCmapUp]=Colormaps_fUS;
xCoor=[1:(handles.nStim+handles.nRecover)]*handles.tIntV;
stimTrial=handles.Stim(handles.nBase-nRest+1:handles.nBase-nRest+(handles.nStim+handles.nRecover));%5
[nz,nx,nDU,nt]= size(handles.PDI); %size(handles.V); modified by Bingxue 
prompt={'it, (5 frames baseline)'};
name='File info';
defaultvalue={'6'};
numinput=inputdlg(prompt,name, 1, defaultvalue);
handles.itSlt=str2num(numinput{1}); % 

prompt={'PDI (3), PDI+vUS(1)'};
name='Process Type';
defaultvalue={num2str(handles.DnUpSlt)};
numinput=inputdlg(prompt,name, 1, defaultvalue);
handles.DnUpSlt=str2num(numinput{1}); % modified by Bingxue

if handles.DnUpSlt==3 % all flow direction, pure PDI; modified by Bingxue
    PDI=handles.PDI(:,:,handles.DnUpSlt,:);    % use fused image; modified by Bingxue 
    PDItrial=reshape(PDI(:,:,handles.nBase-nRest+1:handles.nRpt-nRest),[nz,nx,(handles.nStim+handles.nRecover),handles.nTrial]);%5
    PDImean=mean(PDItrial,4);
    PDIstd=std(PDItrial,1,4);
    PDIratioTrial=PDItrial./mean(PDItrial(:,:,1:nRest,:),3)*100;%5
    PDIratioTrial(isnan(abs(PDIratioTrial)))=100;
    PDIratioTrial(abs(abs(PDIratioTrial)-100)<10)=100;
    PDImRatioTrial=mean(PDIratioTrial,4);
    PDIstdRatioTrial=std(PDIratioTrial,1,4);
    PDIsemRatioTrial=std(PDIratioTrial,1,4)./sqrt(handles.nTrial);
    
    Fig=figure;
    for it=4:9
        subplot(3, 3, it)
%         set(Fig,'Position',[300 400 1500 400]);
        imagesc(PDImRatioTrial(:,:,it));
        colormap(VzCmap);
        colorbar
        caxis([75 125]); %[70 130]
        axis equal tight
        title(['PDI ratil at it=', num2str(it)])
    end
else
    PDI=handles.PDI(:,:,handles.DnUpSlt,:);
    V=handles.V(:,:,handles.DnUpSlt,:);
    Vz=handles.Vz(:,:,handles.DnUpSlt,:);
    
    PDItrial=reshape(PDI(:,:,handles.nBase-nRest+1:handles.nRpt-nRest),[nz,nx,(handles.nStim+handles.nRecover),handles.nTrial]);
    PDImean=mean(PDItrial,4);
    PDIstd=std(PDItrial,1,4);
    PDIratioTrial=PDItrial./mean(PDItrial(:,:,1:nRest,:),3)*100;%5
    PDIratioTrial(isnan(abs(PDIratioTrial)))=100;
    PDIratioTrial(abs(abs(PDIratioTrial)-100)<10)=100;
    PDImRatioTrial=mean(PDIratioTrial,4);
    PDIstdRatioTrial=std(PDIratioTrial,1,4);
    PDIsemRatioTrial=std(PDIratioTrial,1,4)./sqrt(handles.nTrial);
    
    Vtrial=reshape(V(:,:,handles.nBase-nRest+1:handles.nRpt-nRest),[nz,nx,(handles.nStim+handles.nRecover),handles.nTrial]);
    Vmean=mean(Vtrial,4);
    Vstd=std(Vtrial,1,4);
    Vtrial(abs(Vtrial)<2)=1;
    VratioTrial=Vtrial./mean(Vtrial(:,:,1:nRest,:),3)*100;
    VratioTrial(isnan(abs(VratioTrial)))=100;
    VratioTrial(abs(abs(VratioTrial)-100)<10)=100;
%     VratioTrial(abs(abs(VratioTrial)-100)>50)=100;
    VmRatioTrial=mean(VratioTrial,4);
    VstdRatioTrial=std(VratioTrial,1,4);
    VsemRatioTrial=std(VratioTrial,1,4)./sqrt(handles.nTrial);
    
    Vztrial=reshape(Vz(:,:,handles.nBase-nRest+1:handles.nRpt-nRest),[nz,nx,(handles.nStim+handles.nRecover),handles.nTrial]);
    Vzmean=mean(Vztrial,4);
    Vzstd=std(Vztrial,1,4);
    Vztrial(abs(Vztrial)<3)=1;
    VzratioTrial=Vztrial./mean(Vztrial(:,:,1:nRest,:),3)*100;
    VzratioTrial(isnan(abs(VzratioTrial)))=100;
    VzratioTrial(abs(abs(VzratioTrial)-100)<10)=100;
    VzmRatioTrial=mean(VzratioTrial,4);
    VzstdRatioTrial=std(VzratioTrial,1,4);
    VzsemRatioTrial=std(VzratioTrial,1,4)./sqrt(handles.nTrial);
    
    Fig=figure;
    set(Fig,'Position',[150 200 750 200]);
    subplot(1,3,1);
    imagesc(PDImRatioTrial(:,:,handles.itSlt));
    colormap(VzCmap);
    colorbar
    caxis([70 130])
    axis equal tight
    title(['PDI ratil at it=', num2str(handles.itSlt)])
    
    subplot(1,3,2);
    imagesc(VmRatioTrial(:,:,handles.itSlt));
    colormap(VzCmap);
    colorbar
    caxis([70 130])
    axis equal tight
    title(['V ratil at it=', num2str(handles.itSlt)])
    
    subplot(1,3,3);
    imagesc(VzmRatioTrial(:,:,handles.itSlt));
    colormap(VzCmap);
    colorbar
    caxis([70 130])
    axis equal tight
    title(['Vz ratil at it=', num2str(handles.itSlt)])
    
    for it=4:9
        Fig=figure;
%         set(Fig,'Position',[150 200 750 200]); %[300 400 1500 400]
        subplot(1,3,1);
        imagesc(PDImRatioTrial(:,:,it));
        colormap(VzCmap);
        colorbar
        caxis([70 130])
        axis equal tight
        title(['PDI ratil at it=', num2str(it)])
        
        subplot(1,3,2);
        imagesc(VmRatioTrial(:,:,it));
        colormap(VzCmap);
        colorbar
        caxis([70 130])
        axis equal tight
        title(['V ratil at it=', num2str(it)])
        
        subplot(1,3,3);
        imagesc(VzmRatioTrial(:,:,it));
        colormap(VzCmap);
        colorbar
        caxis([70 130])
        axis equal tight
        title(['Vz ratil at it=', num2str(it)])
    end
end

% --- Executes on button press in BTN_ROI_tCourse.
function BTN_ROI_tCourse_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_ROI_tCourse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;
addpath('D:\CODE\DataProcess\BingxueLiu\Code fUS DATA analysis')
addpath('D:\CODE\DataProcess\BingxueLiu\Code fUS DATA analysis\SubFunctions')
nRest = handles.nRest;
%% select ROI
prompt={'number of ROIs'};
name='ROI info';
defaultvalue={num2str(handles.nROI)};
numinput=inputdlg(prompt,name, 1, defaultvalue);
handles.nROI=str2num(numinput{1}); % 

ButtonName = questdlg('Select New ROI?', ...
        'Select', ...
        'Yes', 'No', 'Cancel');
switch ButtonName
    case 'Yes'
        handles.mData=zeros(handles.nROI,handles.nRpt);
        handles.stdData=zeros(handles.nROI,handles.nRpt);
        handles.semData=zeros(handles.nROI,handles.nRpt);
        handles.BW=[];            % ROI selected manually
        handles.mData=[];         % mean value data in chosed ROI across all time
        handles.stdData=[];
        handles.semData=[];
        for iROI=1:handles.nROI
            h=msgbox(['Plese select the ', num2str(iROI), ' ROI...']);
            pause(0.5)
            delete(h)
            [xROI, zROI]=ginput(6);
            xROI=floor(xROI);
            zROI=floor(zROI);
            handles.BW(:,:,iROI)=roipoly(handles.data(:,:,1,1),xROI,zROI);
        end
    case 'No'
                handles.mData=zeros(handles.nROI,handles.nRpt);
        handles.stdData=zeros(handles.nROI,handles.nRpt);
        handles.semData=zeros(handles.nROI,handles.nRpt);
%         handles.BW=[];            % ROI selected manually
        handles.mData=[];         % mean value data in chosed ROI across all time
        handles.stdData=[];
        handles.semData=[];
end

%% High pass filter to remove low temperal fluctuations
%        PDI = handles.PDI;
        handles.cData0 = handles.cData; % save original data w.o. high pass filter
        [~,~,nt] = size(handles.cData);
        [B,A]=butter(4,0.0052/(1/2),'high');    %coefficients for the high pass filter
        nt0 = nt-50; % padding time points
        PDI1(:,:,nt0+1:nt0+nt)= handles.cData;%squeeze(PDI(:,:,3,:));
        PDI1(:,:,1:nt0)=flip(PDI1(:,:,nt0+1:nt0+nt0),3);
        PDI2=filter(B,A,PDI1,[],3);    % 
        handles.cData =PDI2(:,:,nt0+1:end); % High pass filtered
      guidata(hObject, handles); 
%% Normalize trial mean 0 std 1
% handles.cData = (handles.cData-repmat(mean(handles.cData,3),[1,1,nt]))...
%     ./(repmat(sqrt(sum((handles.cData-repmat(mean(handles.cData,3),[1,1,nt])).^2,3)),[1,1,nt]));
% handles.cData = 100*handles.cData+100;
%%
for iROI=1:handles.nROI
    %% Value in the selected ROI
    % baseline pixel index for the ROI t course calculation
    DataBase=handles.BW(:,:,iROI).*mean(handles.cData(:,:,1:handles.nBase),3);
    Locs=abs(DataBase)~=0;
    for iRpt=handles.startRpt:handles.startRpt+handles.nRpt-1
        Data=handles.BW(:,:,iROI).*handles.cData(:,:,iRpt);
        roi=Data(Locs>0);
        handles.mData{iROI}(iRpt)=mean(roi(:)); 
        handles.stdData{iROI}(iRpt)=std(roi(:));
        handles.semData{iROI}(iRpt)=std(roi(:))/sqrt(numel(roi));
        % calculate original data within ROI (no HP)
        Data0 = handles.BW(:,:,iROI).*handles.cData0(:,:,iRpt);
        roi0 = Data0(Locs>0);
        handles.mData0{iROI}(iRpt) = mean(roi0(:));
        handles.stdData0{iROI}(iRpt)=std(roi0(:));
        handles.semData0{iROI}(iRpt)=std(roi0(:))/sqrt(numel(roi0));
    end    
end
guidata(hObject, handles);
for iROI=1:handles.nROI
    %% time course plot for each ROI
    xCoor=[handles.startRpt:handles.startRpt+handles.nRpt-1]*handles.tIntV;
    Color.Shade='b';
    Color.ShadeAlpha=0.3;
    Color.Line='b';
    handles.tCourseFhandle=figure;
     set(handles.tCourseFhandle,'Position',[200 400 1000 200])
    ShadedErrorbar(handles.mData{iROI},handles.semData{iROI},xCoor,Color);
    xlim([xCoor(1),xCoor(end)])
    xlabel('t [s]');
    if handles.DataSlt==0 % PDI
        ylabel('PDI [a.u.]');
    elseif handles.DataSlt==1 % V
        ylabel('V [mm/s]');
    elseif handles.DataSlt==2 % Vz
        ylabel('V_z [mm/s]');
    elseif handles.DataSlt==3 % Vcz
        ylabel('V_c_z [mm/s]');
    end
    hold on;
    plot(xCoor, handles.Stim*abs(max(abs(handles.mData{iROI}))*0.95-min(abs(handles.mData{iROI}))*1.0)+min(abs(handles.mData{iROI}))*1.0,'r','LineWidth',0.01);alpha(0.2);
%     area(1:handles.nRpt, (handles.Stim*abs(max(abs(handles.mData{iROI}))*0.95-min(abs(handles.mData{iROI}))*1.0)+min(abs(handles.mData{iROI}))*1.0).*ones(1,handles.nRpt),'FaceColor','r','EdgeColor','none'); alpha(0.1);
%     ylim([min((handles.mData{iROI})),max((handles.mData{iROI}))])
    grid on
    
    %% Plot trial averaged time course
    % stimulation parttern
    xCoor=[1:(handles.nStim+handles.nRecover)]*handles.tIntV;
    stimTrial=handles.Stim(handles.nBase-nRest+1:handles.nBase-nRest+(handles.nStim+handles.nRecover));
    % obtain trial data
    Trial=reshape(handles.mData{iROI}(handles.nBase-nRest+1:handles.nRpt-nRest),[(handles.nStim+handles.nRecover),handles.nTrial]);
    Trial0 = reshape(handles.mData0{iROI}(handles.nBase-nRest+1:handles.nRpt-nRest),[(handles.nStim+handles.nRecover),handles.nTrial]);
    mTrial(iROI,:)=median(Trial,2)';% mean
    stdTrial(iROI,:)=std(Trial,1,2)';
    semTrial(iROI,:)=(std(Trial,1,2)./sqrt(handles.nTrial))';
    RatioTrial=(Trial-repmat(mean(Trial(1:nRest,:),1),[size(Trial,1),1]))./abs(mean(Trial0(1:nRest,:),1));% mean
    
    mRatioTrial(iROI,:)=median(RatioTrial,2)';   % averaged time course, relative change to 5 basline before stimulus
    stdRatioTrial(iROI,:)=std(RatioTrial,1,2)';
    semRatioTrial(iROI,:)=(std(RatioTrial,1,2)./sqrt(handles.nTrial))';
end
COLOR=[1 0 0;
    0.08 0.17 0.55
    0.31 0.31 0.31];
figure;
for iROI=1:handles.nROI
    Color.Shade=COLOR(iROI,:);
    Color.ShadeAlpha=0.3;
    Color.Line=COLOR(iROI,:);
    %% avearged relative change
    hold on;     ShadedErrorbar(mRatioTrial(iROI,:)*100,semRatioTrial(iROI,:)*100,xCoor,Color);
    title('Trial averaged response')
    ylabel('Relative Change [%]')
    xlabel('t [s]')
end
% hold on, plot(xCoor, stimTrial*(max(abs(mRatioTrial(:)))*100-100)+100,'b')
hold on, plot(xCoor, stimTrial*(max(mRatioTrial(:)))*100,'b'); 

set(gca, 'YGrid', 'on', 'XGrid', 'off')


%% generate Video
PLOTslt = questdlg('Plot CourseVideo?', ...
                         'Select', ...
                         'YES', 'NO', 'Cancel');
switch PLOTslt
    case 'YES'
        [VzCmap, VzCmapDn,VzCmapUp, pdiCmapUp, PhtmCmap]=Colormaps_fUS;
%         Data=squeeze(handles.Vcz(:,:,1,:));
%         Data=squeeze(handles.V(:,:,1,:));
%         Data=squeeze(handles.PDI(:,:,1,:));
        nBase=handles.nBase;
        nPrd=handles.nStim+handles.nRecover;
        
        Data=handles.cData(:,:,:);
        Data=movmean(Data,5,3);  
        ImgBase=squeeze(mean(abs(Data(:,:,1:nBase)),3));
%                 ImgBase=squeeze(mean(abs(Data(:,:,1:end)),3)); % for resting state
%         ImgBase(ImgBase<3)=3;
        %         ImgRLTV=movmean(squeeze(abs(Data(:,:,:))),3,3)./ImgBase;
        ImgRLTV=squeeze(abs(Data(:,:,:)))./ImgBase;
%         ImgRLTV(ImgRLTV>2.5)=1;
        ImgRLTV=movmean(ImgRLTV,5,3);
        mImgRLTV=mean(abs(ImgRLTV),3);
        stdImgRLTV=std(abs(ImgRLTV),1,3);
        
        ImgVBK=squeeze(mean(abs(Data(:,:,1:nBase)),3));
        vBKmsk=zeros(size(ImgVBK));
        vBKmsk(abs(ImgVBK)>2.5)=1;  % ???? modified by Bingxue
%         % 0503 left whisker
%         load('D:\OneDrive\Work\PROJ - FUS\PROJ - US velocimetry\Nature Communication\Figure4-Whisker\0503-CP6-AWAKE-LEFT Whisker\BKmsk.mat')
%         BKmsk1=BKmsk.*vBKmsk;
        
        Background = questdlg('Grey scale background selection?', ...
            'Select', ...
            'PDI', 'ULM', 'Cancel');
        switch Background
            case 'PDI'
                Noise=(mean(handles.eqNoise,3)).^0.3;
                
                if ~isfield(handles,'PDI')
                    [FileName,FilePath]=uigetfile(handles.DefPath);
                    load([FilePath,FileName]);
                    handles.PDI=PDI;
                    eqPDI=mean(squeeze(PDI(:,:,3,:)),3)./Noise;
                    ImgBKGDPDI=log(eqPDI);
                else
                    eqPDI=mean(squeeze(handles.PDI(:,:,3,:)),3)./Noise;
                    ImgBKGDPDI=log(eqPDI);
                end 
                ImgBK=(ImgBKGDPDI-min(ImgBKGDPDI(:)))/(median(ImgBKGDPDI(:))+4*std(ImgBKGDPDI(:))-min(ImgBKGDPDI(:)));
                % 0503 left whisker
%                 load('D:\OneDrive\Work\PROJ - FUS\PROJ - US velocimetry\Nature Communication\Figure4-Whisker\0503-CP6-AWAKE-LEFT Whisker\PDIBK.mat')
%                 ImgBK=PDI;
            case 'ULM'
                %% use ULM as background, PDI correlation map plot
                [nz,nx,~]=size(handles.V);
                if ~isfield(handles,'BB')
                    [FileName,FilePath]=uigetfile(handles.DefPath);
                    load([FilePath,FileName]);
                    handles.BB=BB;
                end
                ImgBK=imresize(abs((handles.BB(:,:,3)).^0.4),[nz,nx]); % use ULM as the background
                ImgBK=(ImgBK-min(ImgBK(:)))/(max(ImgBK(:))-min(ImgBK(:)));
        end
        [nz,nx,nt]=size(ImgRLTV);
        %% make gif
        
        for iSeg=1:ceil(nt/10)
            fig=figure;
            set(gcf,'color','white')
            set(fig,'Position',[300 150 800 500]) % [600 400 800 500] modified by Bingxue
            for it=(iSeg-1)*10+1:min(iSeg*10,nt) %floor(nFrame/20):nFrame
                for i = 1: 1+isfield(handles, 'Acc')          % modified by Bingxue
                    ax(i) = subplot(1+isfield(handles, 'Acc'),1,i);
                end                                         
                subplot(ax(1))                                 % end
                hAxes1 = ax(1); %ax(1);%axes;
               hAxes1.Position = [0.1300, 0.5838, 0.6974, 0.3412];
                imagesc(ImgBK);
%                 caxis(hAxes1,[(1+0.1*sign(min(ImgBK(:))))*min(ImgBK(:)) max(ImgBK(:))*0.9]);
                caxis(hAxes1,[0.05 1])
                colormap(hAxes1,gray)
                colorbar('Ticks',[-10 10], 'TickLabels',{[],[]} );
                axis equal tight
                axis(hAxes1,'off')
                
                hold on;
                hAxes2 = axes;%axes;
                hAxes2.Position = [0.1300, 0.5838, 0.6974, 0.3412];
                RLTVMap=ImgRLTV(:,:,it)*100;
                h2=imagesc(RLTVMap);
%                 set(h2,'AlphaData',abs(RLTVMap)/(max(abs(RLTVMap(:)))/4).*((abs(RLTVMap-100))>10).*((abs(RLTVMap-100))<70).*BKmsk1.*(ImgBK>0.18).*(1-(mImgRLTV>1.15)))
%                 set(h2,'AlphaData',abs(RLTVMap)/(max(abs(RLTVMap(:)))/5).*((abs(RLTVMap-100))>10).*((abs(RLTVMap-100))<70).*BKmsk1.*(1-(mImgRLTV>1.18)).*(1-(stdImgRLTV>0.35)))
%                 set(h2,'AlphaData',abs(RLTVMap)/(min(max(abs(RLTVMap(:))),200)/5).*((abs(RLTVMap-100))>10).*((abs(RLTVMap-100))<70).*vBKmsk.*(1-(stdImgRLTV>0.25)))
                set(h2,'AlphaData',abs(RLTVMap)/(min(max(abs(RLTVMap(:))),200)/5).*((abs(RLTVMap-100))>10).*((abs(RLTVMap-100))<70).*(1-(stdImgRLTV>0.25)))% 0.4 % two Binary Maskes: 1:only show variation larger than 1/10;2:filter out big noise std larger than 0.4; Bingxue  
%                 set(h2,'AlphaData',abs(RLTVMap)/(max(abs(RLTVMap(:)))/5).*((abs(RLTVMap-100))>10).*vBKmsk.*(ImgBK>0.18))
                % colormap(hAxes2,VzCmap);
                % caxis(hAxes2,[-1 1])
                colormap(hAxes2,VzCmap);
                caxis(hAxes2,[50 150])
                colorbar
                hold off
                axis equal tight
                axis(hAxes2,'off')
                
                linkaxes([hAxes1,hAxes2])
                if rem(it-nBase-1,nPrd)<0
                    textshow=text(93,135,['BaseLine',', t=',num2str(it/handles.P.PDIFR),' s']);
                    textshow.Color='white';
                    textshow.FontSize=13;
                elseif rem(it-nBase-1,nPrd)<handles.nStim && rem(it-nBase-1,nPrd)>=0
                    textshow=text(90,135,['Stim. ON',', t=',num2str(it/handles.P.PDIFR),' s']);
                    textshow.Color='red';
                    textshow.FontSize=15;
                else
                    textshow=text(93,135,['Stim. OFF',', t=',num2str(it/handles.P.PDIFR),' s']);
                    textshow.Color='white';
                    textshow.FontSize=13;
                end
                axis off
                if isfield(handles, 'Acc') % modified by Bingxue
                subplot(ax(2));
%                 daccRaw = [0;diff(handles.Acc.accRaw)];
%                 plot((1:handles.Acc.nbin*1e3)/1e3, daccRaw);
%                 hold on
%                 if abs(daccRaw(it*1e3))> 1.3*1e4*abs(mean(daccRaw(:)))
%                     plot(it,daccRaw(it*1e3),'ro');
%                 else
%                     plot(it,daccRaw(it*1e3),'ko');
%                 end   
                 Stim = downsample(handles.Stim0, size(handles.Stim0,2)/handles.Acc.nbin);
                plot(1:handles.Acc.nbin, handles.Acc.accBinVar);
                hold on; plot(handles.Acc.accSumVar*ones(handles.Acc.nbin,1),'-.r');
                hold on; area(1:handles.Acc.nbin, Stim*max(handles.Acc.accBinVar),'FaceColor','r','EdgeColor','none'); alpha(0.25);
                hold on; 
                if ismember(it, handles.Acc.ibin)
                    plot(it,handles.Acc.accBinVar(it),'ro');
                else
                    plot(it,handles.Acc.accBinVar(it),'ko');
                end
                xlabel('Time(Sec)');
                ylabel('Variance of Accel Data');
                end                        % end
                hold off
                drawnow; pause(0.01); % (0.1)
                frames(it)=getframe(gcf);
            end
            close (fig)
            pause(0.2);
        end
        if handles.DataSlt==0
            FileName='PDI-RLTV-';
        elseif handles.DataSlt==1
            FileName='V-RLTV-';
        elseif handles.DataSlt==2
            FileName='Vz-RLTV-';
        elseif handles.DataSlt==3
            FileName='Vcz-RLTV-';
        end
        if handles.DnUpSlt==3
            FileName=[FileName,'ALL'];
        elseif handles.DnUpSlt==1
            FileName=[FileName,'Up'];
        elseif handles.DnUpSlt==2
            FileName=[FileName,'Dn'];
        end
        outfileGIF=[handles.DefPath,FileName,'.gif'];
        newVid = VideoWriter([handles.DefPath,FileName, '.avi']);
        delayTime=0.5; 
        newVid.FrameRate = 1/delayTime;
        newVid.Quality = 100;
        open(newVid)
        for ifile=1:nt
            % save GIF
            im = frame2im(frames(ifile));
            [imind,cm] = rgb2ind(im,256);
            if ifile==1
                imwrite(imind,cm,outfileGIF,'gif','DelayTime',delayTime,'loopcount',inf);
            else
                imwrite(imind,cm,outfileGIF,'gif','DelayTime',delayTime,'writemode','append');
            end
            % save AVI video
            writeVideo(newVid,frames(ifile).cdata)%within the for loop saving one frame at a tim
        end
        close(newVid)
end
% Color.Shade='m';
% Color.ShadeAlpha=0.3;
% Color.Line='m';
% ShadedErrorbar(abs(handles.Vzmean),handles.Vzsem,xCoor,Color);

% axis tight
% grid on
% ylim([0 20]);

guidata(hObject, handles);

% --- Executes on button press in BTN_responseFunction.
function BTN_responseFunction_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_responseFunction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
nRest = handles.nRest;
h=msgbox(['Plese select the responding ROI...']);
pause(0.5)
delete(h)
[xROI, zROI]=ginput(6);
xROI=floor(xROI);
zROI=floor(zROI);
handles.BW=[];
handles.BW(:,:,1)=roipoly(handles.data(:,:,1,1),xROI,zROI);
for iROI=1:1
    %% Value in the selected ROI
    % baseline pixel index for the ROI t course calculation
    DataBase=handles.BW(:,:,iROI).*mean(handles.cData(:,:,1:handles.nBase),3);
    Locs=abs(DataBase)~=0;
    for iRpt=handles.startRpt:handles.startRpt+handles.nRpt-1
        Data=handles.BW(:,:,iROI).*handles.cData(:,:,iRpt);
        roi=Data(Locs>0);
        handles.mData{iROI}(iRpt)=mean(roi(:));
        handles.stdData{iROI}(iRpt)=std(roi(:));
        handles.semData{iROI}(iRpt)=std(roi(:))/sqrt(numel(roi));
                % calculate original data within ROI (no HP)
        Data0 = handles.BW(:,:,iROI).*handles.cData0(:,:,iRpt);
        roi0 = Data0(Locs>0);
        handles.mData0{iROI}(iRpt) = mean(roi0(:));
        handles.stdData0{iROI}(iRpt)=std(roi0(:));
        handles.semData0{iROI}(iRpt)=std(roi0(:))/sqrt(numel(roi0));
    end    
end
guidata(hObject, handles);
for iROI=1:1
    %% time course plot for each ROI
    xCoor=[handles.startRpt:handles.startRpt+handles.nRpt-1]*handles.tIntV;
    
    %% trial averaged time course
    xCoor=[1:(handles.nStim+handles.nRecover)]*handles.tIntV;
    stimTrial=handles.Stim(handles.nBase-nRest+1:handles.nBase-nRest+(handles.nStim+handles.nRecover));% 5
    
    Trial=reshape(handles.mData{iROI}(handles.nBase-nRest+1:handles.nRpt-nRest),[(handles.nStim+handles.nRecover),handles.nTrial]);% 10
    Trial0=reshape(handles.mData0{iROI}(handles.nBase-nRest+1:handles.nRpt-nRest),[(handles.nStim+handles.nRecover),handles.nTrial]);% 10
    mTrial(iROI,:)=median(Trial,2)'; % mean
    stdTrial(iROI,:)=std(Trial,1,2)';
    semTrial(iROI,:)=(std(Trial,1,2)./sqrt(handles.nTrial))';
    RatioTrial=(Trial-repmat(mean(Trial(1:nRest,:),1),[size(Trial,1),1]))./mean(Trial0(1:nRest,:),1);% mean);% 5
    mRatioTrial(iROI,:)=median(RatioTrial,2)'; % mean
    stdRatioTrial(iROI,:)=std(RatioTrial,1,2)';
    semRatioTrial(iROI,:)=(std(RatioTrial,1,2)./sqrt(handles.nTrial))';
end

%% obtain response function
mResponse=mRatioTrial*100;
handles.mResponse = mResponse;
HRFmov=movmean(mResponse,5);% averaged ratio trial after movemean 
lb = [1     0.1     0.1     xCoor(5)]; % [a, b, c, t0]
ub = [200   10      2     xCoor(8)]; 
fun = fittype('a*(t-t0).^2.*exp(-b*abs(t-t0).^c).*normcdf(t,t0,0.006)','indep','t');
Fopts = fitoptions(fun);
Fopts.Lower = lb;
Fopts.Upper = ub;
Fopts.Display='off';
Fopts.TolFun=1e-6;
Fopts.Robust='LAR';
Fopts.StartPoint=[150, 2, 1, xCoor(5)];
est = fit(xCoor',mResponse',fun,Fopts);
HRF=fun(est.a,est.b,est.c,est.t0,xCoor);
R_HRF=(1-sum(abs(mResponse-HRF).^2)./sum(abs((mResponse)-mean(mResponse)).^2));
R_HRFmov=(1-sum(abs(mResponse-HRFmov).^2)./sum(abs((mResponse)-mean(mResponse)).^2));
COLOR=[1 0 0;
    0.08 0.17 0.55
    0.31 0.31 0.31];
tCoor=xCoor-xCoor(5);
figure;
for iROI=1:1
    Color.Shade=COLOR(iROI,:);
    Color.ShadeAlpha=0.3;
    Color.Line=COLOR(iROI,:);
    %% avearged relative change
    hold on;     ShadedErrorbar(mRatioTrial(iROI,:)*100,semRatioTrial(iROI,:)*100,tCoor,Color); 
    title('Trial averaged response')
    ylabel('Relative Change [%]')
    xlabel('t [s]')
end
% hold on, plot(tCoor, stimTrial*(max(abs(mRatioTrial(:)))*100-100)+100,'b')
hold on, plot(tCoor, stimTrial*(max(mRatioTrial(:)))*100,'b');
set(gca, 'YGrid', 'on', 'XGrid', 'off')
hold on, plot(tCoor,HRF,'g','LineWidth',2)
title(['R=',num2str(R_HRF)])
handles.HRF=HRF;
handles.HRFmov=HRFmov;
% hold on, plot(xCoor,HRFmov)
% [muhat,muci] = gamfit(mRatioTrial-min(mRatioTrial));
% HRF=gampdf(xCoor,muhat(1),muhat(2));
% figure,plot(xCoor,mRatioTrial-min(mRatioTrial),'.');
% hold on, plot(xCoor,HRF)
guidata(hObject, handles);

% --- Executes on button press in BTN_StimCorrCoef.
function BTN_StimCorrCoef_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_StimCorrCoef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;
nRest = handles.nRest;
addpath('D:\CODE\DataProcess\BingxueLiu\Code fUS DATA analysis')
addpath('D:\CODE\DataProcess\BingxueLiu\Code fUS DATA analysis\SubFunctions')
[VzCmap, VzCmapDn,VzCmapUp]=Colormaps_fUS;

%%%%%%%%%%%% averaged trial for correlation coefficinet calculation %%%%%%%

Stim=handles.HRF/100;  % use HRF as stim parttern
TrialIndex=ones(handles.nBase-nRest+(handles.nStim+handles.nRecover)*handles.nTrial,1);
TrialIndex(1:handles.nBase-nRest)=0;
% PDI0=squeeze(handles.PDI(:,:,handles.DnUpSlt,TrialIndex>0));
PDI0=squeeze(handles.cData(:,:,TrialIndex>0)); % PDI0: use HP cData 
[nz,nx,nt]=size(PDI0);
PDI1=median(reshape(PDI0,[nz,nx,handles.nStim+handles.nRecover,handles.nTrial]),4); % median value of all trials
BB=ones(3,3,3);
BB(2,2,2)=36;
BB=BB/62;
PDI=convn(PDI1,BB,'same'); % PDI: after spatial smoothing

%% PDI Background
Noise=(mean(handles.eqNoise,3)).^0.35;
eqPDI=mean(handles.cData0,3)./Noise; % use cData0(no HP) as background
ImgBKGDPDI=(log(eqPDI)).^0.3;
ImgBKGDPDI=(ImgBKGDPDI-min(ImgBKGDPDI(:)))/(median(ImgBKGDPDI(:))+4*std(ImgBKGDPDI(:))-min(ImgBKGDPDI(:)));
PLOTslt = questdlg('New Brain Region Outline?', ...
    'Select', ...
    'YES', 'NO', 'Cancel');
switch PLOTslt
    case 'YES'
        handles.BWBrain=[];
        figure;
        hAxes1 = axes;
        imagesc(ImgBKGDPDI);axis image;
        caxis(hAxes1,[0.3 1]);
        colormap(hAxes1,gray)
        [xROI, zROI]=ginput(20);
        xROI=floor(xROI);
        zROI=floor(zROI);
        handles.BWBrain(:,:,1)=roipoly(handles.data(:,:,1,1),xROI,zROI);
end
if ~isfield(handles,'BWBrain')
%     figure;
%     hAxes1 = axes;
%     imagesc(ImgBKGDPDI);
%     caxis(hAxes1,[0.3 1]);
%     colormap(hAxes1,gray)
%     [xROI, zROI]=ginput(20);
%     xROI=floor(xROI);
%     zROI=floor(zROI);
%     handles.BWBrain(:,:,1)=roipoly(handles.data(:,:,1,1),xROI,zROI);
handles.BWBrain(:,:,1) = ones(size(ImgBKGDPDI));
end
%
Stim1 = repmat(Stim, [1, handles.nTrial]); % Stim1 is not trial averaged stimulus

stim=zeros(30,1);
stim(6:10,:)=1;                         % activity is a square window
hrf = hemodynamicResponse(1,[1.5 10 0.5 1 20 0 16]); 
stim=filter(hrf,1,stim);                 % stim is hrf 

coefMapPDI=CoorCoeffMap(PDI0, Stim1,0);  % non trail averaged; for trial averaged, use (PDI1, Stim)
coefMapPDI=CoorCoeffMap(PDI1, Stim,0);
coefMapPDI=CoorCoeffMap(PDI1, stim',0);
handles.coefMapPDI = coefMapPDI;
%% create mask remove specular reflection artifacts
PDIrm = handles.PDIrm;
% %         m1 = squeeze(squeeze(PDIrm(150,100,3,:))); % 150 120 (0528 data)
% %         m2 = squeeze(squeeze(PDIrm(94, 12, 3, :))); % 100 155 (0528 data)
% %         m3 = squeeze(squeeze(PDIrm(44, 154, 3, :)));
% %         m4 = squeeze(squeeze(PDIrm(47, 10, 3, :)));
% %         m5 = squeeze(squeeze(PDIrm(10, 130, 3, :))); %%%% 0520 data
        m1 = squeeze(squeeze(PDIrm(8,96,3,:)));
        m2 = squeeze(squeeze(PDIrm(16, 84, 3, :))); 
        m3 = squeeze(squeeze(PDIrm(45, 168, 3, :)));
        m4 = squeeze(squeeze(PDIrm(38, 2, 3, :)));
        m5 = squeeze(squeeze(PDIrm(3, 5, 3, :))); %%%% 0609 data
        m1 = squeeze(squeeze(PDIrm(14,104,3,:)));
        m2 = squeeze(squeeze(PDIrm(17, 88, 3, :))); 
        m3 = squeeze(squeeze(PDIrm(143, 28, 3, :)));
        m4 = squeeze(squeeze(PDIrm(49, 167, 3, :)));
        m5 = squeeze(squeeze(PDIrm(149, 86, 3, :))); %%%% 0622 data

        corrmap1 = CorrMap(squeeze(PDIrm(:,:,3,:)), m1');
        corrmap2 = CorrMap(squeeze(PDIrm(:,:,3,:)), m2');
        corrmap3 = CorrMap(squeeze(PDIrm(:,:,3,:)), m3');
        corrmap4 = CorrMap(squeeze(PDIrm(:,:,3,:)), m5');
        corrmap5 = CorrMap(squeeze(PDIrm(:,:,3,:)), m5');
                figure; subplot(511);plot(m1);subplot(512);plot(m2);subplot(513);plot(m3);subplot(514);plot(m4);subplot(515);plot(m5)
        figure; subplot(231);imagesc(corrmap1);axis image;subplot(232);imagesc(corrmap2);axis image;
        subplot(233);imagesc(corrmap3);axis image; subplot(234);imagesc(corrmap4);axis image;
        subplot(235);imagesc(corrmap5);axis image;
%         handles.BKmsk =
%         (corrmap1<0.3).*(corrmap2<0.3).*(corrmap3<0.3).*(corrmap4<0.3).*(corrmap5<0.3);%
%         for 0520 data
        handles.BKmsk = (corrmap1<0.3).*(corrmap2<0.3).*(corrmap3<0.3).*(corrmap4<0.6).*(corrmap5<0.6);       
figure; imagesc(handles.BKmsk); axis image;
figure; imagesc(coefMapPDI.*handles.BKmsk); axis image; caxis([0.2,0.8])
handles.BKmsk = 1;
figure;
hAxes1 = axes;
imagesc(abs(ImgBKGDPDI));
%             caxis(hAxes1,[mean(ImgBKGDPDI(:))-1.5*std(ImgBKGDPDI(:)) mean(ImgBKGDPDI(:))+4*std(ImgBKGDPDI(:))]);
caxis(hAxes1,[0.3 1]);
colormap(hAxes1,gray)
colorbar('Ticks',[-10 10], 'TickLabels',{[],[]} );
axis equal tight
%             BKmsk=ImgBKGDPDI+min(ImgBKGDPDI(:));
%             BKmsk=BKmsk/(mean(BKmsk(:))+1*std(BKmsk(:)));
%             BKmsk(BKmsk>1)=1;

% BKmsk=zeros(size(ImgBKGDPDI));
% BKmsk(abs(ImgBKGDPDI)>(mean(ImgBKGDPDI(:))-0.1*std(ImgBKGDPDI(:))))=1;

%                 BKmsk=zeros(size(ImgBKGDPDI));
%   stdPDI = squeeze(stdfilt(handles.PDI(:,:,3,:),true(3))/sqrt(numel(true(3))));
%   mstdPDI = mean(stdPDI,3);
%   BKmsk = mstdPDI>10*median(mstdPDI(:));
%   figure; imagesc(BKmsk);axis image;
%   figure; imagesc((ones(size(BKmsk))-BKmsk).*coefMapPDI); axis image;colorbar;
hold on;
hAxes2 = axes;
coefMap=max(abs(coefMapPDI),[],3);
%              coefMap=max(abs(coefMapPDI(3:152,:,:)),[],3);
h2=imagesc(coefMap);
set(h2,'AlphaData',abs(coefMap)/(max(abs(coefMap(:)))/4).*((abs(coefMap))>0.35).*handles.BKmsk.*handles.BWBrain) % only show coefMap larger than 0.35
% colormap(hAxes2,VzCmap);
% caxis(hAxes2,[-1 1])
colormap(hAxes2,hot);
caxis(hAxes2,[0.3 1])
colorbar
hold off
axis equal tight
axis(hAxes2,'off')
linkaxes([hAxes1,hAxes2])
title('PDI-based activation map')

%%%%
        m1 = squeeze(squeeze(PDIrm(109,198,3,:)));
        m2 = squeeze(squeeze(PDIrm(12, 13, 3, :))); 
        m3 = squeeze(squeeze(PDIrm(17, 228, 3, :)));
        m4 = squeeze(squeeze(PDIrm(79, 4, 3, :)));
        m5 = squeeze(squeeze(PDIrm(16, 125, 3, :))); %%%% 0822 data

        corrmap1 = CorrMap(squeeze(PDIrm(:,:,3,:)), m1');
        corrmap2 = CorrMap(squeeze(PDIrm(:,:,3,:)), m2');
        corrmap3 = CorrMap(squeeze(PDIrm(:,:,3,:)), m3');
        corrmap4 = CorrMap(squeeze(PDIrm(:,:,3,:)), m5');
        corrmap5 = CorrMap(squeeze(PDIrm(:,:,3,:)), m5');
                figure; subplot(511);plot(m1);subplot(512);plot(m2);subplot(513);plot(m3);subplot(514);plot(m4);subplot(515);plot(m5)
        figure; subplot(231);imagesc(corrmap1);axis image;subplot(232);imagesc(corrmap2);axis image;
        subplot(233);imagesc(corrmap3);axis image; subplot(234);imagesc(corrmap4);axis image;
        subplot(235);imagesc(corrmap5);axis image;
%         handles.BKmsk =
%         (corrmap1<0.3).*(corrmap2<0.3).*(corrmap3<0.3).*(corrmap4<0.3).*(corrmap5<0.3);%
%         for 0520 data
        handles.BKmsk = (corrmap1<0.7).*(corrmap2<0.7).*(corrmap3<0.7).*(corrmap4<0.8).*(corrmap5<0.7);       
figure; imagesc(handles.BKmsk); axis image;

figure;
hAxes1 = axes;
imagesc(abs(ImgBKGDPDI));
caxis(hAxes1,[0.3 1]);
colormap(hAxes1,gray)
colorbar('Ticks',[-10 10], 'TickLabels',{[],[]} );
axis equal tight
hold on;
hAxes2 = axes;
coefMap=max(abs(coefMapPDI),[],3);
h2=imagesc(coefMap);
set(h2,'AlphaData',abs(coefMap)*0.9./min(coefMap(:)).*((abs(coefMap))>0.1).*handles.BKmsk.*handles.BWBrain) % only show coefMap larger than 0.35
colormap(hAxes2,hot);
caxis(hAxes2,[0.1 0.6])
colorbar
hold off
axis equal tight
axis(hAxes2,'off')
linkaxes([hAxes1,hAxes2])
title('PDI-based activation map')

guidata(hObject, handles);
% --- Executes on button press in BTN_SeedConnectivity.
function BTN_SeedConnectivity_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_SeedConnectivity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;
addpath('D:\CODE\DataProcess\BingxueLiu\Code fUS DATA analysis')
addpath('D:\CODE\DataProcess\BingxueLiu\Code fUS DATA analysis\SubFunctions')
[VzCmap, VzCmapDn,VzCmapUp]=Colormaps_fUS;
rFrame=1/handles.tIntV;
[nz,nx,nt]=size(handles.cData);
cData=handles.cData;
FCMethod = questdlg('FC analysis method?', ...
                         'Select', ...
                         'Seed-based', 'SVD-based', 'Cancel');
switch FCMethod
    case 'Seed-based'
        xL_ROI=6;
        zL_ROI=6;
        ButtonName = questdlg('Select New ROI?', ...
            'Select', ...
            'Yes', 'No', 'Cancel');
        switch ButtonName
            case 'Yes'
                % selecte ROI
                %         [xROI, zROI]=ginput(6);
                %         xROI=floor(xROI);
                %         zROI=floor(zROI);
                %         handles.BW=roipoly(handles.data,xROI,zROI);
                
                [xROI, zROI]=ginput(1);
                handles.BW=zeros(nz,nx);
                xROI=floor(xROI);
                zROI=floor(zROI);
                handles.BW(zROI+[-zL_ROI:zL_ROI],xROI+[-xL_ROI:xL_ROI])=1;
        end
        %% band pass filtering
        [B,A]=butter(3,[0.02,0.2]/(rFrame/2),'bandpass');
        % [B,A]=butter(3,0.03,'high');
        sData1(:,:,101:100+nt)=cData;
        for iCC=1:100
            sData1(:,:,iCC)=sData1(:,:,201-iCC);
        end
        sData2=filter(B,A,sData1,[],3);    %
        handles.dataHP=sData2(:,:,101:end); % High pass filtered signal
        %% low pass filtering
        % cData=cData-mean(cData,3);
        % [B,A]=butter(3,0.2/(rFrame/2),'low');
        % sData1(:,:,101:100+nt)=cData;
        % for iCC=1:100
        %     sData1(:,:,iCC)=sData1(:,:,201-iCC);
        % end
        % sData2=filter(B,A,sData1,[],3);    %
        % handles.dataHP=sData2(:,:,101:end); % High pass filtered signal
        %% seed ROI signal after HP
        DataFltd=handles.dataHP;
        [nz,nx,nt]=size(DataFltd);
        %% cross correlation coefficient
        dataROI=mean(mean(DataFltd.*handles.BW,1),2);
        dataROI=(dataROI-mean(dataROI));
        dataROI=dataROI./sqrt(mean(dataROI.^2));
        
        cDataN=DataFltd-mean(DataFltd,3);
        cDataN=cDataN./sqrt(mean(cDataN.^2,3));
        
        coefMap=mean(cDataN.*dataROI,3);
        coefMap(find(isnan(coefMap)==1))=0;
        
        z=sqrt(nt-3)*log((1+coefMap)./(1-coefMap))/2;
        zMsk=zeros(size(z));
        zMsk(z>5)=1;
        actMap=coefMap;%.*zMsk;
        %% activation map plot
        figure;
        if handles.DnUpSlt==0 % all flow direction
            ImgShow=handles.ImgShow(:,:,3);
            hAxes1 = axes;
            imagesc(ImgShow);
            caxis(hAxes1,[min(min(ImgShow(:))*1.2, max(ImgShow(:))*1.1) max(min(ImgShow(:))*1.2, max(ImgShow(:))*1.1)]);
            colormap(hAxes1,gray)
            colorbar('Ticks',[-10 10], 'TickLabels',{[],[]} );
            axis equal tight
        else
            hAxes1 = axes;
            ImgShow=abs(handles.ImgShow(:,:,handles.DnUpSlt));
            imagesc(ImgShow);
            %     caxis(hAxes1,[min(min(ImgShow(:))*1.2, max(ImgShow(:))*1.1) max(min(ImgShow(:))*1.2, max(ImgShow(:))*1.1)]);
            colormap(hAxes1,gray)
            colorbar('Ticks',[-10 10], 'TickLabels',{[],[]} );
            axis equal tight
        end
        
        hold on;
        hAxes2 = axes;
        h2=imagesc(actMap);
        % set(h2,'AlphaData',(abs(actMap)).^1)
        AlphaMsk=abs(actMap);
        AlphaMsk(actMap>0.3)=0.7;
        AlphaMsk(actMap>0.4)=0.8;
        AlphaMsk(actMap>0.5)=1;
        set(h2,'AlphaData',AlphaMsk.^1.2)
        % colormap(hAxes2,VzCmap);
        % caxis(hAxes2,[-0.5 0.5])
        colormap(hot);
        caxis(hAxes2,[0.1 1])
        colorbar
        hold off
        axis equal tight
        axis(hAxes2,'off')
        linkaxes([hAxes1,hAxes2])
        %% ROI plot
        if handles.DnUpSlt==0 % all flow direction
            if handles.DataSlt==0 % PDI, all frequency
                ImgShow=handles.ImgShow(:,:,3);
                axes(handles.axes1)
                imagesc(ImgShow);
                colormap(gray);
                colorbar
                caxis([min(min(ImgShow(:))*1.1, max(ImgShow(:))*0.96) max(min(ImgShow(:))*1.1, max(ImgShow(:))*0.96)]);
                handles.cData=squeeze(handles.data(:,:,3,:));
            else
                axes(handles.axes1)
                imagesc(handles.ImgShow(:,:,3));
                colormap(VzCmap);
                caxis(handles.cAxis);
                colorbar
                handles.cData=squeeze(max(abs(handles.data(:,:,1:2,:)),[],3));
            end
        else
            if handles.DataSlt==0 % PDI, all frequency
                ImgShow=handles.ImgShow(:,:,handles.DnUpSlt);
                axes(handles.axes1)
                imagesc(ImgShow);
                colormap(gray);
                colorbar
                caxis([min(min(ImgShow(:))*1.04, max(ImgShow(:))*0.96) max(min(ImgShow(:))*1.04, max(ImgShow(:))*0.96)]);
                handles.cData=squeeze(handles.data(:,:,handles.DnUpSlt,:));
                title(['Up(1)/Down(2): ', num2str(handles.DnUpSlt)])
            else
                [VzCmap, VzCmapDn,VzCmapUp]=Colormaps_fUS;
                axes(handles.axes1)
                imagesc(handles.ImgShow(:,:,handles.DnUpSlt));
                caxis(handles.cAxis);
                colormap(VzCmap)
                colorbar;
                handles.cData=squeeze(handles.data(:,:,handles.DnUpSlt,:));
            end
        end
        axes(handles.axes1)
        hold on
        h3=imagesc(handles.BW); % down flow
        set(h3,'AlphaData',0.5*double(abs(handles.BW)>0.2))
        % colormap(hot);
        % caxis(handles.cAxis);
        colorbar
        hold off
        %% for test - view power spectrum
        % xSlt=61
        % zSlt=10
        % tCoor=[1:nt]*handles.tIntV;
        % fCoor=linspace(-rFrame/2,rFrame/2,nt);
        % oriData=squeeze(handles.cData(zSlt,xSlt,:));
        % fltData=squeeze(handles.dataHP(zSlt,xSlt,:));
        % % oriData=squeeze(mean(mean(handles.cData(:,:,:),1),2));
        % % fltData=squeeze(mean(mean(handles.dataHP(:,:,:),1),2));
        % % eqNoise=squeeze(mean(mean(handles.eqNoise(:,:,:),1),2));
        % oriPSD=fftshift(fft(oriData,nt));
        % fltPSD=fftshift(fft(fltData,nt));
        %
        % xFltData = fltData-mean(fltData);
        % xFltData = xFltData./sqrt(mean(xFltData.^2));
        %
        % Fig=figure;
        % set(Fig,'Position',[400 400 800 800])
        % subplot(4,1,1);hold on; plot(tCoor,oriData);title('Raw data'); xlabel('t [s]'); grid on
        % subplot(4,1,2);hold on; plot(tCoor,fltData);title('BP filtered data'); xlabel('t [s]'); grid on
        % subplot(4,1,3);hold on; plot(tCoor,xFltData);title('xRoss calcu filtered data'); xlabel('t [s]'); grid on
        % subplot(4,1,4); hold on,plot(fCoor,abs(fltPSD));  %hold on; plot(fCoor,abs(oriPSD));
        % ylim([0 1.2*max(abs(fltPSD(:)))]); xlabel('f [HZ]'); grid on
        %
        % % figure;
        % % subplot(3,1,1);hold on; plot(tCoor,oriData);title('Raw data'); xlabel('t [s]')
        % % subplot(3,1,2);hold on; plot(tCoor,fltData);title('BP filtered data'); xlabel('t [s]')
        % % subplot(3,1,3);hold on; plot(fCoor,abs(oriFPS)); hold on,plot(fCoor,abs(fltFPS)); ylim([0 1.2*max(abs(fltFPS(:)))]); xlabel('f [HZ]')
        %% %%%%%%%%%%%%%%%%%%%%% for literature data test %%%%%%%%%%%%%%%%%%%%%%%%
        % fech = 1/1.5;
        % % [B,A] = butter(3,0.05/fech*2,'high');
        % [B,A] = butter(3,[0.05 0.2]/fech*2,'bandpass');
        % nt=144;sData1(:,:,51:50+nt)=cData;
        % for iCC=1:50
        %     sData1(:,:,iCC)=sData1(:,:,101-iCC);
        % end
        % sData2=filter(B,A,sData1,[],3);    %
        % DataFltd=sData2(:,:,51:end); % High pass filtered signal
        % %% JT FC processing code
        % [nz,nx,nt]=size(DataFltd);
        % dataROI=mean(mean(DataFltd.*handles.BW,1),2);
        % dataROI=(dataROI-mean(dataROI));
        % dataROI=dataROI./sqrt(mean(dataROI.^2));
        %
        % cDataN=DataFltd-mean(DataFltd,3);
        % cDataN=cDataN./sqrt(mean(cDataN.^2,3));
        %
        % coefMap=mean(cDataN.*dataROI,3);
        % coefMap(find(isnan(coefMap)==1))=0;
        % z=sqrt(nt-3)*log((1+coefMap)./(1-coefMap))/2;
        % zMsk=zeros(size(z));
        % zMsk(z>3)=1;
        % actMap=coefMap.*zMsk;
        %
        % %% literature FC processing code
        % % DopplerM = mean(DataFltd, 3);
        % % DopplerN = bsxfun(@minus,DataFltd, DopplerM);
        % % DopplerN = bsxfun(@rdivide,DopplerN, sqrt(mean(DopplerN.^2,3)));
        % %
        % % dataROI=mean(mean(DataFltd.*handles.BW,1),2);
        % % dataROI=(dataROI-mean(dataROI));
        % % dataROI=dataROI./sqrt(mean(dataROI.^2));
        % %
        % % coefMap=mean(bsxfun(@times,DopplerN,dataROI),3);
        % % z=sqrt(nt-3)*log((1+coefMap)./(1-coefMap))/2;
        % % zMsk=zeros(size(z));
        % % zMsk(z>2.9)=1;
        % % actMap=coefMap.*zMsk;
        % %% FC plot
        % figure;
        % if handles.DnUpSlt==0 % all flow direction
        %     ImgShow=handles.ImgShow(:,:,3);
        %     hAxes1 = axes;
        %     imagesc(ImgShow);
        % %     caxis(hAxes1,[min(min(ImgShow(:))*1.2, max(ImgShow(:))*1.1) max(min(ImgShow(:))*1.2, max(ImgShow(:))*1.1)]);
        %     caxis(hAxes1,[-1.8 0.5])
        %     colormap(hAxes1,gray)
        %     colorbar('Ticks',[-10 10], 'TickLabels',{[],[]} );
        % %     axis equal tight
        % else
        %     hAxes1 = axes;
        %     ImgShow=(handles.ImgShow(:,:,handles.DnUpSlt));
        %     imagesc(ImgShow);
        % %     caxis(hAxes1,[min(min(ImgShow(:))*1.2, max(ImgShow(:))*1.1) max(min(ImgShow(:))*1.2, max(ImgShow(:))*1.1)]);
        %     colormap(hAxes1,gray)
        %     colorbar('Ticks',[-10 10], 'TickLabels',{[],[]} );
        % %     axis equal tight
        % end
        %
        % hold on;
        % hAxes2 = axes;
        % h2=imagesc(actMap);
        % set(h2,'AlphaData',(abs(actMap)).^1.1)
        % % colormap(hAxes2,VzCmap);
        % % caxis(hAxes2,[-0.5 0.5])
        % colormap(hot);
        % caxis(hAxes2,[0 0.8])
        % colorbar
        % hold off
        % % axis equal tight
        % axis(hAxes2,'off')
        % linkaxes([hAxes1,hAxes2])
        % %% ROI plot
        % if handles.DnUpSlt==0 % all flow direction
        %     if handles.DataSlt==0 % PDI, all frequency
        %         ImgShow=handles.ImgShow(:,:,3);
        %         axes(handles.axes1)
        %         imagesc(ImgShow);
        %         colormap(gray);
        %         colorbar
        %         caxis([-1.8 0.5]);
        %         handles.cData=squeeze(handles.data(:,:,3,:));
        %     else
        %         axes(handles.axes1)
        %         imagesc(handles.ImgShow(:,:,3));
        %         colormap(VzCmap);
        %         caxis(handles.cAxis);
        %         colorbar
        %         handles.cData=squeeze(max(abs(handles.data(:,:,1:2,:)),[],3));
        %     end
        % else
        %     if handles.DataSlt==0 % PDI, all frequency
        %         ImgShow=handles.ImgShow(:,:,handles.DnUpSlt);
        %         axes(handles.axes1)
        %         imagesc(ImgShow);
        %         colormap(gray);
        %         colorbar
        %         caxis([min(min(ImgShow(:))*1.04, max(ImgShow(:))*0.96) max(min(ImgShow(:))*1.04, max(ImgShow(:))*0.96)]);
        %         handles.cData=squeeze(handles.data(:,:,handles.DnUpSlt,:));
        %         title(['Up(1)/Down(2): ', num2str(handles.DnUpSlt)])
        %     else
        %         [VzCmap, VzCmapDn,VzCmapUp]=Colormaps_fUS;
        %         axes(handles.axes1)
        %         imagesc(handles.ImgShow(:,:,handles.DnUpSlt));
        %         caxis(handles.cAxis);
        %         colormap(VzCmap)
        %         colorbar;
        %         handles.cData=squeeze(handles.data(:,:,handles.DnUpSlt,:));
        %     end
        % end
        % axes(handles.axes1)
        % hold on
        % h3=imagesc(handles.BW); % down flow
        % set(h3,'AlphaData',0.3*double(abs(handles.BW)>0.2))
        % % colormap(hot);
        % % caxis(handles.cAxis);
        % colorbar
        % hold off
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'SVD-based'
        [nz,nx,nSVD]=size(cData);
        rData=reshape(cData,[nz*nx,nSVD]);
        Sn=rData./mean(rData,2); % temporal normalization
        Sm=Sn-mean(Sn,2); % temporal 0 centering
        S=Sm;
        % S=rData./max(rData,[],2);
        S_COVt=(S'*S);
        [V,D]=eig(S_COVt); % V is the right singular Vector of S/eigenvector; D is the eigenvalue/square of Singular value
        for it=1:nSVD
            Ddiag(it)=abs(sqrt(D(it,it)));
        end
        Ddiag=20*log10(Ddiag/max(Ddiag)); % singular value in db
        [Ddesc, Idesc]=sort(Ddiag,'descend');
        % figure,plot(Ddesc);
        for it=1:nSVD
            Vdesc(:,it)=V(:,Idesc(it));
        end
        UDelta=S*Vdesc;
        U=reshape(UDelta,[nz,nx,nt]);
        %%
        [VzCmap, VzCmapDn,VzCmapUp]=Colormaps_fUS;
        Fig=figure;
        set(Fig,'Position',[400 400 1600 600])
        for r=1:6
            subplot(2,3,r);imagesc(U(:,:,r));colormap(VzCmap);
            caxis([-1 1]); colorbar; axis equal tight;
            title(['r=',num2str(r)])
        end
end
guidata(hObject, handles);


% --- Executes on button press in BTN_RESET.
function BTN_RESET_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_RESET (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.DefPath='G:\PROJ-R-Stroke';
addpath('D:\CODE\DataProcess\BingxueLiu\Code fUS DATA analysis')
addpath('D:\CODE\DataProcess\BingxueLiu\Code fUS DATA analysis\SubFunctions')
handles.DnUpSlt=1;
handles.DataSlt=1;
handles.nROI=1;
guidata(hObject, handles);
% --- Executes on button press in BTN_SAVE.
function BTN_SAVE_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_SAVE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% clc;
if handles.DataSlt==0 % PDI
    SaveName='PDI-tCourse-Corr-';
    iSaveName=[SaveName,'1'];
    iFile=1;
    while exist([handles.DefPath,iSaveName,'.mat'])==2
        iFile=iFile+1;
        iSaveName=[SaveName,num2str(iFile)];
    end
end
% time course
tc.mData=handles.mData;
tc.stdData=handles.stdData;
tc.semData=handles.semData;
tc.BW=handles.BW;
P=handles.P;
% pearson correlation
pc.cData = handles.cData;
pc.cData0 = handles.cData0;
pc.coefMapPDI = handles.coefMapPDI;
pc.HRF = handles.HRF;
pc.mResponse = handles.mResponse;
pc.BKmsk = handles.BKmsk;

Samples = handles.Samples;
Stim = handles.Stim;

save ([handles.DefPath,iSaveName,'.mat'],'tc','pc','P','Samples','Stim');

% saveas(handles.tCourseFhandle,[handles.DefPath,iSaveName,'.bmp'],'bmp');
% saveas(handles.tCourseFhandle,[handles.DefPath,iSaveName,'.fig'],'fig');
guidata(hObject, handles);


% --- Executes on button press in BTN_relativeCHG2Base.
function BTN_relativeCHG2Base_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_relativeCHG2Base (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=handles.data;

%%%%% handles.PDI0: original PDI 
%     handles.PDI : PDI signal after motion correction
%     handles.
