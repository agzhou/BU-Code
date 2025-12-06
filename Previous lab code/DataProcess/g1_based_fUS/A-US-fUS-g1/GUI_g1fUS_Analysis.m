function varargout = GUI_g1fUS_Analysis(varargin)
% GUI_G1FUS_ANALYSIS MATLAB code for GUI_g1fUS_Analysis.fig
%      GUI_G1FUS_ANALYSIS, by itself, creates a new GUI_G1FUS_ANALYSIS or raises the existing
%      singleton*.
%
%      H = GUI_G1FUS_ANALYSIS returns the handle to a new GUI_G1FUS_ANALYSIS or the handle to
%      the existing singleton*.
%
%      GUI_G1FUS_ANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_G1FUS_ANALYSIS.M with the given input arguments.
%
%      GUI_G1FUS_ANALYSIS('Property','Value',...) creates a new GUI_G1FUS_ANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_g1fUS_Analysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_g1fUS_Analysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_g1fUS_Analysis

% Last Modified by GUIDE v2.5 07-May-2021 15:28:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_g1fUS_Analysis_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_g1fUS_Analysis_OutputFcn, ...
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


% --- Executes just before GUI_g1fUS_Analysis is made visible.
function GUI_g1fUS_Analysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_g1fUS_Analysis (see VARARGIN)

% Choose default command line output for GUI_g1fUS_Analysis
handles.output = hObject;
handles.CodePath=pwd;
addpath(handles.CodePath);
addpath([handles.CodePath, '\SubFunctions'])
handles.DefPath='G:\PROJ-D-vUS';
handles.DefProbe='L22-14v';
handles.DefCtnFreq=18.5; % default center frequency, MHz
handles.DefPitch=0.1;    % default transducer element pitch, mm
handles.DefSpeed=1540;   % default sound speed, m/s or mm/ms
handles.DefCutFreq=25;   % default cutoff frequency for Power Doppler processing
handles.SVDRank=[25 1000];   % SVD rank
handles.newSlect=1;
handles.DataType=1; 
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_g1fUS_Analysis wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_g1fUS_Analysis_OutputFcn(hObject, eventdata, handles) 
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
addpath(handles.CodePath);
addpath([handles.CodePath, '\SubFunctions'])
clc
[filename,datapath]=uigetfile(handles.DefPath);
handles.DefPath=datapath;
guidata(hObject, handles);
disp('Loading data...');
load ([datapath, filename]);
disp('Data was loaded...');
handles.IQ=IQ; % 
handles.P=P;
handles.DAQinfo.rFrame=P.CCFR; % Hz
handles.DAQinfo.f0=P.frequency*1e3; % Hz
handles.DAQinfo.C=P.vSound; % m/s

[handles.nz, handles.nx, handles.nt]=size(handles.IQ);
fileinfo=strsplit(filename(1:end-4),'-');
handles.Agl=str2num(fileinfo{2});
handles.nAgl=str2num(fileinfo{3});
handles.fCC=str2num(fileinfo{4});
handles.nCC=str2num(fileinfo{5});
handles.iPlane=str2num(fileinfo{6});

prompt={'Transducer Center frequency (MHz): ',...
    'Transducer element pitch (mm): ',...
    'Sound speed (m/s or mm/ms): ',...
    'Angle:',...
    'nAngle:',...
    'CC frame rate:',...
    'nCC frames:',...
    'Processing plane:'};
name='File info';
defaultvalue={num2str(handles.DefCtnFreq),...
    num2str(handles.DefPitch),...
    num2str(handles.DefSpeed),...
    num2str(handles.Agl),...
    num2str(handles.nAgl),...
    num2str(handles.fCC),...
    num2str(handles.nCC),...
    num2str(handles.iPlane)};
numinput=inputdlg(prompt,name, 1, defaultvalue);
handles.DefCtnFreq=str2num(numinput{1});
handles.DefPitch=str2num(numinput{2});
handles.DefSpeed=str2num(numinput{3});
handles.Agl=str2num(numinput{4});
handles.nAgl=str2num(numinput{5});
handles.fCC=str2num(numinput{6}); % CC frame rate, Hz
handles.nCC=str2num(numinput{7});
handles.SVDRank(2)=handles.nCC;
handles.iPlane=str2num(numinput{8});
handles.P=P;
handles.filenameBase=filename(1:end-6);
axes(handles.axes1)
imagesc(abs((handles.IQ(:,:,100)))); 
colormap(handles.axes1, hot)

% figure,
% imagesc(abs((handles.dataRAW(:,:,10)))); 
% colormap(hot)
set(handles.Disp_FileName,'string',filename);
guidata(hObject, handles);

% --- Executes on button press in BTN_IQ2sIQ.
function BTN_IQ2sIQ_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_IQ2sIQ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
addpath(handles.CodePath);
addpath([handles.CodePath, '\SubFunctions'])
clc;

prompt={'SVD Rank (low):', ['SVD Rank (High):(Max Rank: ',num2str(handles.nCC),')'],'High pass cutoff frequency (Hz)',...
    ['nCC_process (nCC total: ',num2str(handles.nCC),')']};
name='Power Doppler data processing';
defaultvalue={num2str(handles.SVDRank(1)), num2str(handles.SVDRank(2)), num2str(handles.DefCutFreq),num2str(handles.nCC)};
numinput=inputdlg(prompt,name, 1, defaultvalue);
handles.SVDRank=[str2num(numinput{1}),str2num(numinput{2})];
handles.DefCutFreq=str2num(numinput{3});
nCC_proc=str2num(numinput{4});
PRSSinfo.SVDrank=handles.SVDRank;
PRSSinfo.HPfC=handles.DefCutFreq;
PRSSinfo.NEQ=0;
PRSSinfo.rFrame=handles.fCC;
[handles.sIQ, handles.sIQHP, handles.sIQHHP, handles.eqNoise]=IQ2sIQ(handles.IQ,PRSSinfo);
axes(handles.axes1)
imagesc(abs((handles.sIQ(:,:,100)))); axis image;
colormap(handles.axes1, hot)
guidata(hObject, handles);


% --- Executes on button press in BTN_PDIfUS.
function BTN_PDIfUS_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_PDIfUS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
addpath(handles.CodePath);
addpath([handles.CodePath, '\SubFunctions'])
clc;
[handles.PDI]=sIQ2PDI(handles.sIQHP); % SVD + HPfilter = PDI

%[handles.PDI]=sIQ2PDI(handles.sIQ);
% [handles.PDIHP]=sIQ2PDI_GPU(handles.sIQHP);
% [handles.PDIHHP]=sIQ2PDI_GPU(handles.sIQHHP);
axes(handles.axes1)
imagesc(log(abs((handles.PDI(:,:,3))))); colorbar; axis image;
colormap(handles.axes1, hot)
guidata(hObject, handles);

% --- Executes on button press in BTN_sIQ2GG.
function BTN_sIQ2GG_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_sIQ2GG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
addpath(handles.CodePath);
addpath([handles.CodePath, '\SubFunctions'])
clc;
prompt={'g1 start', 'g1 nt','g1 nTau'};
name='Power Doppler data processing';
defaultvalue={num2str(1),num2str(handles.nCC), '100'};
numinput=inputdlg(prompt,name, 1, defaultvalue);
handles.PRSinfo.g1StartT=str2num(numinput{1});
handles.PRSinfo.g1nT=str2num(numinput{2});
handles.PRSinfo.g1nTau=str2num(numinput{3});
handles.GG= sIQ2GG_GPU(handles.sIQ, handles.PRSinfo);
% handles.GG= sIQ2GG(handles.sIQ, handles.PRSinfo);
%% SHOW GG(:,:,iTau)
iTau=str2num(get(handles.Input_iTau,'string'));
axes(handles.axes1)
imagesc((abs((handles.GG(:,:,iTau))))); colorbar; axis image;
colormap(handles.axes1, jet)
title(['iTau=',num2str(iTau)])
caxis([0 1])
guidata(hObject, handles);

% --- Executes on button press in BTN_g1fUS.
function BTN_g1fUS_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_g1fUS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Ag1=abs(handles.GG(:,:,1))-min(abs(handles.GG(:,:,1:10)),[],3);
handles.Ag12=abs(handles.GG(:,:,1))-min(abs(handles.GG(:,:,:)),[],3);
% handles.InTg1=sum((handles.GG(:,:,1:10))-min(handles.GG(:,:,1:10),[],3),3);
handles.InTg1=abs(sum((handles.GG(:,:,1:10)),3));
% handles.imInTg1=sum(abs(imag(handles.GG(:,:,1:10))),3);
handles.imInTg1=abs(sum(abs((imag(handles.GG(:,:,1:10)))),3));
handles.reInTg1=abs(sum(abs((real(handles.GG(:,:,1:10)))),3));

PDI=handles.PDI(:,:,3)./handles.eqNoise.^1.8;

figure;
imagesc(handles.imInTg1(:,20:end-20));axis image; colorbar;colormap(hot);axis off;

figure;
imagesc(log(PDI(:,20:end-20)));axis image; colorbar;colormap(hot);axis off;

Fig=figure;
set(Fig,'Position',[100 100 1600 800]);
subplot(2,3,1)
imagesc(handles.Ag1);
colormap(hot)
title('g1(1)-min(g1(1:10)')
subplot(2,3,3)
imagesc((abs(handles.InTg1)));
colormap(hot)
title('sum(abs(g1(1:10)))')


subplot(2,3,2)
imagesc((handles.reInTg1))
colormap(hot)
title('sum(real(g1(1:10)))')
subplot(2,3,5)
imagesc((abs(handles.imInTg1)));
colormap(hot)
title('sum(imag(g1(1:10)))')


subplot(2,3,4)
imagesc(handles.Ag12);
colormap(hot)
title('g1(1)-min(g1(:)')
subplot(2,3,6)
imagesc(log(abs(PDI)));
colormap(hot)
title('PDI)')
guidata(hObject, handles);

function BTN_PixSelectPlot_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_PixSelectPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
addpath(handles.CodePath);
addpath([handles.CodePath, '\SubFunctions'])
clc;
prompt={'New Select(0: No; 1: Yes)', 'DataType(0:Phantom; 1: invivo)'};
name='Select';
defaultvalue={num2str(handles.newSlect),num2str(handles.DataType)};
numinput=inputdlg(prompt,name, 1, defaultvalue);
handles.newSlect=str2num(numinput{1});
handles.DataType=str2num(numinput{2});
if handles.newSlect==1
    axes(handles.axes1);
    [handles.slt(2), handles.slt(1)]=ginput(1); % [x z]
    handles.slt=round(handles.slt);
end
guidata(hObject, handles);

% handles.slt = [40,60];

%% Frequency analysis
IQ=(handles.IQ(handles.slt(1),handles.slt(2),:));
sIQ=(handles.sIQ(handles.slt(1),handles.slt(2),:));
sIQHP=(handles.sIQHP(handles.slt(1),handles.slt(2),:));

IQAVG=squeeze((IQ));
nt=length(IQAVG);
fCoor=linspace(-handles.fCC/2,handles.fCC/2,nt)';
fIQ=fftshift(fft(IQAVG));

sIQAVG=squeeze(((sIQ)));
nt=length(sIQAVG);
FsIQ=fftshift(fft(sIQAVG));

sIQHPAVG=squeeze((sIQHP));
PDISVDHP=mean(abs(sIQHPAVG).^2);
nt=length(sIQHPAVG);
FsIQHP=fftshift(fft(sIQHPAVG));

Pneg=sum(abs(FsIQ(1:floor(end/2)))-median(abs(FsIQ)));
Ppos=sum(abs(FsIQ(floor(end/2)+1:end))-median(abs(FsIQ)));
rNeg=Pneg/(Pneg+Ppos);
rPos=Ppos/(Pneg+Ppos);
tCoor=linspace(1,nt,nt)/handles.fCC*1000;

fig=figure;
set(fig,'Position',[200 400 800 650])
subplot(2,1,1),plot(fCoor,abs(fIQ),'k')
hold on; plot(fCoor,abs(FsIQ),'r')
hold on, plot(fCoor,abs(FsIQHP),'b')
%     xlim([-800 800])
ylim([0 max(abs(FsIQ))])
xlabel('Frequency [Hz]')
ylabel('Power Spectrum Density')
% legend({'Raw','SVD','HP'})
title(['Pneg=',num2str(rNeg),'; Ppos=',num2str(rPos)])

subplot(2,1,2),plot(tCoor,(angle(squeeze(IQ))),'k');
hold on, plot(tCoor,(angle(squeeze(sIQ))),'r');
hold on,plot(tCoor,(angle(squeeze(sIQHP))),'b');
title('phase(IQ)');% ylim([0,1])
xlabel('time lag, [ms]')
legend({'Raw','SVD','SVDHP'});

%% sIQ analysis
fig=figure;
set(fig,'Position',[200 400 400 650])
plot(tCoor,sIQHP); 
title('sIQ'); 
xlabel('time lag, [ms]')


%% GG analysis
handles.tauCoor=linspace(1,handles.PRSinfo.g1nTau,handles.PRSinfo.g1nTau)/handles.fCC; % s
gg=squeeze(handles.GG(handles.slt(1),handles.slt(2),:));
Fig=figure;
set(Fig,'Position',[300 200 400 1000]);
subplot(5,1,2);
plot(handles.tauCoor*1e3,squeeze(real(gg)),'k')
xlabel('tau [ms]')
grid on
title(['GG-real,', '[', num2str(handles.slt),']'])

subplot(5,1,3);
plot(handles.tauCoor*1e3,squeeze(imag(gg)),'b')
xlabel('tau [ms]')
grid on
title(['GG-imag'])

subplot(5,1,5);
plot(squeeze((gg)),'k')
xlim([-1 1]);ylim([-1 1])
grid on
title(['GG-complex'])

subplot(5,1,1);
plot(handles.tauCoor*1e3,squeeze(abs(gg)),'r')
xlabel('tau [ms]')
ylim([0 1])
grid on
title(['GG-magnitude'])

subplot(5,1,4);
plot(handles.tauCoor*1e3,squeeze(angle(gg)),'k')
xlabel('tau [ms]'); ylim([-pi, pi]);
grid on
title(['GG-phase'])

%% plot multiple pixels gg
prompt={'Multiple Pixes? (0: No; 1: Yes)'};
name='Select';
numinput=inputdlg(prompt,name, 1, {'1'});
npix = str2num(numinput{1});

if handles.newSlect==1
    axes(handles.axes1);
    mslt = ginput(npix); % [x z]
    mslt=round(mslt);
end
mslt = fliplr(mslt);
Fig=figure;
set(Fig,'Position',[300 400 1200 600]);
for ipix = 1: npix
ggs(:,ipix)=squeeze(handles.GG(mslt(ipix,1),mslt(ipix,2),:));
subplot(2,2,1);
plot(handles.tauCoor*1e3,squeeze(real(ggs(:,ipix))))
xlabel('tau [ms]')
grid on
title(['GG-real'])
hold on

subplot(2,2,2);
plot(handles.tauCoor*1e3,squeeze(imag(ggs(:,ipix))))
xlabel('tau [ms]')
grid on
title(['GG-imag'])
hold on

subplot(2,2,3);
plot(squeeze((ggs(:,ipix))))
xlim([-1 1]);ylim([-1 1])
grid on
title(['GG-complex'])
hold on
legend(num2str(mslt))

subplot(2,2,4);
plot(handles.tauCoor*1e3,squeeze(abs(ggs(:,ipix))))
xlabel('tau [ms]')
ylim([0 1])
grid on
title(['GG-magnitude'])
hold on
end

    



% --- Executes on button press in BTN_PixSelectPlot.
% function BTN_PixSelectPlot_Callback(hObject, eventdata, handles)
% % hObject    handle to BTN_PixSelectPlot (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% addpath(handles.CodePath);
% addpath([handles.CodePath, '\SubFunctions'])
% clc;
% prompt={'New Select(0: No; 1: Yes)', 'DataType(0:Phantom; 1: invivo)'};
% name='Select';
% defaultvalue={num2str(handles.newSlect),num2str(handles.DataType)};
% numinput=inputdlg(prompt,name, 1, defaultvalue);
% handles.newSlect=str2num(numinput{1});
% handles.DataType=str2num(numinput{2});
% if handles.newSlect==1
%     axes(handles.axes1);
%     [handles.slt(2), handles.slt(1)]=ginput(1); % [x z]
%     handles.slt=round(handles.slt);
% end
% guidata(hObject, handles);
% 
%     %% Frequency analysis
%     IQ=(handles.IQ(handles.slt(1),handles.slt(2),:));
%     sIQ=(handles.sIQ(handles.slt(1),handles.slt(2),:));
%     sIQHP=(handles.sIQHP(handles.slt(1),handles.slt(2),:));
%     
%     IQAVG=squeeze((IQ));
%     nt=length(IQAVG);
%     fCoor=linspace(-handles.fCC/2,handles.fCC/2,nt)';
%     fIQ=fftshift(fft(IQAVG));
%     
%     sIQAVG=squeeze(((sIQ)));
%     nt=length(sIQAVG);
%     FsIQ=fftshift(fft(sIQAVG));
%     
%     sIQHPAVG=squeeze((sIQHP));
%     PDISVDHP=mean(abs(sIQHPAVG).^2);
%     nt=length(sIQHPAVG);
%     FsIQHP=fftshift(fft(sIQHPAVG));
%     
%     Pneg=sum(abs(FsIQ(1:floor(end/2)))-median(abs(FsIQ)));
%     Ppos=sum(abs(FsIQ(floor(end/2)+1:end))-median(abs(FsIQ)));
%     rNeg=Pneg/(Pneg+Ppos);
%     rPos=Ppos/(Pneg+Ppos);
%     tCoor=linspace(1,nt,nt)/handles.fCC*1000;
%     
%     fig=figure;
%     set(fig,'Position',[1200 400 800 650])
%     subplot(2,1,1),plot(fCoor,abs(fIQ),'k')
%     hold on; plot(fCoor,abs(FsIQ),'r')
%     hold on, plot(fCoor,abs(FsIQHP),'b')
% %     xlim([-800 800])
%     ylim([0 max(abs(FsIQ))])
%     xlabel('Frequency [Hz]')
%     ylabel('Power Spectrum Density')
%     % legend({'Raw','SVD','HP'})
%     title(['Pneg=',num2str(rNeg),'; Ppos=',num2str(rPos)])
%     
%     subplot(2,1,2),plot(tCoor,(angle(squeeze(IQ))),'k');
%     hold on, plot(tCoor,(angle(squeeze(sIQ))),'r');
%     hold on,plot(tCoor,(angle(squeeze(sIQHP))),'b');
%     title('phase(IQ)');% ylim([0,1])
%     xlabel('time lag, [ms]')
%     legend({'Raw','SVD','SVDHP'})
%     
% if handles.DataType==0 % phantom data select and processing
%     %% g1 analysis for phantom data
%     PRSSinfo.g1StartT=1;
%     PRSSinfo.g1nT=handles.nCC;
%     PRSSinfo.g1nTau=100;
%     PRSSinfo.rfnScale=1;
%     PRSSinfo.C=1540;                    % sound speed, m/s
%     PRSSinfo.FWHM=[125 100]*1e-6;        % (X, Z) spatial resolution, Full Width at Half Maximum of point spread function, m
%     PRSSinfo.rFrame=handles.fCC;               % sIQ frame rate, Hz
%     PRSSinfo.f0=16.625E6;
%     PRSSinfo.SVDrank=[3 PRSSinfo.g1nT];
%     handles.tauCoor=linspace(1,PRSSinfo.g1nTau,PRSSinfo.g1nTau)/handles.fCC; % s
%     
%     % gg=(handles.GG(handles.slt(1),handles.slt(2),:));
%     gg=sIQ2GG(sIQ, PRSSinfo);
%     PRSSinfo.Dim=size(gg);
%     [g1Vz]=GG2Vz(squeeze(gg).', PRSSinfo, 10)*1e3; % mm/s
%     Vcz=(ColorDoppler(sIQ,PRSSinfo)); % color Doppler, all frequency
%     [Mf, Vx, Vz, V, pVz, R, Ms, CR, gfit]=sIQ2vUS_SV_GPU(sIQ, PRSSinfo);
%     Fig=figure;
%     set(Fig,'Position',[300 400 1200 600]);
%     subplot(2,2,1);
%     plot(handles.tauCoor*1e3,squeeze(real(gg)),'.b')
%     hold on, plot(handles.tauCoor*1e3,squeeze(real(gfit)),'r')
%     xlabel('tau [ms]')
%     grid on
%     title(['GG-real, g1Vz0=', num2str(g1Vz),', Vz=', num2str(Vz),', Vcz=', num2str(Vcz), ' mm/s'])
%     
%     subplot(2,2,2);
%     plot(handles.tauCoor*1e3,squeeze(imag(gg)),'.b')
%     hold on, plot(handles.tauCoor*1e3, squeeze(imag(gfit)),'r')
%     xlabel('tau [ms]')
%     grid on
%     title(['GG-imag, Vt=', num2str(Vx),' mm/s, pVz=',num2str(pVz)])
%     
%     subplot(2,2,3);
%     plot(squeeze(gg),'.b')
%     hold on, plot(squeeze(gfit),'r')
%     xlim([-1 1]);ylim([-1 1])
%     grid on
%     title(['GG-complex, V=', num2str(V), 'mm/s; R=',num2str(R)])
%     
%     subplot(2,2,4);
%     plot(handles.tauCoor*1e3,squeeze(abs(gg)),'.b')
%     hold on; plot(handles.tauCoor*1e3,squeeze(abs(gfit)),'r')
%     xlabel('tau [ms]')
%     ylim([0 1])
%     grid on
%     title(['GG-magnitude, Ms=', num2str(Ms),', Mf=',num2str(Mf)])
% else % invivo data select and vUS processing
%     %% g1 analysis for invivo data
%     PRSSinfo.g1StartT=1;
%     PRSSinfo.g1nT=handles.nCC;
%     PRSSinfo.g1nTau=100;
%     PRSSinfo.rfnScale=1;
%     PRSSinfo.C=1540;                    % sound speed, m/s
%     PRSSinfo.FWHM=[125 90]*1e-6;        % (X, Z) spatial resolution, Full Width at Half Maximum of point spread function, m
%     PRSSinfo.rFrame=handles.fCC;               % sIQ frame rate, Hz
%     PRSSinfo.f0=16.625E6;
%     PRSSinfo.useMsk=0; % 1: use ULM data as spatial mask; 0: no spatial mask
%         
%     handles.tauCoor=linspace(1,PRSSinfo.g1nTau,PRSSinfo.g1nTau)/handles.fCC; % s
%     fIQ=(fft(sIQ,nt,3)); % no fft shift
%     for iNP=1:2
%         iFIQ=zeros(size(fIQ));
%         iFIQ(:,:,(floor(nt/2))*(iNP-1)+1:floor(nt/2)*iNP)=fIQ(:,:,(floor(nt/2))*(iNP-1)+1:floor(nt/2)*iNP);
%         iIQ=(ifft(iFIQ,nt,3));
%         gg(:,:,:,iNP)=sIQ2GG(iIQ, PRSSinfo); % g1 of p or n frequency signal
%     end
%     gg(:,:,:,3)=sIQ2GG(sIQ, PRSSinfo);
%     PRSSinfo.MpVz=0.8;
%     [Mf, Vz, V, pVz, Vcz, R, CR, Vx, Ms, pnRatio,gfit]=sIQ2vUS_NPDV_GPU(sIQ, PRSSinfo);
%     Fig=figure;
%     set(Fig,'Position',[300 400 1200 600]);
%     subplot(2,2,1);
%     plot(handles.tauCoor*1e3,squeeze(real(gg(:,:,:,3))),'.k')
%     hold on, plot(handles.tauCoor*1e3,squeeze(real(gfit(:,:,:,3))),'k')
%     hold on, plot(handles.tauCoor*1e3,squeeze(real(gg(:,:,:,1))),'.b')
%     hold on, plot(handles.tauCoor*1e3,squeeze(real(gfit(:,:,:,1))),'b')
%     hold on, plot(handles.tauCoor*1e3,squeeze(real(gg(:,:,:,2))),'.r')
%     hold on, plot(handles.tauCoor*1e3,squeeze(real(gfit(:,:,:,2))),'r')
%     xlabel('tau [ms]')
%     grid on
%     title(['GG-real, Vz=', num2str(Vz),' mm/s, [', num2str(handles.slt),']'])
%     
%     subplot(2,2,2);
%     plot(handles.tauCoor*1e3,squeeze(imag(gg(:,:,:,3))),'.k')
%     hold on, plot(handles.tauCoor*1e3,squeeze(imag(gfit(:,:,:,3))),'k')
%     hold on, plot(handles.tauCoor*1e3,squeeze(imag(gg(:,:,:,1))),'.b')
%     hold on, plot(handles.tauCoor*1e3,squeeze(imag(gfit(:,:,:,1))),'b')
%     hold on, plot(handles.tauCoor*1e3,squeeze(imag(gg(:,:,:,2))),'.r')
%     hold on, plot(handles.tauCoor*1e3,squeeze(imag(gfit(:,:,:,2))),'r')
%     xlabel('tau [ms]')
%     grid on
%     title(['GG-imag, Vt=', num2str(Vx),' mm/s, pVz=',num2str(pVz)])
%     
%     subplot(2,2,3);
%     plot(squeeze((gg(:,:,:,3))),'.k')
%     hold on, plot(squeeze((gfit(:,:,:,3))),'k')
%     hold on, plot(squeeze((gg(:,:,:,1))),'.b')
%     hold on, plot(squeeze((gfit(:,:,:,1))),'b')
%     hold on, plot(squeeze((gg(:,:,:,2))),'.r')
%     hold on, plot(squeeze((gfit(:,:,:,2))),'r')
%     xlim([-1 1]);ylim([-1 1])
%     grid on
%     title(['GG-complex, V=', num2str(V), 'mm/s; R=',num2str(R)])
%     
%     subplot(2,2,4);
%     plot(handles.tauCoor*1e3,squeeze(abs(gg(:,:,:,3))),'.k')
%     hold on, plot(handles.tauCoor*1e3,squeeze(abs(gfit(:,:,:,3))),'k')
%     hold on, plot(handles.tauCoor*1e3,squeeze(abs(gg(:,:,:,1))),'.b')
%     hold on, plot(handles.tauCoor*1e3,squeeze(abs(gfit(:,:,:,1))),'b')
%     hold on, plot(handles.tauCoor*1e3,squeeze(abs(gg(:,:,:,2))),'.r')
%     hold on, plot(handles.tauCoor*1e3,squeeze(abs(gfit(:,:,:,2))),'r')
%     xlabel('tau [ms]')
%     ylim([0 1])
%     grid on
%     title(['GG-magnitude, Ms=', num2str(Ms),', Mf=',num2str(Mf)])
%     
% end




% --- Executes on button press in BTN_Reset.
function BTN_Reset_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_Reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
addpath(handles.CodePath);
addpath([handles.CodePath, '\SubFunctions'])
clc;
handles.output = hObject;
handles.CodePath=pwd;
addpath(handles.CodePath);
addpath([handles.CodePath, '\SubFunctions'])
handles.DefPath='G:\PROJ-D-vUS';
handles.DefProbe='L22-14v';
handles.DefCtnFreq=18.5; % default center frequency, MHz
handles.DefPitch=0.1;    % default transducer element pitch, mm
handles.DefSpeed=1540;   % default sound speed, m/s or mm/ms
handles.DefCutFreq=25;   % default cutoff frequency for Power Doppler processing
handles.SVDRank=[25 1000];   % SVD rank
handles.newSlect=1;
handles.DataType=1; 
% Update handles structure
guidata(hObject, handles);

% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
addpath(handles.CodePath);
addpath([handles.CodePath, '\SubFunctions'])
clc;
%% Show GG at different iTau
set(hObject,'SliderStep',[1/(handles.PRSinfo.g1nTau-1), 3/(handles.PRSinfo.g1nTau-1)])
set(hObject,'Max',handles.PRSinfo.g1nTau)
iTau=handles.PRSinfo.g1nTau-min(round(get(hObject,'Value')),handles.PRSinfo.g1nTau-1);
set(handles.Input_iTau,'string',iTau);

axes(handles.axes1)
imagesc((abs((handles.GG(:,:,iTau))))); colorbar; axis image;
colormap(handles.axes1, jet)
caxis([0 1])
title(['Tau=',num2str(iTau/handles.fCC*1e3),'ms'])

%% Show IQ/sIQ at different t
% set(hObject,'SliderStep',[1/(handles.nCC-1), 3/(handles.nCC-1)])
% set(hObject,'Max',handles.nCC)
% iT=handles.nCC-min(round(get(hObject,'Value')),handles.nCC-1);
% set(handles.Input_iTau,'string',iT);
% 
% axes(handles.axes1)
% imagesc(log(abs((handles.IQ(:,:,iT))))); colorbar; 
% colormap(handles.axes1, gray)
% % caxis([0 1])
% title(['t=',num2str(iT/handles.fCC*1e3),'ms'])

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function Input_iTau_Callback(hObject, eventdata, handles)
% hObject    handle to Input_iTau (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Input_iTau as text
%        str2double(get(hObject,'String')) returns contents of Input_iTau as a double
iTau=str2num(get(handles.Input_iTau,'string'));
axes(handles.axes1)
imagesc((abs((handles.GG(:,:,iTau))))); colorbar; axis image;
colormap(handles.axes1, jet)
title(['iTau=',num2str(iTau)])
caxis([0 1])

% --- Executes during object creation, after setting all properties.
function Input_iTau_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Input_iTau (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Disp_FileName_Callback(hObject, eventdata, handles)
% hObject    handle to Disp_FileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Disp_FileName as text
%        str2double(get(hObject,'String')) returns contents of Disp_FileName as a double


% --- Executes during object creation, after setting all properties.
function Disp_FileName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Disp_FileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


