function varargout = GUI_CardiacRespiratory(varargin)
% GUI_CARDIACRESPIRATORY MATLAB code for GUI_CardiacRespiratory.fig
%      GUI_CARDIACRESPIRATORY, by itself, creates a new GUI_CARDIACRESPIRATORY or raises the existing
%      singleton*.
%
%      H = GUI_CARDIACRESPIRATORY returns the handle to a new GUI_CARDIACRESPIRATORY or the handle to
%      the existing singleton*.
%
%      GUI_CARDIACRESPIRATORY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_CARDIACRESPIRATORY.M with the given input arguments.
%
%      GUI_CARDIACRESPIRATORY('Property','Value',...) creates a new GUI_CARDIACRESPIRATORY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_CardiacRespiratory_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_CardiacRespiratory_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_CardiacRespiratory

% Last Modified by GUIDE v2.5 18-Jun-2018 12:10:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_CardiacRespiratory_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_CardiacRespiratory_OutputFcn, ...
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


% --- Executes just before GUI_CardiacRespiratory is made visible.
function GUI_CardiacRespiratory_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_CardiacRespiratory (see VARARGIN)

% Choose default command line output for GUI_CardiacRespiratory
handles.output = hObject;
addpath('D:\OneDrive\Work\PROJ - FUS\CODE\Functions');
addpath('D:\OneDrive\Work\PROJ - FUS\CODE\GUI');
handles.DefPath='H:\PJ - USI';
handles.DefProbe='L22-14v';
handles.DefCtnFreq=18.5; % default center frequency, MHz
handles.DefPitch=0.1;    % default transducer element pitch, mm
handles.DefSpeed=1540;   % default sound speed, m/s or mm/ms
handles.DefCutFreq=50;   % default cutoff frequency for Power Doppler processing
handles.iTau=20;
handles.Filt=0;
handles.GGSource=2;
handles.RankLow=30;   % SVD rank low cutoff
handles.RankHigh=180;  % SVD rank high cutoff
handles.nTau=50;
handles.nP=20;
handles.nPslt=5;
handles.DefnCC_proc=200;
handles.newSlt='Y';
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_CardiacRespiratory wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_CardiacRespiratory_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in BTN_load.
function BTN_load_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,datapath]=uigetfile(handles.DefPath);
handles.DefPath=datapath;
guidata(hObject, handles);
%% Load file and get file information
disp('Loading data...');
load ([datapath, filename]);
% IQ0=IQData{1};
% handles.IQ=squeeze(IQ0(:,:,1,1,:));
% cIQ=hilbert(squeeze(sum(IQ,3))); % Coherence Compounding IQ data, then Hilbert transform to get complext Coherence compounded IQ data (envelop+phase)
cIQ=IQ;
handles.IQ=cIQ; % 
%handles.Img=Img; % 
[handles.nz, handles.nx, handles.nt]=size(handles.IQ);
fileinfo=strsplit(filename(1:end-4),'-');
handles.Agl=str2num(fileinfo{2});
handles.nAgl=str2num(fileinfo{3});
handles.fCC=str2num(fileinfo{4});
handles.nCC=str2num(fileinfo{5});
handles.iPlane=str2num(fileinfo{8});
% handles.P.zCoor=linspace(P.zCoor(1),P.zCoor(1)+P.zCoor(end),length(P.zCoor)); 
% handles.P.xCoor=P.xCoor;
handles.P=P;
disp('Data Loaded!');
guidata(hObject, handles);

% --- Executes on button press in BTN_PDI.
function BTN_PDI_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_PDI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% addpath('D:\OneDrive\Work\PROJ - FUS\CODE\Functions');
% addpath('D:\OneDrive\Work\PROJ - FUS\CODE\GUI');
addpath('Z:\US-CODE\20190331\Functions');
addpath('Z:\US-CODE\20190331\GUI');
clc;
prompt={'SVD Rank (low):', ['SVD Rank (High):(Max Rank: ',num2str(handles.nCC),')'],...
    ['nCC_process (nCC total: ',num2str(handles.nCC),')']};
name='Power Doppler data processing';
defaultvalue={num2str(handles.RankLow), num2str(handles.RankHigh), num2str(min(handles.DefnCC_proc,handles.nCC))};
numinput=inputdlg(prompt,name, 1, defaultvalue);
handles.RankLow=str2num(numinput{1});
handles.RankHigh=str2num(numinput{2});
nCC_proc=str2num(numinput{3});
handles.DefnCC_proc=nCC_proc;
guidata(hObject, handles);
IQR=handles.IQ(:,:,1:nCC_proc);
% IQR=handles.Img(:,:,1:nCC_proc);
%% SVD process 2 (eigen-to-SVD use MATLAB)
rank=[handles.RankLow:handles.RankHigh];
[nz,nx,nt]=size(IQR);
S=reshape(IQR,[nz*nx,nt]);
S_COVt=(S'*S);
[V,D]=eig(S_COVt); % V is the right singular Vector of S/eigenvector; D is the eigenvalue/square of Singular value
for it=1:nt 
    Ddiag0(it)=abs(sqrt(D(it,it)));
end
Ddiag=20*log10(Ddiag0/max(Ddiag0)); % singular value in db
[Ddesc, Idesc]=sort(Ddiag0,'descend');
% figure,plot(Ddesc);
for it=1:nt
    Vdesc(:,it)=V(:,Idesc(it));
end
%% SVD singular value, temporal singular vector, and spatial singular vector check
U=S*Vdesc;
% %% plot singular value curve and PSD of temporal singular vector
% figure,plot(Ddesc)
% ylabel('Singular value dB')
% xlabel('Degree of singular vector')
% 
% FV=fftshift(fft(Vdesc,nt,1),1);
% fCoor=linspace(-handles.P.CCFR/2,handles.P.CCFR/2,nt);
% svdCoor=1:nt;
% figure,imagesc(svdCoor,fCoor,abs(FV))
% xlabel('iSingular vector')
% ylabel('f [HZ]')
% 
% FVD=abs(FV).*repmat(Ddesc,[nt 1]);
% FVDdb=20*log10(FVD/max(FVD(:)));
% figure,imagesc(svdCoor,fCoor,FVDdb)
% xlabel('iSingular vector')
% ylabel('f [HZ]')
% xlim([0 500])
% ylim([-500 500])
%% temporal singular vector
tCoor=linspace(1/handles.P.CCFR,nt/handles.P.CCFR,nt)*1e3;
iD=25;
figure,
subplot(2,1,1);plot(tCoor,real(Vdesc(:,iD)));
% xlabel('t [ms]')
ylabel('Re')
title(['iSVD=',num2str(iD)])
subplot(2,1,2);plot(tCoor,imag(Vdesc(:,iD)));
xlabel('t [ms]')
ylabel('Imag')
%% spatial singular vector

Ursp=reshape(U,[nz,nx,nt]);
figure,imagesc(handles.P.xCoor,handles.P.zCoor,abs(squeeze(Ursp(:,:,iD))))
title(['iD=',num2str(iD)])
xlabel('x [mm]');
ylabel('z [mm]');
axis equal tight
colorbar
%%
Vrank=zeros(size(Vdesc)); VrankBulk=Vrank;
Vrank(:,rank)=Vdesc(:,rank);
VrankBulk(:,1:handles.RankLow-1)=Vdesc(:,1:handles.RankLow-1);%1
Vnoise=zeros(size(Vdesc));
Vnoise(:,end)=Vdesc(:,min(300,nt));
UDelta=S*Vdesc;
sBlood0=reshape(UDelta*Vrank',[nz,nx,nt]);
sBulk=reshape(UDelta*VrankBulk',[nz,nx,nt]);
% sBlood=sBlood0./repmat(std(abs(sBlood0),1,2),[1,nx]);
%%%% Noise equalization 
sNoise=reshape(UDelta*Vnoise',[nz,nx,nt]);
sNoiseMed=medfilt2(abs(squeeze(mean(sNoise,3))),[30 30],'symmetric');
sNoiseMedNorm=sNoiseMed/min(sNoiseMed(:));
sBlood=sBlood0./repmat(sNoiseMedNorm,[1,1,nt]);
%% variance of blood signal
stdBlood=abs(std(sBlood,1,3));

% fig1=figure;
% set(fig1,'Position',[300 400 900 300]);
% subplot(1,2,1);imagesc(stdBlood);
% title(['std(sBlood),sBlood=SVD(IQ), Rank=[',num2str([handles.RankLow handles.RankHigh]),']'])
% caxis([median(stdBlood(:))/1.3, max(stdBlood(:))*0.5])
% colormap(hot);colorbar

%% power doppler
PDI=mean(abs(sBlood).^2,3); 
PDIdb=10*log10(PDI./max(PDI(:)));

axes(handles.axes1)
imagesc(PDIdb);
[nz,nx]=size(PDIdb);
caxis([min(PDIdb(:))*(1+sign(min(PDIdb(:)))*0.1), max(PDIdb(:))*(1-sign(max(PDIdb(:)))*0.2)]);
colormap(handles.axes1, hot);
colorbar;
xlim([1 nx])
ylim([1 nz])

PDIrfn=imresize(PDIdb,5);
zCoor=linspace(handles.P.zCoor(1),handles.P.zCoor(end),numel(handles.P.zCoor)*5);
xCoor=linspace(handles.P.xCoor(1),handles.P.xCoor(end),numel(handles.P.xCoor)*5);
figure,
imagesc(xCoor,zCoor,PDIdb);
caxis([min(PDIdb(:))*(1+sign(min(PDIdb(:)))*0.1), max(PDIdb(:))*(1-sign(max(PDIdb(:)))*0.3)]);
colormap(hot);
colorbar;
axis tight equal
title(['SVD Filtered'])
xlabel('x [mm]');
ylabel('z [mm]');
figure,
imagesc(PDIdb);
caxis([min(PDIdb(:))*(1+sign(min(PDIdb(:)))*0.1), max(PDIdb(:))*(1-sign(max(PDIdb(:)))*0.3)]);
colormap(hot);
colorbar;
axis tight equal
title(['SVD Filtered'])



handles.sBlood=sBlood0;
handles.sBulk=sBulk;
handles.PDI=PDIdb;
guidata(hObject, handles);

% --- Executes on button press in BTN_CRmapIQ.
function BTN_CRmapIQ_Callback(hObject, eventdata, handles)
% Calculate cardiac and respiratory map based on IQ (raw data)
% hObject    handle to BTN_CRmapIQ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% addpath('D:\OneDrive\Work\PROJ - FUS\CODE\Functions');
% addpath('D:\OneDrive\Work\PROJ - FUS\CODE\GUI');
addpath('Z:\US-CODE\20190331\Functions');
addpath('Z:\US-CODE\20190331\GUI');
clc;
IQ=handles.IQ;
[nz,nx,nt]=size(IQ);
tCoor=linspace(1/handles.fCC,nt/handles.fCC,nt)*1e3;
% 2.1 Fourier transfom of g1(iTau)

Lfft=2^nextpow2(nt);
FphaseIQ=(fft(squeeze(angle(IQ(:,:,:))-mean(angle(IQ(:,:,:)),3)),Lfft,3));
FrealIQ=(fft(squeeze(real(IQ(:,:,:))-mean(real(IQ(:,:,:)),3)),Lfft,3));
FabsIQ=(fft(squeeze(abs(IQ(:,:,:))-mean(abs(IQ(:,:,:)),3)),Lfft,3));
fCoor=linspace(0, handles.fCC/2,Lfft/2);
fCoorStep=mean(diff(fCoor));
% 2.2 locate the cardiac rate and respiratory rate
fig1=figure;
subplot(3,1,1);
plot(fCoor,abs(squeeze(mean(mean(FphaseIQ(:,:,1:Lfft/2),1),2))))
xlabel('Frequency [HZ]')
ylabel('Power (phase)')
subplot(3,1,2);
plot(fCoor,abs(squeeze(mean(mean(FrealIQ(:,:,1:Lfft/2),1),2))))
xlabel('Frequency [HZ]')
ylabel('Power (real)')
subplot(3,1,3);
plot(fCoor,abs(squeeze(mean(mean(FabsIQ(:,:,1:Lfft/2),1),2))))
xlabel('Frequency [HZ]')
ylabel('Power (magnitude)')

prompt={'Cardiac Frequency', 'Cardiac Frequency Range (CF+/- CFR)',...
    'Respiratory Frequency', 'Respiratory Frequency Range (RF+/- RFR)'};
name='CF and RF processing';
defaultvalue={'8','0','1.2','0'};
numinput=inputdlg(prompt,name, 1, defaultvalue);
CF=str2num(numinput{1});
CFR=str2num(numinput{2});
RF=str2num(numinput{3});
RFR=str2num(numinput{4});
close (fig1)


% 2.3 get cardiac and respiratory maps 
% CFrange=round((CF-CFR)/fCoorStep):round((CF+CFR)/fCoorStep);
% Cardiac=squeeze(max(abs(FabsIQ(:,:,CFrange)),[],3));
% RFrange=round((RF-RFR)/fCoorStep):round((RF+RFR)/fCoorStep);
% Resp=squeeze(max(abs(FrealIQ(:,:,RFrange)),[],3));
% % 2.4 plot cardiac and respiratory maps 
% fig=figure;
% set(fig,'Position',[200 400 1200 300])
% subplot(1,2,1);imagesc(handles.P.xCoor,handles.P.zCoor,Cardiac);
% axis tight equal
% title(['Cardiac map - IQ'])
% xlabel('x [mm]');
% ylabel('z [mm]');
% subplot(1,2,2);imagesc(handles.P.xCoor,handles.P.zCoor,Resp);
% axis tight equal
% title(['Respiratory map - IQ'])
% xlabel('x [mm]');
% ylabel('z [mm]');

%% 2.4 bandstop filtering to remove the cardiac signal
FimagIQ=(fft(squeeze(imag(IQ(:,:,:))-mean(imag(IQ(:,:,:)),3)),Lfft,3));
FrealIQ=(fft(squeeze(real(IQ(:,:,:))-mean(real(IQ(:,:,:)),3)),Lfft,3));
CFR=1;
% CFrange=[round((CF-CFR)/fCoorStep):round((CF+CFR)/fCoorStep),round((CF*2-CFR)/fCoorStep):round((CF*2+CFR)/fCoorStep)];
CFrange=[1:round((CF+CFR)/fCoorStep),round((CF*2-CFR)/fCoorStep):round((CF*2+CFR)/fCoorStep),round((CF*3-CFR)/fCoorStep):round((CF*3+CFR)/fCoorStep),round((CF*4-CFR)/fCoorStep):round((CF*4+CFR)/fCoorStep)];

FimagIQ(:,:,CFrange)=0;
IQrmvC=ifft(FimagIQ(:,:,1:512),1000,3);



guidata(hObject, handles);


% --- Executes on button press in BTN_g1.
function BTN_g1_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_g1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% addpath('D:\OneDrive\Work\PROJ - FUS\CODE\Functions');
% addpath('D:\OneDrive\Work\PROJ - FUS\CODE\GUI');
addpath('Z:\US-CODE\20190331\Functions');
addpath('Z:\US-CODE\20190331\GUI');
clc;
IQ=handles.IQ;

prompt={'ntau',...
    'tStart',...
    ['nt (nCC total: ',num2str(handles.nCC),')']};
name='g1 select and calculate';
defaultvalue={num2str(handles.nTau),...
    '1',...
    num2str(handles.nCC-handles.nTau)};
numinput=inputdlg(prompt,name, 1, defaultvalue);
nTau=str2num(numinput{1});
tStart=str2num(numinput{2});
nt=str2num(numinput{3});
disp('Calculating g1...')
handles.GG=IQ2g1(IQ,tStart,nt,nTau);
disp('g1 is calculated!')
guidata(hObject, handles);

% --- Executes on button press in BTN_CRmapG1.
function BTN_CRmapG1_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_CRmapG1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% addpath('D:\OneDrive\Work\PROJ - FUS\CODE\Functions');
% addpath('D:\OneDrive\Work\PROJ - FUS\CODE\GUI');
addpath('Z:\US-CODE\20190331\Functions');
addpath('Z:\US-CODE\20190331\GUI');
clc;

GG=handles.GG;
[nz,nx,nTau]=size(GG);
tCoor=linspace(1/handles.fCC,nTau/handles.fCC,nTau)*1e3;
% 2.1 Fourier transfom of g1(iTau)
Lfft=2^nextpow2(nTau);
FphaseGG=(fft(squeeze(angle(GG(:,:,:))-mean(angle(GG(:,:,:)),3)),Lfft,3));
FrealGG=(fft(squeeze(real(GG(:,:,:))-mean(real(GG(:,:,:)),3)),Lfft,3));
FabsGG=(fft(squeeze(abs(GG(:,:,:))-mean(abs(GG(:,:,:)),3)),Lfft,3));
fCoor=linspace(0, handles.fCC/2,Lfft/2);
fCoorStep=mean(diff(fCoor));
% 2.2 locate the cardiac rate and respiratory rate
fig1=figure;
subplot(3,1,1);
plot(fCoor,abs(squeeze(mean(mean(FphaseGG(:,:,1:Lfft/2),1),2))))
xlabel('Frequency [HZ]')
ylabel('Power (phase)')
subplot(3,1,2);
plot(fCoor,abs(squeeze(mean(mean(FrealGG(:,:,1:Lfft/2),1),2))))
xlabel('Frequency [HZ]')
ylabel('Power (real)')
subplot(3,1,3);
plot(fCoor,abs(squeeze(mean(mean(FabsGG(:,:,1:Lfft/2),1),2))))
xlabel('Frequency [HZ]')
ylabel('Power (magnitude)')

prompt={'Cardiac Frequency', 'Cardiac Frequency Range (CF+/- CFR)',...
    'Respiratory Frequency', 'Respiratory Frequency Range (RF+/- RFR)'};
name='CF and RF processing';
defaultvalue={'8','0','1.2','0'};
numinput=inputdlg(prompt,name, 1, defaultvalue);
CF=str2num(numinput{1});
CFR=str2num(numinput{2});
RF=str2num(numinput{3});
RFR=str2num(numinput{4});
close (fig1)
% 2.3 get cardiac and respiratory maps 
CFrange=ceil((CF-CFR)/fCoorStep):ceil((CF+CFR)/fCoorStep);
Cardiac=squeeze(max(abs(FabsGG(:,:,CFrange)),[],3));
RFrange=ceil((RF-RFR)/fCoorStep):ceil((RF+RFR)/fCoorStep);
Resp=squeeze(max(abs(FphaseGG(:,:,RFrange)),[],3));
% 2.4 plot cardiac and respiratory maps 
fig=figure;
set(fig,'Position',[200 400 1200 300])
subplot(1,2,1);imagesc(handles.P.xCoor,handles.P.zCoor,Cardiac);
axis tight equal
title(['Cardiac map - GG'])
xlabel('x [mm]');
ylabel('z [mm]');
subplot(1,2,2);imagesc(handles.P.xCoor,handles.P.zCoor,Resp);
axis tight equal
title(['Respiratory map - GG'])
xlabel('x [mm]');
ylabel('z [mm]');

% --- Executes on button press in BTN_xRossCorr.
function BTN_xRossCorr_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_xRossCorr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% addpath('D:\OneDrive\Work\PROJ - FUS\CODE\Functions');
% addpath('D:\OneDrive\Work\PROJ - FUS\CODE\GUI');
addpath('Z:\US-CODE\20190331\Functions');
addpath('Z:\US-CODE\20190331\GUI');
clc;
IQ0=handles.IQ;
IQ=IQ0(107:110,46,:);
[nz,nx,nt]=size(IQ(:,:,:));
refIQ=squeeze(IQ(:,:,1));
for it=1:nt
    iIQ=squeeze(IQ(:,:,it));
    xRoss(it)=sum(iIQ(:).*conj(refIQ(:)))/(sqrt(sum(abs(iIQ(:)).^2))*sqrt(sum(abs(refIQ(:)).^2)));
end
tCoor=linspace(1/handles.fCC,nt/handles.fCC,nt)*1e3;
figure;
subplot(1,2,1);
plot(tCoor,abs(xRoss))
xlabel('t [ms]');
ylabel('magnitude')
subplot(1,2,2);
plot(tCoor,(angle(xRoss)))
xlabel('t [ms]');
ylabel('Phase')

%% g1
nTau=200;
g1=squeeze(IQ2g1(IQ,1,nt,nTau));
tCoorG1=linspace(1/handles.fCC,nTau/handles.fCC,nTau)*1e3;
figure;
subplot(1,2,1);
plot(tCoorG1,abs(g1))
xlabel('t [ms]');
ylabel('magnitude, g1')
subplot(1,2,2);
plot(tCoorG1,(angle(g1)))
xlabel('t [ms]');
ylabel('Phase,g1')

handles.xRoss=xRoss;
guidata(hObject, handles);


% --- Executes on button press in BTN_IQroi.
function BTN_IQroi_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_IQroi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% addpath('D:\OneDrive\Work\PROJ - FUS\CODE\Functions');
% addpath('D:\OneDrive\Work\PROJ - FUS\CODE\GUI');
addpath('Z:\US-CODE\20190331\Functions');
addpath('Z:\US-CODE\20190331\GUI');
clc;
PDI=handles.PDI;
IQ=handles.sBulk;%
[nz,nx,nt]=size(IQ);
tCoor=linspace(1,nt,nt)/handles.fCC*1000;
axes(handles.axes1)
imagesc(PDI)
caxis(MyCaxis(PDI));
%% select ROI to plot IQ(t) and cross correlation
ButtonName = questdlg('Select New ROI?', ...
        'Select', ...
        'Yes', 'No', 'Cancel');
switch ButtonName
    case 'Yes'
        [loc_x, loc_y]=ginput(6);  % no image rotation,
        BW=roipoly(PDI,loc_x,loc_y);
        handles.BW=BW;
    case 'No'
        BW = handles.BW;
end
DPI_roi=PDI.*BW;
axes(handles.axes1)
hold on
OlapPlot = imagesc(DPI_roi);colormap(jet);
alpha(OlapPlot,double((DPI_roi)>-25))
IQROI=IQ.*repmat(BW,[1, 1,nt]);
% dataROI=IQ(40:140,25:95,:);
% dataROI=IQ;
meanROI=squeeze(mean(mean(IQROI,1),2));

refIQ=squeeze(IQROI(:,:,1));
tic
for it=1:nt
    iIQ=squeeze(IQROI(:,:,it));
    xRossROI(it)=sum(iIQ(:).*conj(refIQ(:)))/(sqrt(sum((iIQ(:)).^2))*sqrt(sum((refIQ(:)).^2)));
end
toc

figure,
subplot(2,1,1);
yyaxis left
plot(tCoor,abs(meanROI))
xlabel('t [ms]')
ylabel('Magnitude IQ');
yyaxis right
plot(tCoor,abs(xRossROI))
ylabel('Magnitude xCoef')
%xlim([0 1000])
subplot(2,1,2);
yyaxis left
plot(tCoor,angle(meanROI))
xlabel('t [ms]')
ylabel('Phase IQ');
yyaxis right
plot(tCoor,(angle(xRossROI)))
xlabel('t [ms]')
ylabel('Phase xCoef')
%xlim([0 1000])

% figure,
% subplot(2,1,1)
% plot(tCoor,abs(meanROI));
% xlabel('t [ms]')
% ylabel('magnitude')
% title('IQ')
% subplot(2,1,2)
% plot(tCoor,angle(meanROI));
% xlabel('t [ms]')
% ylabel('phase')
% 
% 
% figure,
% subplot(2,1,1);
% plot(tCoor,abs(xRossROI))
% xlabel('t [ms]')
% ylabel('Magnitude')
% title('Cross correlation Coeff')
% subplot(2,1,2);
% plot(tCoor,angle((xRossROI)))
% xlabel('t [ms]')
% ylabel('Angle')
% title('Cross correlation Coeff')

%% IQ and cross correlation coefficient at different depths
% loc_z=40:20:140;
% HW=5;
% for idepth=1:length(loc_z)
%     IQ_depth(:,:,:,idepth)=IQ(loc_z(idepth)-HW:loc_z(idepth)+HW,150:200,:);%50:110
%     IQ_depth_mean(:,idepth)=squeeze(mean(mean(IQ_depth(:,:,:,idepth),1),2));
%     IQ_depth_mean_norm(:,idepth)=(IQ_depth_mean(:,idepth)-min(IQ_depth_mean(:,idepth)))./(max(IQ_depth_mean(:,idepth))-min(IQ_depth_mean(:,idepth)));
%     
%     refIQ=squeeze(IQ_depth(:,:,1,idepth));
%     for it=1:nt
%         iIQ=squeeze(IQ_depth(:,:,it,idepth));
%         xRoss(it,idepth)=sum(iIQ(:).*conj(refIQ(:)))/(sqrt(sum((iIQ(:)).^2))*sqrt(sum((refIQ(:)).^2)));
%     end
% 
%     figure,
%     subplot(2,1,1);
%     yyaxis left
%     plot(tCoor,abs(IQ_depth_mean(:,idepth)))
%     xlabel('t [ms]')
%     ylabel('Magnitude IQ');
%     yyaxis right
%     plot(tCoor,abs(xRoss(:,idepth)))
%     ylabel('Magnitude xCoef')
%     title(['IQ, iDepth=',num2str(loc_z(idepth)-HW),'-',num2str(loc_z(idepth)+HW)])
%     %xlim([0 1000])
%     subplot(2,1,2);
%     yyaxis left
%     plot(tCoor,angle(IQ_depth_mean(:,idepth)))
%     xlabel('t [ms]')
%     ylabel('Phase IQ');
%     yyaxis right
%     plot(tCoor,(angle(xRoss(:,idepth))))
%     xlabel('t [ms]')
%     ylabel('Phase xCoef')
%     %xlim([0 1000])
%     
%     figure,
%     subplot(2,1,1);
%     plot(tCoor,abs(xRoss(:,idepth)))
%     xlabel('t [ms]')
%     ylabel('Magnitude')
%     title('Cross correlation Coeff')
%     title(['Cross correlation Coeff, iDepth=',num2str(loc_z(idepth)-HW),'-',num2str(loc_z(idepth)+HW)])
%     subplot(2,1,2);
%     plot(tCoor,abs(real(xRoss(:,idepth))))
%     xlabel('t [ms]')
%     ylabel('Angle')
%     title('Cross correlation Coeff')
% end

% figure,
% subplot(2,1,1);plot(tCoor,abs(IQ_depth_mean))
% xlabel('t [ms]')
% ylabel('Magnitude IQ');
% subplot(2,1,2);plot(tCoor,angle(IQ_depth_mean))
% xlabel('t [ms]')
% ylabel('Phase IQ');
% % figure,plot(tCoor,abs(IQ_depth_mean_norm))
% % xlabel('t [ms]')
% % ylabel('Magnitude normalized IQ');
% 
% for idepth=1:length(loc_z)
%     refIQ=squeeze(IQ_depth(:,:,1,idepth));
%     for it=1:nt
%         iIQ=squeeze(IQ_depth(:,:,it,idepth));
%         xRoss(it,idepth)=sum(iIQ(:).*conj(refIQ(:)))/(sqrt(sum(abs(iIQ(:)).^2))*sqrt(sum(abs(refIQ(:)).^2)));
%     end
% end
% figure,
% subplot(2,1,1);
% plot(tCoor,abs(xRoss))
% xlabel('t [ms]')
% ylabel('Magnitude')
% title('Cross correlation Coeff')
% subplot(2,1,2);
% plot(tCoor,abs(real(xRoss)))
% xlabel('t [ms]')
% ylabel('Angle')
% title('Cross correlation Coeff')

guidata(hObject, handles);
% --- Executes on button press in BTN_CRrmv.
function BTN_CRrmv_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_CRrmv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% addpath('D:\OneDrive\Work\PROJ - FUS\CODE\Functions');
% addpath('D:\OneDrive\Work\PROJ - FUS\CODE\GUI');
addpath('Z:\US-CODE\20190331\Functions');
addpath('Z:\US-CODE\20190331\GUI');
clc;
PDI=handles.PDI;
IQ=handles.IQ;
[nz,nx,nt]=size(IQ);

axes(handles.axes1)
imagesc(PDI)
xlim([1 nx])
ylim([1 nz])
caxis(MyCaxis(PDI));
%% select ROI
% select ROI to plot IQ(t) and cross correlation
% [loc_x, loc_y]=ginput(6);  % 
% BW=roipoly(PDI,loc_x,loc_y);
% handles.BW=BW;
% DPI_roi=PDI.*BW;
% axes(handles.axes1)
% hold on
% OlapPlot = imagesc(DPI_roi);colormap(jet);
% alpha(OlapPlot,double((DPI_roi)>-25))
%% 0 single pixel comparision 
% IQROI0=(IQ(1:160,21:250,:));
% IQsvd=handles.sBlood(1:160,21:250,:);
% IQbulkSvd=handles.sBulk(1:160,21:250,:);
IQROI0=IQ;
IQsvd=handles.sBlood;
IQbulkSvd=handles.sBulk;

[nz,nx,nt0]=size(IQROI0);
tCoor=linspace(1,nt,nt)/handles.fCC*1000;
% %% 0.0 select pixel 
[xSlt, zSlt]=ginput(1);  % 
iz=round(zSlt)-0
ix=round(xSlt)-20
% iz=17;
% ix=49;
%% 0.1 SVD results
iIQ0=IQROI0(iz,ix,:);
SVDrsd=IQsvd(iz,ix,:);
SVDrmv=IQbulkSvd(iz,ix,:);
SVDrmvSpec=fftshift(fft(SVDrmv,nt0,3));
%% 0.2 data prep for FFT and HP filtering
nTextra=floor(40e-3*handles.P.CCFR);
iIQ=zeros(1,1,nt+2*nTextra);
iIQ(:,:,1:nTextra)=2*iIQ0(:,:,1)-flip(iIQ0(:,:,1:nTextra),3);
iIQ(:,:,nTextra+1:end-nTextra)=iIQ0;
iIQ(:,:,end-nTextra+1:end)=2*iIQ0(:,:,end)-flip(iIQ0(:,:,end-nTextra+1:end),3);
[nz,nx,nt]=size(iIQ);
%% 0.3 FFT-based filtering
FIQ0=(fft(((iIQ)-mean((iIQ),3)),nt,3));
fCoor=linspace(0,handles.P.CCFR/2,nt/2);
fCoorStep=mean(diff(fCoor));
CF=8;
CFR=10;
% CFrange=[round((CF-CFR)/fCoorStep):round((CF+CFR)/fCoorStep),round((CF*2-CFR)/fCoorStep):round((CF*2+CFR)/fCoorStep)];
% CFrange=[1:round((CF+CFR)/fCoorStep),round((CF*2-CFR)/fCoorStep):round((CF*2+CFR)/fCoorStep),round((CF*3-CFR)/fCoorStep):round((CF*3+CFR)/fCoorStep),round((CF*4-CFR)/fCoorStep):round((CF*4+CFR)/fCoorStep)];
% CFrange=[round((CF-CFR)/fCoorStep):round((CF+CFR)/fCoorStep),round((CF*2)/fCoorStep):round((CF*2)/fCoorStep)];
% CFrange=[round((CF-CFR)/fCoorStep):round((CF+CFR)/fCoorStep)];
% CFrange=[1:round((CF+CFR)/fCoorStep)];
CFrange=[1:round((CF+CFR)/fCoorStep),round((CF*2)/fCoorStep)+[-1:1],round((CF*3)/fCoorStep)+[-2:1],...
    round((CF*4)/fCoorStep)+[-2:1],round((CF*5)/fCoorStep)+[-2:1], round((CF*6)/fCoorStep)+[-2:1],...
    round((CF*7)/fCoorStep)+[-2:1],round((CF*8)/fCoorStep)+[-2:1]];

CFrangeAll=[CFrange,nt-CFrange+2];
CFrangeAll(CFrangeAll>nt)=[];

FIQ=zeros(size(FIQ0));
FIQ(:,:,CFrangeAll)=FIQ0(:,:,CFrangeAll);
temp=ifft(FIQ,nt,3);
FFTrmv=temp(:,:,nTextra+1:end-nTextra)+mean((iIQ));

FIQ=FIQ0;
FIQ(:,:,CFrangeAll)=0;
temp=ifft(FIQ,nt,3);
FFTrsd=temp(:,:,nTextra+1:end-nTextra);
FFTrmvSpec=fftshift(fft(FFTrmv,nt0,3));
%% 0.3 Butterworth HP
fCut=20;
[B,A]=butter(2,fCut/(handles.fCC/2),'high');    %coefficients for the high pass filter
iIQfilt0=filter(B,A,iIQ,[],3);    % blood signal (filtering in the time dimension)
HPrsd=iIQfilt0(:,:,nTextra+1:end-nTextra);
HPrmv=iIQ0-HPrsd;
HPrmvSpec=fftshift(fft(HPrmv,nt0,3));

iIQ0Spec=fftshift(fft(iIQ0,nt0,3));
fCoor=linspace(-handles.P.CCFR/2,handles.P.CCFR/2,nt0);
%% 0.4 plot and compare
% Real and imaginary
% fig=figure;
% set(fig,'Position',[200 400 1200 600])
% axReal(1)=subplot(2,3,1);
% plot(tCoor,squeeze(real(iIQ0)))
% hold on, plot(tCoor,squeeze(real(SVDrmv)))
% ylabel('real')
% xlabel('t [ms]')
% legend({'Raw','SVD'})
% title('RAW vs SVDrmv')
% axReal(2)=subplot(2,3,2);
% plot(tCoor,squeeze(real(iIQ0)))
% hold on, plot(tCoor,squeeze(real(FFTrmv)))
% ylabel('real')
% xlabel('t [ms]')
% legend({'Raw','FFT'})
% title('RAW vs FFTrmv')
% axReal(3)=subplot(2,3,3);
% plot(tCoor,squeeze(real(iIQ0)))
% hold on, plot(tCoor,squeeze(real(HPrmv)))
% ylabel('real')
% xlabel('t [ms]')
% legend({'Raw','HP'})
% title('RAW vs HPrmv')
% 
% axImag(1)=subplot(2,3,4);
% plot(tCoor,squeeze(imag(iIQ0)))
% hold on, plot(tCoor,squeeze(imag(SVDrmv)))
% ylabel('imag')
% xlabel('t [ms]')
% legend({'Raw','SVDrmv'})
% 
% axImag(2)=subplot(2,3,5);
% plot(tCoor,squeeze(imag(iIQ0)))
% hold on, plot(tCoor,squeeze(imag(FFTrmv)))
% ylabel('imag')
% xlabel('t [ms]')
% legend({'Raw','FFTrmv'})
% axImag(3)=subplot(2,3,6);
% plot(tCoor,squeeze(imag(iIQ0)))
% hold on, plot(tCoor,squeeze(imag(HPrmv)))
% ylabel('imag')
% xlabel('t [ms]')
% legend({'Raw','HPrmv'})
% linkaxes(axReal)
% linkaxes(axImag)

% magnitude and phase
fig=figure;
set(fig,'Position',[200 400 1200 600])
axabs(1)=subplot(2,3,1);
plot(tCoor,squeeze(abs(iIQ0)))
hold on, plot(tCoor,squeeze(abs(SVDrmv)))
ylabel('Magnitude')
xlabel('t [ms]')
legend({'Raw','SVD'})
title('RAW vs SVDrmv')
axabs(2)=subplot(2,3,2);
plot(tCoor,squeeze(abs(iIQ0)))
hold on, plot(tCoor,squeeze(abs(FFTrmv)))
ylabel('Magnitude')
xlabel('t [ms]')
legend({'Raw','FFT'})
title('RAW vs FFTrmv')
axabs(3)=subplot(2,3,3);
plot(tCoor,squeeze(abs(iIQ0)))
hold on, plot(tCoor,squeeze(abs(HPrmv)))
ylabel('Magnitude')
xlabel('t [ms]')
legend({'Raw','HP'})
title('RAW vs HPrmv')

axangle(1)=subplot(2,3,4);
plot(tCoor,squeeze(angle(iIQ0)))
hold on, plot(tCoor,squeeze(angle(SVDrmv)))
ylabel('Phase')
xlabel('t [ms]')
legend({'Raw','SVDrmv'})

axangle(2)=subplot(2,3,5);
plot(tCoor,squeeze(angle(iIQ0)))
hold on, plot(tCoor,squeeze(angle(FFTrmv)))
ylabel('Phase')
xlabel('t [ms]')
legend({'Raw','FFTrmv'})
axangle(3)=subplot(2,3,6);
plot(tCoor,squeeze(angle(iIQ0)))
hold on, plot(tCoor,squeeze(angle(HPrmv)))
ylabel('Phase')
xlabel('t [ms]')
legend({'Raw','HPrmv'})
linkaxes(axabs)
linkaxes(axangle)

% spectrum of removded signal
fig=figure;
set(fig,'Position',[200 400 1200 600])
axabs(1)=subplot(4,1,1);
plot(fCoor,squeeze(abs(iIQ0Spec)))
ylabel('Power')
xlabel('f [Hz]')
title('Spectrum raw')

axabs(2)=subplot(4,1,2);
plot(fCoor,squeeze(abs(SVDrmvSpec)))
ylabel('Power')
xlabel('f [Hz]')
title('Spectrum SVD')

axabs(3)=subplot(4,1,3);
plot(fCoor,squeeze(abs(FFTrmvSpec)))
ylabel('Power')
xlabel('f [Hz]')
title('Spectrum FFT')

axabs(4)=subplot(4,1,4);
plot(fCoor,squeeze(abs(HPrmvSpec)))
ylabel('Power')
xlabel('f [Hz]')
title('Spectrum HP')
linkaxes(axabs)

% plot residual signal
fig=figure;
set(fig,'Position',[300 200 1200 600])
axRsdabs(1)=subplot(2,3,1);
plot(tCoor,squeeze(abs(SVDrsd)))
ylabel('magnitude')
xlabel('t [ms]')
title('SVD')
axRsdabs(2)=subplot(2,3,2);
 plot(tCoor,squeeze(abs(FFTrsd)))
ylabel('magnitude')
xlabel('t [ms]')
title('FFT')
axRsdabs(3)=subplot(2,3,3);
 plot(tCoor,squeeze(abs(HPrsd)))
ylabel('magnitude')
xlabel('t [ms]')
title('HP')
axRsdangle(1)=subplot(2,3,4);
 plot(tCoor,squeeze(angle(SVDrsd)))
ylabel('Phase')
xlabel('t [ms]')
title('SVD')
axRsdangle(2)=subplot(2,3,5);
 plot(tCoor,squeeze(angle(FFTrsd)))
ylabel('Phase')
xlabel('t [ms]')
title('FFT')
axRsdangle(3)=subplot(2,3,6);
 plot(tCoor,squeeze(angle(HPrsd)))
ylabel('Phase')
xlabel('t [ms]')
title('HP')
linkaxes(axRsdabs)
linkaxes(axRsdangle)
%% 1 cardiac and respiratory remove using fft-based filtering
% IQROI0=(IQ(30:160,51:200,1:1000));
% [nz,nx,nt]=size(IQROI0);
% tCoor=linspace(1,nt,nt)/handles.fCC*1000;
% IQROI=zeros(nz,nx,nt+100);
% IQROI(:,:,1:50)=2*IQROI0(:,:,1)-flip(IQROI0(:,:,1:50),3);
% IQROI(:,:,51:end-50)=IQROI0;
% IQROI(:,:,end-49:end)=2*IQROI0(:,:,end)-flip(IQROI0(:,:,end-49:end),3);
% [nz,nx,nt]=size(IQROI);
% %% 1.1 single pixel check
% [xSlt, zSlt]=ginput(1);  % 
% iz=round(zSlt)-29
% ix=round(xSlt)-50
% % iz=19;
% % ix=33;
% FIQ0(iz,ix,:)=(fft(((IQROI(iz,ix,:))-mean((IQROI(iz,ix,:)),3)),nt,3));
% FrealIQ(iz,ix,:)=(fft((real(IQROI(iz,ix,:))-mean(real(IQROI(iz,ix,:)),3)),nt,3));
% fCoor=linspace(0,handles.P.CCFR/2,nt/2);
% fCoorStep=mean(diff(fCoor));
% CFR=10;
% CF=8;
% % CFrange=[round((CF-CFR)/fCoorStep):round((CF+CFR)/fCoorStep),round((CF*2-CFR)/fCoorStep):round((CF*2+CFR)/fCoorStep)];
% % CFrange=[1:round((CF+CFR)/fCoorStep),round((CF*2-CFR)/fCoorStep):round((CF*2+CFR)/fCoorStep),round((CF*3-CFR)/fCoorStep):round((CF*3+CFR)/fCoorStep),round((CF*4-CFR)/fCoorStep):round((CF*4+CFR)/fCoorStep)];
% % CFrange=[round((CF-CFR)/fCoorStep):round((CF+CFR)/fCoorStep),round((CF*2)/fCoorStep):round((CF*2)/fCoorStep)];
% % CFrange=[round((CF-CFR)/fCoorStep):round((CF+CFR)/fCoorStep)];
% % CFrange=[1:round((CF+CFR)/fCoorStep)];
% CFrange=[1:round((CF+CFR)/fCoorStep),round((CF*2)/fCoorStep)+[-1:1],round((CF*3)/fCoorStep)+[-2:1],...
%     round((CF*4)/fCoorStep)+[-2:1],round((CF*5)/fCoorStep)+[-2:1], round((CF*6)/fCoorStep)+[-2:1],...
%     round((CF*7)/fCoorStep)+[-2:1],round((CF*8)/fCoorStep)+[-2:1]];
% 
% CFrangeAll=[CFrange,nt-CFrange+2];
% CFrangeAll(CFrangeAll>nt)=[];
% FIQ=zeros(size(FIQ0));
% FIQ(:,:,CFrangeAll)=FIQ0(:,:,CFrangeAll);
% % FIQ=(FIQ0);
% % FIQ(:,:,CFrangeAll)=0;
% temp=ifft(FIQ,nt,3);
% IQrmvC=temp(:,:,51:end-50);
% 
% % figure,
% % subplot(3,1,1);
% % plot(abs(squeeze(FIQ0(iz,ix,:))))
% % hold on, plot(abs(squeeze(FIQ(iz,ix,:))))
% % ylabel('Spectrum')
% % subplot(3,1,2)
% % plot(tCoor,squeeze(real(IQROI0(iz,ix,:)))-mean(squeeze(real(IQROI0(iz,ix,:)))))
% % hold on, plot(tCoor,squeeze(real(IQrmvC(iz,ix,:))))
% % ylabel('real')
% % legend({'Raw','Filt'})
% % subplot(3,1,3)
% % plot(tCoor,squeeze(imag(IQROI0(iz,ix,:)))-mean(squeeze(imag(IQROI0(iz,ix,:)))))
% % hold on, plot(tCoor,squeeze(imag(IQrmvC(iz,ix,:))))
% % ylabel('imag')
% % legend({'Raw','Filt'})
% 
% figure
% subplot(2,1,1)
% plot(tCoor,squeeze(real(IQROI0(iz,ix,:)))-squeeze(real(IQrmvC(iz,ix,:)))-mean(squeeze(real(IQROI0(iz,ix,:)))-squeeze(real(IQrmvC(iz,ix,:)))))
% ylabel('real')
% legend({'residual'})
% title('FFT filter')
% subplot(2,1,2)
% plot(tCoor,squeeze(imag(IQROI0(iz,ix,:)))-squeeze(imag(IQrmvC(iz,ix,:)))-mean(squeeze(imag(IQROI0(iz,ix,:)))-squeeze(imag(IQrmvC(iz,ix,:)))))
% ylabel('imag')
% legend({'residual'})
%% 1.2 whole volume check
% FIQ0=(fft(((IQROI)-mean((IQROI),3)),nt,3));
% fCoor=linspace(0,handles.P.CCFR/2,nt/2);
% fCoorStep=mean(diff(fCoor));
% CFR=10;
% CF=8;
% CFrange=[1:round((CF+CFR)/fCoorStep),round((CF*2)/fCoorStep)+[-1:1],round((CF*3)/fCoorStep)+[-2:1],...
%     round((CF*4)/fCoorStep)+[-2:1],round((CF*5)/fCoorStep)+[-2:1], round((CF*6)/fCoorStep)+[-2:1],...
%     round((CF*7)/fCoorStep)+[-2:1],round((CF*8)/fCoorStep)+[-2:1]];
% 
% CFrangeAll=[CFrange,nt-CFrange+2];
% CFrangeAll(CFrangeAll>nt)=[];
% % FIQ=zeros(size(FIQ0));
% % FIQ(:,:,CFrangeAll)=FIQ0(:,:,CFrangeAll);
% FIQ=(FIQ0);
% FIQ(:,:,CFrangeAll)=0;
% temp=ifft(FIQ,nt,3);
% IQrmv=temp(:,:,21:end-20);
%  
% PDI=mean(abs(IQrmv).^2,3);
% PDIdb=10*log10(PDI./max(PDI(:)));
% PDIrfn=imresize(PDIdb,5);
% figure
% imagesc(PDIdb);
% caxis([min(PDIdb(:))*(1+sign(min(PDIdb(:)))*0.1), max(PDIdb(:))*(1-sign(max(PDIdb(:)))*0.3)]);
% colormap(hot)
%% 2 high pass based removal
% fCut=50;
% %% 2.1 single pixel check
% % [xSlt, zSlt]=ginput(1);  % 
% % iz=round(zSlt)-29
% % ix=round(xSlt)-50
% % iz=19;
% % ix=33;
% 
% iIQ0=IQROI(iz,ix,:);
% [B,A]=butter(4,fCut/(handles.fCC/2),'high');    %coefficients for the high pass filter
% iIQfilt0=filter(B,A,iIQ0,[],3);    % blood signal (filtering in the time dimension)
% iIQ=iIQ0(:,:,51:end-50);
% iIQfilt=iIQfilt0(:,:,51:end-50);

% figure,
% subplot(2,1,1)
% plot(tCoor,squeeze(real(iIQ)))
% hold on, plot(tCoor,squeeze(real(iIQ))-squeeze(real(iIQfilt)))
% ylabel('real')
% legend({'Raw','HPFilt'})
% subplot(2,1,2)
% plot(tCoor,squeeze(imag(iIQ)))
% hold on, plot(tCoor,squeeze(imag(iIQ))-squeeze(imag(iIQfilt)))
% ylabel('imag')
% legend({'Raw','HPFilt'})

% figure,
% subplot(2,1,1)
%  plot(tCoor,squeeze(real(iIQfilt)))
% ylabel('real')
% title('Butterworth HP')
% subplot(2,1,2)
% plot(tCoor,squeeze(imag(iIQfilt)))
% ylabel('imag')
% %% 3 SVD based removal
% IQsvd=handles.sBlood(30:160,51:200,1:1000);
% IQbulkSvd=handles.sBulk(30:160,51:200,1:1000);
% %% 3.1 single pixel check
% % [xSlt, zSlt]=ginput(1);  % 
% % iz=round(zSlt)-29
% % ix=round(xSlt)-50
% % iz=19;
% % ix=33;
% 
% iIQ0=IQROI0(iz,ix,:);
% iIQsvd=IQsvd(iz,ix,:);
% iIQsvdRmv=IQbulkSvd(iz,ix,:);

% figure,
% subplot(2,1,1)
% plot(tCoor,squeeze(real(iIQ0)))
% hold on, plot(tCoor,squeeze(real(iIQsvdRmv)))
% ylabel('real')
% legend({'Raw','SVD'})
% subplot(2,1,2)
% plot(tCoor,squeeze(imag(iIQ0)))
% hold on, plot(tCoor,squeeze(imag(iIQsvdRmv)))
% ylabel('imag')
% legend({'Raw','SVD'})
% 
% figure,
% subplot(2,1,1)
%  plot(tCoor,squeeze(real(iIQsvd)))
% ylabel('real')
% title('SVD')
% subplot(2,1,2)
% plot(tCoor,squeeze(imag(iIQsvd)))
% ylabel('imag')

%% 4 cardiac and respiratory remove using localized moving window averaging
% IQROI=(IQ(30:160,51:200,1:100));
% [nz,nx,nt]=size(IQROI);
% tCoor=linspace(1,nt,nt)/handles.fCC*1000;
% Lwindow=2;
% 
% % IQmavg=MoveAvg(IQROI,Lwindow,3);
% IQmavg=movmean(IQROI,Lwindow,3);
% IQrmv=IQROI-IQmavg;
% 
% [xSlt, zSlt]=ginput(1);  % 
% iz=round(zSlt)-29
% ix=round(xSlt)-50
% 
% zSlt=iz;
% xSlt=ix;
% fig=figure;
% set(fig,'Position',[400 200 500 700])
% subplot(2,2,1)
% plot(tCoor,abs(squeeze(IQROI(zSlt,xSlt,:))))
% hold on; plot(tCoor,abs(squeeze(IQmavg(zSlt,xSlt,:))))
% xlabel('t [ms]')
% ylabel('magnitude')
% legend({'RAW','MAVG'})
% subplot(2,2,2)
% plot(tCoor,angle(squeeze(IQROI(zSlt,xSlt,:))))
% hold on; plot(tCoor,angle(squeeze(IQmavg(zSlt,xSlt,:))))
% xlabel('t [ms]')
% ylabel('phase')
% 
% subplot(2,2,3)
% plot(tCoor,real(squeeze(IQROI(zSlt,xSlt,:))))
% hold on; plot(tCoor,real(squeeze(IQmavg(zSlt,xSlt,:))))
% xlabel('t [ms]')
% ylabel('real')
% legend({'RAW','MAVG'})
% % ylabel('real')
% % ylabel('real')
% subplot(2,2,4)
% plot(tCoor,imag(squeeze(IQROI(zSlt,xSlt,:))))
% hold on; plot(tCoor,imag(squeeze(IQmavg(zSlt,xSlt,:))))
% xlabel('t [ms]')
% ylabel('image')
% 
% figure
% subplot(3,1,1)
% plot(tCoor,real(squeeze(IQROI(zSlt,xSlt,:)))-real(squeeze(IQmavg(zSlt,xSlt,:))))
% xlabel('t [ms]')
% ylabel('real')
% legend({'MAVG'})
% % ylabel('real')
% % ylabel('real')
% subplot(3,1,2)
% plot(tCoor,imag(squeeze(IQROI(zSlt,xSlt,:)))-imag(squeeze(IQmavg(zSlt,xSlt,:))))
% xlabel('t [ms]')
% ylabel('image')
% subplot(3,1,3)
% plot(tCoor,abs(squeeze(IQROI(zSlt,xSlt,:)))-abs(squeeze(IQmavg(zSlt,xSlt,:))))
% xlabel('t [ms]')
% ylabel('abs')
% 
% PDI=mean(abs(IQrmv).^2,3);
% PDIdb=10*log10(PDI./max(PDI(:)));
% PDIrfn=imresize(PDIdb,5);
% figure
% imagesc(PDIdb);
% caxis([min(PDIdb(:))*(1+sign(min(PDIdb(:)))*0.1), max(PDIdb(:))*(1-sign(max(PDIdb(:)))*0.3)]);
% colormap(hot)
%% Cardiac and respiratory remove using regression
%% global phase regression
% IQROI=(IQ(30:160,51:200,1:100));
% [nz,nx,nt]=size(IQROI);
% tCoor=linspace(1,nt,nt)/handles.fCC*1000;
% % refIQ=squeeze(IQROI(:,:,1));
% % for it=1:nt
% %     iIQ=squeeze(IQROI(:,:,it));
% %     xCoefRaw(it)=sum(iIQ(:).*conj(refIQ(:)))/(sqrt(sum((iIQ(:)).^2))*sqrt(sum((refIQ(:)).^2)));
% % end
% % handles.xCoefRaw=xCoefRaw;
% % guidata(hObject, handles);
% xCoefRaw=handles.xCoefRaw;
% % % PLOT xCoef
% tCoor=linspace(1,nt,nt)/handles.fCC*1000;
% % figure;
% % subplot(2,1,1);
% % plot(tCoor,real(xCoefRaw));
% % ylabel('xCoef real')
% % xlabel('t [ms]');
% % subplot(2,1,2);
% % plot(tCoor,imag(xCoefRaw));
% % ylabel('xCoef Imag')
% % xlabel('t [ms]');
% 
% % figure;
% % subplot(2,1,1);
% % plot(tCoor,abs(xCoefRaw));
% % ylabel('xCoef magnitude')
% % xlabel('t [ms]');
% % subplot(2,1,2);
% % plot(tCoor,angle(xCoefRaw));
% % ylabel('xCoef Phase')
% % xlabel('t [ms]');
% 
% 
% % Caridac phase compensation
% Greal=shiftdim(real(xCoefRaw)-mean(real(xCoefRaw)));
% GImag=shiftdim(imag(xCoefRaw)-mean(imag(xCoefRaw)));
% GPhase=shiftdim(angle(xCoefRaw)-mean(angle(xCoefRaw)));
% [xSlt, zSlt]=ginput(1);  % 
% iz=round(zSlt)-29
% ix=round(xSlt)-50
% 
% % for iz=1:nz
% %     for ix=1:nx
%         ipixIQ=squeeze(angle(IQROI(iz,ix,:)));
%         ipixIQdiff=diff(ipixIQ);
%         if max(abs(ipixIQdiff))>6
%             clear Iwrap
%             Ind=find((abs(ipixIQdiff))>6);
%             if mod(numel(Ind),2)~=0
%                 Ind(numel(Ind)+1)=length(ipixIQ);
%             end
%             Iwrap(:,1)=Ind(1:2:end)+1;
%             Iwrap(:,2)=Ind(2:2:end);
%             [nChange,~]=size(Iwrap);
%             for iChange=1:nChange
%                 ipixIQ(Iwrap(iChange,1):Iwrap(iChange,2))=ipixIQ(Iwrap(iChange,1):Iwrap(iChange,2))-sign(ipixIQdiff(Ind(1)))*2*pi;
%             end
%         end
%         IQROIphase(iz,ix,:)=ipixIQ;
%         ipixIQmean=mean((ipixIQ));
%         ipixIQ0=ipixIQ-ipixIQmean;
%         %% fit phase
%         fitE=@(c) sum(abs(c(1)*(GPhase)-ipixIQ0).^2);
%         fitC0=(max(ipixIQ0)-min(ipixIQ0))/(max(GPhase)-min(GPhase));
%         [fitC, fVal]=fminsearch(fitE,fitC0,optimset('Display','off'));
%         CReal(iz,ix,1)=fitC(1);
%         IQbulkPhase(iz,ix,:)=fitC(1)*GPhase+ipixIQmean;
%         IQphaseCp(iz,ix,:)=(IQROIphase(iz,ix,:))-IQbulkPhase(iz,ix,:);
%         IQROIphaseCP(iz,ix,:)=IQROI(iz,ix,:).*exp(-1i*IQbulkPhase(iz,ix,:));
% %     end
% % end
% % IQROIphaseCP=IQROI.*exp(-1i*IQbulkPhase);
% zSlt=iz;
% xSlt=ix;
% % figure,
% % subplot(2,1,1)
% % plot(tCoor,(squeeze(IQROIphase(zSlt,xSlt,:))))
% % hold on; plot(tCoor,(squeeze(IQbulkPhase(zSlt,xSlt,:))))
% % xlabel('t [ms]')
% % ylabel('Phase')
% % legend({'RAW','Fit'})
% % % ylabel('real')
% % subplot(2,1,2)
% % plot(tCoor,(squeeze(IQphaseCp(zSlt,xSlt,:))))
% % xlabel('t [ms]')
% % ylabel('phase CPed')
% figure,
% subplot(2,1,1)
% plot(tCoor,real(squeeze(IQROI(zSlt,xSlt,:)))-mean(real(squeeze(IQROI(zSlt,xSlt,:)))))
% hold on; plot(tCoor,real(squeeze(IQROIphaseCP(zSlt,xSlt,:)))-mean(real(squeeze(IQROIphaseCP(zSlt,xSlt,:)))))
% xlabel('t [ms]')
% ylabel('real')
% legend({'RAW','CP'})
% % ylabel('real')
% subplot(2,1,2)
% plot(tCoor,imag(squeeze(IQROI(zSlt,xSlt,:)))-mean(imag(squeeze(IQROI(zSlt,xSlt,:)))))
% hold on; plot(tCoor,imag(squeeze(IQROIphaseCP(zSlt,xSlt,:)))-mean(imag(squeeze(IQROIphaseCP(zSlt,xSlt,:)))))
% xlabel('t [ms]')
% ylabel('image')
% 
% figure,
% subplot(2,1,1)
% plot(tCoor,abs(squeeze(IQROI(zSlt,xSlt,:)))-mean(abs(squeeze(IQROI(zSlt,xSlt,:)))))
% hold on; plot(tCoor,abs(squeeze(IQROIphaseCP(zSlt,xSlt,:)))-mean(abs(squeeze(IQROIphaseCP(zSlt,xSlt,:)))))
% xlabel('t [ms]')
% ylabel('magnitude')
% legend({'RAW','CP'})
% % ylabel('real')
% subplot(2,1,2)
% plot(tCoor,angle(squeeze(IQROI(zSlt,xSlt,:)))-mean(angle(squeeze(IQROI(zSlt,xSlt,:)))))
% hold on; plot(tCoor,angle(squeeze(IQROIphaseCP(zSlt,xSlt,:)))-mean(angle(squeeze(IQROIphaseCP(zSlt,xSlt,:)))))
% xlabel('t [ms]')
% ylabel('phase')
%% subROI-based real regression
% nLateralROI=2;
% nPixPerLateralROI=nx/nLateralROI;
% nPixPerAxialROI=2;
% nAxialROI=ceil(nz/nPixPerAxialROI);
% for iLateral=1:nLateralROI
%     for iAxial=1:nAxialROI
%         izROIstart=(iAxial-1)*nPixPerAxialROI+1;
%         izROIend=min(iAxial*nPixPerAxialROI,nz);
%         ixROIstart=(iLateral-1)*nPixPerLateralROI+1;
%         ixROIend=min(iLateral*nPixPerLateralROI,nx);
%         iIQROI=IQROI(izROIstart:izROIend,ixROIstart:ixROIend,:);
%         iROIabs=squeeze(mean(mean(abs(iIQROI),1),2));
%         iROI0abs=iROIabs-mean(iROIabs);
%         for ix=ixROIstart:ixROIend
%             for iz=izROIstart:izROIend
%                 ipixIQ=squeeze(real(IQROI(iz,ix,:)));
%                 ipixIQmean=mean((ipixIQ));
%                 ipixIQ0=ipixIQ-ipixIQmean;
%                 %% fit magnitude
%                 ipixIQ0abs=(ipixIQ0);
%                 fitE=@(c) sum(abs(c(1)*(iROI0abs)-ipixIQ0abs).^2);
%                 fitC0=(max(ipixIQ0abs)-min(ipixIQ0abs))/(max(iROI0abs)-min(iROI0abs));
%                 [fitC, fVal]=fminsearch(fitE,fitC0,optimset('Display','off'));
%                 CReal(iz,ix,1)=fitC(1);
%                 IQbulkMag(iz,ix,:)=fitC(1)*iROI0abs+ipixIQmean;
%                 
% %                 %% fit real
% %                 ipixIQ0Real=real(ipixIQ0);
% %                 fitE=@(c) sum(abs(c(1)*(iROI0Abs)-ipixIQ0Real).^2);
% %                 fitC0=(max(ipixIQ0Real)-min(ipixIQ0Real))/(max(iROI0Abs)-min(iROI0Abs));
% %                 [fitC, fVal]=fminsearch(fitE,fitC0,optimset('Display','off'));
% %                 CReal(iz,ix,1)=fitC(1);
% %                 IQbulkReal(iz,ix,:)=fitC(1)*iROI0Abs;
% %                 %% fit imaginary 
% %                 ipixIQ0Imag=imag(ipixIQ0);
% %                 fitE=@(c) sum(abs(c(1)*(iROI0Imag)-ipixIQ0Imag).^2);
% %                 fitC0=(max(ipixIQ0Imag)-min(ipixIQ0Imag))/(max(iROI0Imag)-min(iROI0Imag));
% %                 [fitC, fVal]=fminsearch(fitE,fitC0,optimset('Display','off'));
% %                 CImag(iz,ix,1)=fitC(1);
% %                 IQbulkImag(iz,ix,:)=fitC(1)*iROI0Imag;
%                 
%             end
%         end
%     end
% end
% IQresid=abs(IQROI)-IQbulkMag;
% 
% zSlt=iz;
% xSlt=ix;
% figure,
% subplot(2,1,1)
% plot(tCoor,abs(squeeze(IQROI(zSlt,xSlt,:))))
% hold on; plot(tCoor,(squeeze(IQbulkMag(zSlt,xSlt,:))))
% xlabel('t [ms]')
% ylabel('magnitude')
% % ylabel('real')
% subplot(2,1,2)
% plot(tCoor,angle(squeeze(IQROI(zSlt,xSlt,:))))
% hold on; plot(tCoor,angle(squeeze(IQbulkMag(zSlt,xSlt,:))))
% xlabel('t [ms]')
% ylabel('phase')

% figure,
% subplot(2,1,1)
% plot(tCoor,real(squeeze(IQ(zSlt,xSlt,:)-IQbulkReal(zSlt,xSlt,:)))+real(ipixIQmean))
% xlabel('t [ms]')
% ylabel('magnitude')
% % ylabel('real')
% subplot(2,1,2)
% plot(tCoor,angle(squeeze(IQ(zSlt,xSlt,:))))
% xlabel('t [ms]')
% ylabel('phase')
        
%% 
% figure,
% subplot(2,1,1)
% plot(tCoor,real(meanROI))
% xlabel('t [ms]')
% % ylabel('magnitude')
% ylabel('real')
% subplot(2,1,2)
% plot(tCoor,angle(meanROI))
% xlabel('t [ms]')
% ylabel('phase')

% --- Executes on button press in BTN_BSI.
function BTN_BSI_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_BSI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% addpath('D:\OneDrive\Work\PROJ - FUS\CODE\Functions');
% addpath('D:\OneDrive\Work\PROJ - FUS\CODE\GUI');
addpath('Z:\US-CODE\20190331\Functions');
addpath('Z:\US-CODE\20190331\GUI');
clc;
IQ0=handles.sBulk;% handles.sBulk
P=handles.P;
[~,~,nt]=size(IQ0);
ButtonName = questdlg('Select New ROI?', ...
        'Select', ...
        'Yes', 'No', 'Cancel');
switch ButtonName
    case 'Yes'
        [loc_x, loc_y]=ginput(1);  % no image rotation,
        Loc = [round(loc_y),round(loc_x)];
        handles.Loc=Loc;
    case 'No'
        Loc = handles.Loc;
end
 IQ=IQ0(Loc(1)-2: Loc(1)+2, Loc(2)-2:Loc(2)+2, :);
% ROI selection reference
% 3*3 roi 0.02s*9000 = 3min; 10*10 roi 0.033s*9000 = 5min; 30*30 roi
% 0.048s*9000 = 7min;
% IQ=abs(IQ); % image shift compensation for image intensity/magnitude
[nz,nx,nt]=size(IQ);

nRfnImg= 20;%10 single pixel 50um --> single pixel 5 um
nRfnCoef=30;
nOrgPix=0.4; %0.4 maximum detectable shift 0.4*10*5um == 0.4*50um = 20 um 
IQRfn=imresize(squeeze(IQ(:,:,round(nt/2))),nRfnImg);
RefIQ=IQRfn(nRfnImg+1:end-nRfnImg,nRfnImg+1:end-nRfnImg);

tic
for it=1:nt

    iIQ0=squeeze(IQ(:,:,it));
    iIQrfn=imresize(iIQ0,nRfnImg);
    for iShiftDx=-nOrgPix*nRfnImg:nOrgPix*nRfnImg
%         iStartX=max(iShiftDx,1);iEndX=min(nx*nRfnImg+iShiftDx,nx*nRfnImg);
%         iNx=iEndX-iStartX+1;
        %tic
        for iShiftDz=-nOrgPix*nRfnImg:nOrgPix*nRfnImg
            iIQ=iIQrfn(nRfnImg+iShiftDz+1:nRfnImg+iShiftDz+(nz-2)*nRfnImg,nRfnImg+1+iShiftDx:nRfnImg+iShiftDx+(nx-2)*nRfnImg);
            Xcoef(iShiftDz+nOrgPix*nRfnImg+1,iShiftDx+nOrgPix*nRfnImg+1,it)=sum(iIQ(:).*conj(RefIQ(:)))/(sqrt(sum(abs(iIQ(:)).^2))*sqrt(sum(abs(RefIQ(:)).^2)));
        end
        %toc
    end
%     %%
%     iIQ0=squeeze(IQ(:,:,it));
%     iIQrfn=imresize(iIQ0,nRfnImg);
%     for iShiftDx=-nOrgPix*nRfnImg:nOrgPix*nRfnImg
%         iStartX=max(iShiftDx,1);iEndX=min(nx*nRfnImg+iShiftDx,nx*nRfnImg);
%         iNx=iEndX-iStartX+1;
%         for iShiftDz=-nOrgPix*nRfnImg:nOrgPix*nRfnImg
%             iIQ1=circshift(iIQrfn,[iShiftDz,iShiftDx]);
%             iStartZ=max(iShiftDz,1);iEndZ=min(nz*nRfnImg+iShiftDz,nz*nRfnImg);
%             iNz=iEndZ-iStartZ+1;
%            
% %             iRefIQ=refIQ(1:iNz,1:iNx);
%             iRefIQ=refIQ(iStartZ:iEndZ,iStartX:iEndX);
%             iIQ=iIQ1(iStartZ:iEndZ,iStartX:iEndX);
%             Xcoef(iShiftDz+nOrgPix*nRfnImg+1,iShiftDx+nOrgPix*nRfnImg+1,it)=sum(iIQ(:).*conj(iRefIQ(:)))/(sqrt(sum(abs(iIQ(:)).^2))*sqrt(sum(abs(iRefIQ(:)).^2)));
%         end
%     end
    
end
XcoefOrg=squeeze(Xcoef(nOrgPix*nRfnImg+1,nOrgPix*nRfnImg+1,:));
XcoefMax=squeeze(max(max(Xcoef,[],1),[],2));
tCoor=linspace(1,nt,nt)/handles.fCC*1000;
toc
figure;
subplot(2,1,1);
plot(tCoor,abs(XcoefOrg));
hold on,plot(tCoor,abs(XcoefMax))
xlabel('t [ms]');
ylabel('magnitude')
legend({'XCoef at 0 shift'; 'XCoefMax' })
subplot(2,1,2);
plot(tCoor,angle(XcoefOrg));
hold on,plot(tCoor,angle(XcoefMax))
xlabel('t [ms]');
ylabel('Phase')
legend({'XCoef at 0 shift'; 'XCoefMax' })

% plot xCoef map
iXcoef=abs(squeeze(Xcoef(:,:,1)));
zCoorCoef=[-nOrgPix*nRfnImg:nOrgPix*nRfnImg]*P.dzImg*1e3/nRfnImg;
xCoorCoef=[-nOrgPix*nRfnImg:nOrgPix*nRfnImg]*P.dxImg*1e3/nRfnImg;
figure,imagesc(zCoorCoef,xCoorCoef,iXcoef)
colormap(jet)
colorbar
xlabel('X shift [um]')
ylabel('Z shift [um]')
axis equal tight
%caxis([0.95,1]);

%% refine cross correlation map and find the image shift
zCoefRfn=linspace(-nOrgPix*nRfnImg,nOrgPix*nRfnImg,(2*nOrgPix*nRfnImg+1)*nRfnCoef);
xCoefRfn=linspace(-nOrgPix*nRfnImg,nOrgPix*nRfnImg,(2*nOrgPix*nRfnImg+1)*nRfnCoef);
tic
[nc,~,~] = size(Xcoef);
XcoefRfn = zeros([nc*nRfnCoef, nc*nRfnCoef, nt]);
for it=1:nt
    iXcoef=abs(squeeze(Xcoef(:,:,it)));
    iXcoefFilt=imgaussfilt(iXcoef,2);
    iXcoefRfn=imresize(iXcoefFilt,nRfnCoef);
%     
%     figure,imagesc(zCoefRfn,xCoefRfn,iXcoefRfn)
%     colormap(jet)
%     colorbar
%     xlabel('X shift [pix]')
%     ylabel('Z shift [pix]')
%     axis equal tight
%     % the summit point
%     Temp=zeros(size(xRossTrfn));
%     Temp(xRossTrfn>max(xRossTrfn(:))-0.005)=xRossTrfn(xRossTrfn>max(xRossTrfn(:))-0.005);
%     figure,imagesc(z,x,a)
    [zMax(it),xMax(it)]=find(iXcoefRfn==max(iXcoefRfn(:)));
    XcoefRfn(:,:,it) = iXcoefRfn; 
end
toc

iXcoef=abs(squeeze(XcoefRfn(:,:,1)));
zCoorCoef= linspace(-nOrgPix*nRfnImg*P.dzImg*1e3/nRfnImg,nOrgPix*nRfnImg*P.dzImg*1e3/nRfnImg,nc*nRfnCoef );
xCoorCoef= linspace(-nOrgPix*nRfnImg*P.dxImg*1e3/nRfnImg,nOrgPix*nRfnImg*P.dxImg*1e3/nRfnImg,nc*nRfnCoef);
figure,imagesc(zCoorCoef,xCoorCoef,iXcoef)
colormap(jet)
colorbar
xlabel('X shift [um]')
ylabel('Z shift [um]')
axis equal tight
%caxis([0.95,1]);


zShiftPix=zMax-(nOrgPix*nRfnImg+1/2)*nRfnCoef;
xShiftPix=xMax-(nOrgPix*nRfnImg+1/2)*nRfnCoef;
zShiftPeak = zCoorCoef(zMax);
xShiftPeak = xCoorCoef(xMax);
% zShiftPeak=zCoefRfn(zMax)*P.dzImg*1e3/nRfnImg;
% xShiftPeak=xCoefRfn(xMax)*P.dxImg*1e3/nRfnImg;

% figure;
% subplot(2,1,1);
% plot(tCoor,zShiftPix);
% xlabel('t [ms]');
% ylabel('zShift, [Pix]')
% title('Image Shift')
% subplot(2,1,2);
% plot(tCoor,xShiftPix);
% xlabel('t [ms]');
% ylabel('xShift, [Pix]')

% % plot video
% addpath D:\US_TissueMotion\Functions
% ImgInfo.savePathName = ['D:\US_TissueMotion\', 'XcoefRfn_ROI15.mat'];
% ImgInfo.t = 500/1000; 
% make3Dmovie_singleROI(abs(XcoefRfn),zShiftPix,xShiftPix,ImgInfo);

% plot xCoef and recovered image shift
figure;
subplot(2,1,1);
yyaxis left
plot(tCoor,abs(XcoefOrg))
xlabel('t [ms]');
ylabel('xCoef, magnitude')
yyaxis right
plot(tCoor,sqrt(abs(zShiftPeak).^2+abs(xShiftPeak).^2));
ylabel('total Shift, [um]')

subplot(2,1,2);
yyaxis left
plot(tCoor,angle(XcoefOrg))
xlabel('t [ms]');
ylabel('xCoef, Phase')
yyaxis right
plot(tCoor,zShiftPeak);
hold on, plot(tCoor,xShiftPeak);
legend({'xCoef','z','x'})
ylabel('z and x Shift, [um]')

guidata(hObject, handles);


%% calculate all ROIs
clear zShiftPeak xShiftPeak;
ROIx=[50:50:200];%5
ROIz=[30:40:150];
nxROI=length(ROIx);
nzROI=length(ROIz);
for ixROI=1:nxROI-1
    for izROI=1:nzROI-1
        IQ=IQ0(ROIz(izROI):ROIz(izROI)+5,ROIx(ixROI):ROIx(ixROI)+5,:);
        % IQ=abs(IQ); % image shift compensation for image intensity/magnitude
        [nz,nx,nt]=size(IQ);
        nRfnImg=20;
        nRfnCoef=30;
        nOrgPix=0.4;
        IQRfn=imresize(squeeze(IQ(:,:,1)),nRfnImg);
        RefIQ=IQRfn(nRfnImg+1:end-nRfnImg,nRfnImg+1:end-nRfnImg);
        for it=1:nt
            iIQ0=squeeze(IQ(:,:,it));
            iIQrfn=imresize(iIQ0,nRfnImg);
            for iShiftDx=-nOrgPix*nRfnImg:nOrgPix*nRfnImg
                iStartX=max(iShiftDx,1);iEndX=min(nx*nRfnImg+iShiftDx,nx*nRfnImg);
                iNx=iEndX-iStartX+1;
                for iShiftDz=-nOrgPix*nRfnImg:nOrgPix*nRfnImg
                    iIQ=iIQrfn(nRfnImg+1+iShiftDz:nRfnImg+iShiftDz+(nz-2)*nRfnImg,nRfnImg+1+iShiftDx:nRfnImg+iShiftDx+(nx-2)*nRfnImg);
                    Xcoef(iShiftDz+nOrgPix*nRfnImg+1,iShiftDx+nOrgPix*nRfnImg+1,it)=sum(iIQ(:).*conj(RefIQ(:)))/(sqrt(sum(abs(iIQ(:)).^2))*sqrt(sum(abs(RefIQ(:)).^2)));
                end
            end
        end
        % XcoefOrg=squeeze(Xcoef(nOrgPix*nRfnImg+1,nOrgPix*nRfnImg+1,:));
        XcoefOrg=squeeze(max(max(Xcoef,[],1),[],2));
        % plot xCoef map
        iXcoef=abs(squeeze(Xcoef(:,:,1)));
        zCoorCoef=[-nOrgPix*nRfnImg:nOrgPix*nRfnImg]*P.dzImg*1e3/nRfnImg;
        xCoorCoef=[-nOrgPix*nRfnImg:nOrgPix*nRfnImg]*P.dxImg*1e3/nRfnImg;
        figure,imagesc(zCoorCoef,xCoorCoef,iXcoef)
        colormap(jet)
        colorbar
        xlabel('X shift [um]')
        ylabel('Z shift [um]')
        axis equal tight
        %% refine cross correlation map and find the image shift
        zCoefRfn=linspace(-nOrgPix*nRfnImg,nOrgPix*nRfnImg,(2*nOrgPix*nRfnImg+1)*nRfnCoef);
        xCoefRfn=linspace(-nOrgPix*nRfnImg,nOrgPix*nRfnImg,(2*nOrgPix*nRfnImg+1)*nRfnCoef);
        for it=1:100
            iXcoef=abs(squeeze(Xcoef(:,:,it)));
            iXcoefFilt=imgaussfilt(iXcoef,2);
            iXcoefRfn=imresize(iXcoefFilt,nRfnCoef);
            [zMax(it),xMax(it)]=find(iXcoefRfn==max(iXcoefRfn(:)));
        end
        zShiftPix=zMax-(nOrgPix*nRfnImg+1/2)*nRfnCoef;
        xShiftPix=xMax-(nOrgPix*nRfnImg+1/2)*nRfnCoef;
        zShiftPeak(izROI,ixROI,:)=zCoefRfn(zMax)*P.dzImg*1e3/nRfnImg;
        xShiftPeak(izROI,ixROI,:)=xCoefRfn(xMax)*P.dxImg*1e3/nRfnImg;
        
        zShiftPeak0 = squeeze(zShiftPeak(izROI,ixROI,:));
        xShiftPeak0 = squeeze(xShiftPeak(izROI,ixROI,:));
        
        figure; 
        subplot(2,1,1);
        yyaxis left
        plot(tCoor,abs(XcoefOrg))
        xlabel('t [ms]');
        ylabel('xCoef, magnitude')
        yyaxis right
        plot(tCoor,sqrt(abs(zShiftPeak0).^2+abs(xShiftPeak0).^2));
        ylabel('total Shift, [um]')

        subplot(2,1,2);
        yyaxis left
        plot(tCoor,angle(XcoefOrg))
        xlabel('t [ms]');
        ylabel('xCoef, Phase')
        yyaxis right
        plot(tCoor,zShiftPeak0);
        hold on, plot(tCoor,xShiftPeak0);
        legend({'xCoef','z','x'})
        ylabel('z and x Shift, [um]')
    end
    ixROI
end

%% compensate for the image shift
IQ=IQ0(51:80,50:110,:);
[nz,nx,nt]=size(IQ);
IQc=zeros(nz,nx,nt);
nPixPerDepth=30;
nDepth=floor(nz/nPixPerDepth);
tic
for it=1:100
    for iDepth=1:nDepth
        iIQ0=squeeze(IQ((iDepth-1)*nPixPerDepth+1:min(iDepth*nPixPerDepth,nz),:,it));
        iIQrfn=imresize(iIQ0,nRfnImg*nRfnCoef);
        iIQshift=circshift(iIQrfn,[zShiftPix(it),xShiftPix(it)]);
        iIQc=imresize(iIQshift,[nPixPerDepth,nx]);
        IQc((iDepth-1)*nPixPerDepth+1:min(iDepth*nPixPerDepth,nz),:,it)=iIQc;
    end
end
toc

%% calculate and compare xCoef between before and after image shift compensated 
% before compensation
IQxCoef=IQ0(51:140,50:180,:);
[nz,nx,nt]=size(IQxCoef);
refIQ=squeeze(IQxCoef(:,:,1));
for it=1:100
    iIQ=squeeze(IQxCoef(:,:,it));
    xCoefRaw(it)=sum(iIQ(:).*conj(refIQ(:)))/(sqrt(sum(abs(iIQ(:)).^2))*sqrt(sum(abs(refIQ(:)).^2)));
end
% after compensation
IQxCoef=IQc(:,:,:);
% IQxCoef=IQxCoef(:,:,1:100).*repmat(reshape(exp(-1i*angle(xCoefRaw)),[1 1 numel(xCoefRaw)]),[nz,nx, 1]);


refIQ=squeeze(IQxCoef(:,:,1));
for it=1:100
    iIQ=squeeze(IQxCoef(:,:,it));
    xCoefCmp(it)=sum(iIQ(:).*conj(refIQ(:)))/(sqrt(sum(abs(iIQ(:)).^2))*sqrt(sum(abs(refIQ(:)).^2)));
end
figure;
subplot(2,1,1);
plot(tCoor,abs(xCoefRaw));
hold on, plot(tCoor,abs(xCoefCmp));
xlabel('t [ms]');
ylabel('xCoef magnitude')
legend({'Raw', 'Shift Cmp'})
subplot(2,1,2);
plot(tCoor,angle(xCoefRaw));
hold on,plot(tCoor,angle(xCoefCmp));
xlabel('t [ms]');
ylabel('xCoef Phase')
legend({'Raw', 'Shift Cmp'})

% %%
% 
%     
% iXcoef=abs(squeeze(Xcoef(:,:,14)));
% % xRossTfilt=medfilt2(xRossT,[5 5]);
% iXcoefFilt=imgaussfilt(iXcoef,2);
% iXcoefRfn=imresize(iXcoefFilt,10);
% zCoorCoef=linspace(-nOrgPix*nRfnImg,nOrgPix*nRfnImg,(2*nOrgPix*nRfnImg+1)*10);
% xCoorCoef=linspace(-nOrgPix*nRfnImg,nOrgPix*nRfnImg,(2*nOrgPix*nRfnImg+1)*10);
% 
% figure,imagesc(zCoorCoef,xCoorCoef,iXcoefRfn)
% colormap(jet)
% colorbar
% xlabel('X shift [pix]')
% ylabel('Z shift [pix]')
% axis equal tight
% 
% % tCoor=linspace(1/handles.fCC,nt/handles.fCC,nt)*1e3;
% 
% %% global phase fluctuation
% clear iIQ refIQ;
% refIQ=squeeze(IQ(:,:,1));
% for it=1:100
%     iIQ=squeeze(IQ(:,:,it));
%     GPF(it)=-angle(sum(iIQ(:).*conj(refIQ(:))));
%     xRossGPF(it)=sum(iIQ(:)*exp(1i*GPF(it)).*conj(refIQ(:)))/(sqrt(sum(abs(iIQ(:)).^2))*sqrt(sum(abs(refIQ(:)).^2)));
% end
% handles.xRoss=Xcoef;
% guidata(hObject, handles);
% figure;
% subplot(2,1,1);
% plot(tCoor,abs(squeeze(Xcoef(nOrgPix*nRfnImg+1,nOrgPix*nRfnImg+1,:))));
% hold on,plot(tCoor,abs(xRossGPF))
% xlabel('t [ms]');
% ylabel('magnitude')
% legend({'XCoef at 0 shift'; 'XCoefGPF' })
% subplot(2,1,2);
% plot(tCoor,real(squeeze(Xcoef(nOrgPix*nRfnImg+1,nOrgPix*nRfnImg+1,:))));
% hold on,plot(tCoor,real(xRossGPF))
% xlabel('t [ms]');
% ylabel('Phase')
% legend({'XCoef at 0 shift'; 'XCoefGPF' })


