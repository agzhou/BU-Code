function varargout = GUI_SVDAnalysis(varargin)
% GUI_SVDANALYSIS MATLAB code for GUI_SVDAnalysis.fig
%      GUI_SVDANALYSIS, by itself, creates a new GUI_SVDANALYSIS or raises the existing
%      singleton*.
%
%      H = GUI_SVDANALYSIS returns the handle to a new GUI_SVDANALYSIS or the handle to
%      the existing singleton*.
%
%      GUI_SVDANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_SVDANALYSIS.M with the given input arguments.
%
%      GUI_SVDANALYSIS('Property','Value',...) creates a new GUI_SVDANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_SVDAnalysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_SVDAnalysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_SVDAnalysis

% Last Modified by GUIDE v2.5 15-Apr-2021 00:43:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_SVDAnalysis_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_SVDAnalysis_OutputFcn, ...
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


% --- Executes just before GUI_SVDAnalysis is made visible.
function GUI_SVDAnalysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_SVDAnalysis (see VARARGIN)

% Choose default command line output for GUI_SVDAnalysis
handles.output = hObject;
handles.CodePath=pwd;
addpath(handles.CodePath);
addpath([handles.CodePath, '\SubFunctions'])
handles.DefPath='F\';
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

% UIWAIT makes GUI_SVDAnalysis wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_SVDAnalysis_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Reset_pushbutton.
function Reset_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Reset_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.IQ = handles.IQ0;
handles.dt = 0;
guidata(hObject, handles);





function Input_LowRank_Callback(hObject, eventdata, handles)
% hObject    handle to Input_LowRank (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Input_LowRank as text
%        str2double(get(hObject,'String')) returns contents of Input_LowRank as a double


% --- Executes during object creation, after setting all properties.
function Input_LowRank_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Input_LowRank (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Input_HighRank_Callback(hObject, eventdata, handles)
% hObject    handle to Input_HighRank (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Input_HighRank as text
%        str2double(get(hObject,'String')) returns contents of Input_HighRank as a double


% --- Executes during object creation, after setting all properties.
function Input_HighRank_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Input_HighRank (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in LoadIQData_pushbutton.
function LoadIQData_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to LoadIQData_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,datapath]=uigetfile(handles.DefPath);
handles.DefPath=datapath;
guidata(hObject, handles);
disp('Loading data...');
load ([datapath, filename]);
disp('Data was loaded...');
handles.IQ=IQ; 
% handles.IQ=mIQ; % for multiple IQ combined data;
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
handles.dt = 0;
guidata(hObject, handles);


% --- Executes on button press in IQ2dIQ_pushbutton.
function IQ2dIQ_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to IQ2dIQ_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
IQ = handles.IQ;
dt = str2num(get(handles.Input_dt,'string'));
dIQ = IQ(:,:,(dt+1):end) - IQ(:,:,1:(end-dt));
disp('dIQ calculation finished!');
if dt ~= 0
handles.dIQ = dIQ;
end

handles.dt = dt;
% [handles.nz, handles.nx, handles.nt]=size(handles.IQ);
guidata(hObject, handles);



% --- Executes on button press in PlotSVD_pushbutton.
function PlotSVD_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to PlotSVD_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.dt==0
dataRAW = handles.IQ;
else
    dataRAW = handles.dIQ;
end

fCC = handles.fCC;
% DefCutFreq = handles.DefCutFreq;
IQR = dataRAW(:,:,1:handles.nCC);
[nz,nx,nt]=size(IQR);
% rank=[RankLow:RankHigh];
S=reshape(IQR,[nz*nx,nt]);
S_COVt=(S'*S);
[V,D]=eig(S_COVt); % V is the right singular Vector of S/eigenvector; D is the eigenvalue/square of Singular value
for it=1:nt 
    Ddiag(it)=abs(sqrt(D(it,it)));
end
Ddiag=20*log10(Ddiag/max(Ddiag)); % singular value in db
[Ddesc, Idesc]=sort(Ddiag,'descend');

for it=1:nt
    Vdesc(:,it)=V(:,Idesc(it)); % Vdesc is the right singluar matrix in SVD (has little numerical error)
end

%% plot singular values
figure;
yyaxis left; plot(Ddesc);
DDesc = cumsum(sort(sqrt(diag(D)),'descend'));
yyaxis right; plot(20*log10(DDesc/max(DDesc)));
% hold on; xline(rank(1)); 
xlabel('iNumber of Singular Values'); ylabel('Singular Values (dB)');

%% plot specturm of temporal singular vectors
for i = 1:nt
%     Vdesc0(:,i) = Vdesc(:,i).*conj(mean(Vdesc(:,i)))./abs(mean(Vdesc(:,i))).^2;
Vdesc0(:,i) = (Vdesc(:,i)-min(Vdesc(:,i)))./(max(Vdesc(:,i))-min(Vdesc(:,i)));
end
FVdesc = fftshift(fft(Vdesc0,[],1),1);
fCoor = linspace(-fCC/2,fCC/2,nt);

figure; imagesc([1:nt],fCoor, 20*log10(abs(FVdesc)/max(max(abs(FVdesc)))));colormap(jet);axis square;
% hold on; line([rank(1),rank(1)],[fCoor(1),fCoor(end)],'Color','black','LineWidth',1);
% hold on; line([1,nt],[-DefCutFreq, -DefCutFreq],'Color','black','LineWidth',1);
% hold on; line([1,nt],[DefCutFreq, DefCutFreq],'Color','black','LineWidth',1);
xlabel('iNumber of Singular Vectors'); ylabel('Frequency [Hz]'); 
title('Spectrum of Temperal Singular Vectors');colorbar;

figure; imagesc(abs(Vdesc)); axis image; colorbar; title('Temperal Singular Vectors');
axes(handles.axes2);
iNumber=str2num(get(handles.Input_iNumber,'string'));
    plot(real(Vdesc(:,iNumber)),'LineWidth',1);title('Temperal Singular Vectors Real Part');
    legend(strsplit(num2str(iNumber)),'Location','eastoutside');
axes(handles.axes3);
    plot(imag(Vdesc(:,iNumber)),'LineWidth',1);title('Temperal Singular Vectors Imag Part');
    legend(strsplit(num2str(iNumber)),'Location','eastoutside');
    
%% plot spatial singular Vectors
prompt={'Plot Spatial Singular Vectors?(0: No; 1: Yes)'};
name='plot U & V';
defaultvalue={'0'};%num2str(handles.newSlect),num2str(handles.DataType)};
numinput=inputdlg(prompt,name, 1, defaultvalue);
handles.plotUSlect=str2num(numinput{1});
if handles.plotUSlect==1
    disp('Calculating SVD...');
    tic
    [UU,DD,VV]=svd(S);
    for ith = 1:handles.nt
    Ui(:,:,ith) = reshape(UU(:,ith),[nz,nx]);
    end
    toc
    axes(handles.axes1);
    iNumber=str2num(get(handles.Input_iNumber,'string'));
    imagesc((abs(Ui(:,:,iNumber))).^0.4);  colormap(gray);colorbar;%axis image;
    title([num2str(iNumber),' th Spatial Singular Vector']);
    handles.Ui = Ui;
end

handles.Vdesc = Vdesc;
handles.Ddesc = Ddesc;
handles.FVdesc = FVdesc;
handles.fCoor = fCoor;
handles.S = S;
guidata(hObject, handles);




% --- Executes on button press in sIQ2PDI_pushbutton.
function sIQ2PDI_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to sIQ2PDI_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Vdesc = handles.Vdesc;
S = handles.S;
if handles.dt ==0
[nz,nx,nt] = size(handles.IQ);
else
    [nz,nx,nt] = size(handles.dIQ);
end
nt = handles.nCC;
LowRank = str2num(get(handles.Input_LowRank,'string'));
HighRank = str2num(get(handles.Input_HighRank,'string'));
rank = [LowRank: HighRank];

Vrank=zeros(size(Vdesc));
Vrank(:,rank)=Vdesc(:,rank);
Vnoise=zeros(size(Vdesc));
Vnoise(:,end-0.25*nt:end)=Vdesc(:,end-0.25*nt:end);%end-50:
UDelta=S*Vdesc;
sBlood0=reshape(UDelta*Vrank',[nz,nx,nt]);
% sBlood=sBlood0./repmat(std(abs(sBlood0),1,2),[1,nx]);
%%%% Noise equalization 
sNoise=reshape(UDelta*Vnoise',[nz,nx,nt]);
B=ones(100, 100);%[50,50]
% sNoiseMed=convn(abs(squeeze(mean(sNoise,3))),B,'same');
sNoiseMed=medfilt2(abs(squeeze(mean(sNoise,3))),[30,30],'symmetric');%[150,150]
sNoiseMedNorm=sNoiseMed/min(sNoiseMed(:));

prompt={'Import eqNoise?(0: No; 1: Yes)'};
name='Noise equalization';
defaultvalue={'0'};%num2str(handles.newSlect),num2str(handles.DataType)};
numinput=inputdlg(prompt,name, 1, defaultvalue);
handles.eqNoiseSlect=str2num(numinput{1});
if handles.eqNoiseSlect == 0
    handles.eqNoise = sNoiseMedNorm;
else
    [filename,datapath]=uigetfile(handles.DefPath);
    disp('Loading eqNoise...');
    load ([datapath, filename]);
    disp('eqNoise was loaded...');
    handles.eqNoise = mean(eqNoise,3);
end
guidata(hObject, handles);

%% high pass filter
[B,A]=butter(4,30/(handles.fCC/2),'high');    %coefficients for the high pass filter
sIQ1(:,:,101:100+nt)=sBlood0;
sIQ1(:,:,1:100)=flip(sIQ1(:,:,101:200),3);
sIQ2=filter(B,A,sIQ1,[],3);    % blood signal (filtering in the time dimension)
sIQHP=sIQ2(:,:,101:end); % High pass filtered sIQ
clear sIQ1 sIQ2

%% noise equalization
sBlood=sBlood0./repmat(handles.eqNoise,[1,1,nt]);
sIQHP = sIQHP./repmat(handles.eqNoise,[1,1,nt]);

PDISVD=mean(abs(sBlood).^2,3); 
PDISVDdb=10*log10(PDISVD./max(PDISVD(:))); % SVD-based PD image in dB

PDISVDHP=mean(abs(sIQHP).^2,3); 
PDISVDHPdb=10*log10(PDISVDHP./max(PDISVDHP(:))); % SVD-based PD image in dB

axes(handles.axes1);
imagesc(PDISVDHP.^0.35);axis image;colormap(hot);
title(['PDI',' Rank ',num2str(num2str(LowRank)),'-', num2str(num2str(HighRank))]);
% 
% %%%%%% mIQ %%%%%%%%%
% sIQ = reshape(sIQHP,[nz,nx,200,32]);
% sIQ = sIQ./repmat(sNoiseMedNorm,[1,1,200,32]);
% PDI=squeeze(mean(abs(sIQ).^2,3));
% figure; imagesc(PDI(:,:,1).^0.35);axis image;colormap(hot);
% title(['PDI',' Rank ',num2str(num2str(LowRank)),'-', num2str(num2str(HighRank))]);
% figure; imagesc(PDI(:,:,2).^0.35);axis image;colormap(hot);
% title(['PDI',' Rank ',num2str(num2str(LowRank)),'-', num2str(num2str(HighRank))]);
% figure; imagesc(PDI(:,:,3).^0.35);axis image;colormap(hot);
% title(['PDI',' Rank ',num2str(num2str(LowRank)),'-', num2str(num2str(HighRank))]);
% figure; imagesc(PDI(:,:,4).^0.35);axis image;colormap(hot);
% title(['PDI',' Rank ',num2str(num2str(LowRank)),'-', num2str(num2str(HighRank))]);
% figure; imagesc(PDI(:,:,5).^0.35);axis image;colormap(hot);
% title(['PDI',' Rank ',num2str(num2str(LowRank)),'-', num2str(num2str(HighRank))]);
% figure; plot(squeeze(mean(mean(PDI,1),2)));
% 
% hrf = hemodynamicResponse(1,[2 16 0.5 1 20 0]);
% stim = zeros(32,1);
% stim(1:2,:)=1;
% stimhrf = filter(hrf,1,stim);
% 
% coefmap = CoorCoeffMap(PDI, stimhrf', 0); 
% 
% %%%%%


if handles.dt==0
handles.sBlood0 = sBlood0;
handles.sBlood = sBlood;
handles.sIQHP = sIQHP;
else 
    handles.dsBlood0 = sBlood0;
handles.dsBlood = sBlood;
handles.dsIQHP = sIQHP;
end
handles.rank = rank;
handles.PDISVD = PDISVD;
handles.PDISVDdb = PDISVDdb;
handles.PDISVDHP = PDISVDHP;
handles.PDISVDHPdb = PDISVDHPdb;
guidata(hObject, handles);

% --- Executes on button press in PixSelectPlot_pushbotton.
function PixSelectPlot_pushbotton_Callback(hObject, eventdata, handles)
% hObject    handle to PixSelectPlot_pushbotton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
IQ = handles.IQ;
sIQ = handles.sBlood;
sIQHP = handles.sIQHP;
% if handles.dt ~= 0
dIQ = handles.dIQ;
dIQ = cat(3,dIQ(:,:,1),dIQ);
dsIQ = handles.dsBlood;
dsIQHP = handles.dsIQHP;
dsIQ = cat(3,dsIQ(:,:,1),dsIQ);
dsIQHP = cat(3,dsIQHP(:,:,1),dsIQHP);
% end
fCC = handles.fCC;
fCoor = linspace(-fCC/2,fCC/2,size(IQ,3));
axes(handles.axes1);   
imagesc(handles.PDISVDdb); axis image;colormap(gray); 
    [slt(2), slt(1)]=ginput(1); % [x z]
    slt=round(slt);

    IQAVG = squeeze(IQ(slt(1),slt(2),:));
    fIQ = fftshift(fft(IQAVG));

    dIQAVG=squeeze(dIQ(slt(1),slt(2),:));
    nt=length(IQAVG);
    fdIQ=fftshift(fft(dIQAVG));
    
    sIQAVG=squeeze(sIQ(slt(1),slt(2),:));
    nt=length(sIQAVG);
    FsIQ=fftshift(fft(sIQAVG));
    
    dsIQAVG=squeeze(dsIQ(slt(1),slt(2),:));
    nt=length(dsIQAVG);
    FdsIQ=fftshift(fft(dsIQAVG));
    
    Pneg=sum(abs(FsIQ(1:floor(end/2)))-median(abs(FsIQ)));
    Ppos=sum(abs(FsIQ(floor(end/2)+1:end))-median(abs(FsIQ)));
    rNeg=Pneg/(Pneg+Ppos);
    rPos=Ppos/(Pneg+Ppos);
    tCoor=linspace(1,nt,nt)/fCC*1000;
    
    figure;
    figure; subplot(221);
    plot(tCoor,abs(IQAVG),'k-');
    subplot(223);  plot(tCoor, abs(dIQAVG),'b-');
    subplot(222); plot(tCoor,abs(sIQAVG),'r-');
    subplot(224);  plot(tCoor,abs(dsIQAVG),'g-');
    
    figure;
    subplot(211);
    plot(tCoor,abs(IQAVG),'k-');
    hold on; yyaxis right; plot(tCoor, abs(dIQAVG),'b-');
    hold on;  plot(tCoor,abs(sIQAVG),'r-');
    hold on;  plot(tCoor,abs(dsIQAVG),'g-');
        xlabel('Time [ms]')
    ylabel('Pressure amplitude')
    legend({'IQ','dIQ','sIQ','dsIQ'})
    title(['Time course @ ROI ', '[',num2str(slt(1)),',',num2str(slt(2)),']']);
 
    
    subplot(2,1,2),plot(fCoor,abs(fIQ),'k')
    hold on; plot(fCoor, abs(fdIQ),'b')
    hold on; plot(fCoor,abs(FsIQ),'r')
    hold on; plot(fCoor,abs(FdsIQ),'g')
%     xlim([-800 800])
    ylim([0 mean(abs(fIQ))])
    xlabel('Frequency [Hz]')
    ylabel('Power Spectrum Density')
    legend({'IQ','dIQ','sIQ','dsIQ'})
    title(['Pneg=',num2str(rNeg),'; Ppos=',num2str(rPos)])
    
    
    figure;
    subplot(141); imagesc(fCoor,[1:handles.nz],abs(squeeze(fftshift(fft(IQ(:,slt(2),:))))));title('fIQ');
    subplot(142); imagesc(fCoor,[1:handles.nz],abs(squeeze(fftshift(fft(dIQ(:,slt(2),:))))));title('fdIQ');
    subplot(143); imagesc(fCoor,[1:handles.nz],abs(squeeze(fftshift(fft(sIQ(:,slt(2),:))))));title('fsIQ');
    subplot(144); imagesc(fCoor,[1:handles.nz],abs(squeeze(fftshift(fft(dsIQ(:,slt(2),:))))));title('fdsIQ')
   
    
        figure;
    subplot(141); imagesc(abs(squeeze(IQ(:,slt(2),:))));title('IQ');
    subplot(142); imagesc(abs(squeeze(dIQ(:,slt(2),:))));title('dIQ');
    subplot(143); imagesc(abs(squeeze(sIQ(:,slt(2),:))));title('sIQ');
    subplot(144); imagesc(abs(squeeze(dsIQ(:,slt(2),:))));title('dsIQ');
  
    


% --- Executes on slider movement.
function iNumber_slider_Callback(hObject, eventdata, handles)
% hObject    handle to iNumber_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
set(hObject,'SliderStep',[1/(handles.nt-1), 10/(handles.nt-1)])
set(hObject,'Min',1)
set(hObject,'Max',handles.nt)
iNumber=round(get(hObject,'Value'));
set(handles.Input_iNumber,'string',iNumber);

%% plot temperal singular vectors
Vdesc = handles.Vdesc;
axes(handles.axes2);
    plot(real(Vdesc(:,iNumber)),'LineWidth',1);title('Temperal Singular Vectors Real Part');
    legend(strsplit(num2str(iNumber)),'Location','eastoutside');
axes(handles.axes3);
    plot(imag(Vdesc(:,iNumber)),'LineWidth',1);title('Temperal Singular Vectors Imag Part');
    legend(strsplit(num2str(iNumber)),'Location','eastoutside');
    
%% plot spatial singular Vectors
if handles.plotUSlect==1
    Ui = handles.Ui;
    axes(handles.axes1);
    imagesc((abs(Ui(:,:,iNumber))).^0.4); axis image; colormap(gray);colorbar;
    title([num2str(iNumber),' th Spatial Singular Vector']); 
end


% --- Executes during object creation, after setting all properties.
function iNumber_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iNumber_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function Input_iNumber_Callback(hObject, eventdata, handles)
% hObject    handle to Input_iNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iNumber=str2num(get(handles.Input_iNumber,'string'));
%% plot temperal singular vectors
Vdesc = handles.Vdesc;
axes(handles.axes2);
iNumber=str2num(get(handles.Input_iNumber,'string'));
    plot(real(Vdesc(:,iNumber)),'LineWidth',1);title('Temperal Singular Vectors Real Part');
    legend(strsplit(num2str(iNumber)),'Location','eastoutside');
axes(handles.axes3);
    plot(imag(Vdesc(:,iNumber)),'LineWidth',1);title('Temperal Singular Vectors Imag Part');
    legend(strsplit(num2str(iNumber)),'Location','eastoutside');
    
%% plot spatial singular Vectors
if handles.plotUSlect==1
    Ui = handles.Ui;
    axes(handles.axes1);
    iNumber=str2num(get(handles.Input_iNumber,'string'));
    imagesc((abs(Ui(:,:,iNumber))).^0.4); axis image; colormap(gray);colorbar;
    title([num2str(iNumber),' th Spatial Singular Vector']); 
end


% Hints: get(hObject,'String') returns contents of Input_iNumber as text
%        str2double(get(hObject,'String')) returns contents of Input_iNumber as a double


% --- Executes during object creation, after setting all properties.
function Input_iNumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Input_iNumber (see GCBO)
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


% --------------------------------------------------------------------
function CopyImage_Callback(hObject, eventdata, handles)
% hObject    handle to CopyImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

f1 = figure % Open a new figure with handle f1
s = copyobj(handles.axes1,f1) % Copy axes object h into figure f1
set(gcf,'position',[100,100,800,500]);colorbar; colormap(hot);



function Input_dt_Callback(hObject, eventdata, handles)
% hObject    handle to Input_dt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Input_dt as text
%        str2double(get(hObject,'String')) returns contents of Input_dt as a double


% --- Executes during object creation, after setting all properties.
function Input_dt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Input_dt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
