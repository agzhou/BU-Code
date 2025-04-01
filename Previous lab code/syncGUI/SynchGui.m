% USE FOR DMITRY SPECKLE RECORDER
% SWITCH BRANCH FROM IOSI TRIGGER OUT TO FS AI6 / ACH6 IN
% MANUALLY START SPECKLE RECORDER AFTER STARTING ACQUISITION
% SPECKLE RECORDER LINE 4 OUTPUT EXPOSURE ACTIVE
function varargout = SynchGui(varargin)
%SYNCHGUI MATLAB code file for SynchGui.fig
%      SYNCHGUI, by itself, creates a new SYNCHGUI or raises the existing
%      singleton*.
%
%      H = SYNCHGUI returns the handle to a new SYNCHGUI or the handle to
%      the existing singleton*.
%
%      SYNCHGUI('Property','Value',...) creates a new SYNCHGUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to SynchGui_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      SYNCHGUI('CALLBACK') and SYNCHGUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in SYNCHGUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's
%      Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SynchGui

% Last Modified by GUIDE v2.5 23-Feb-2018 13:29:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SynchGui_OpeningFcn, ...
                   'gui_OutputFcn',  @SynchGui_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before SynchGui is made visible.
function SynchGui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for SynchGui
handles.output = hObject;
daqreset;
global Data
Data.frameReadOut = [];
Data.ts = [];
Data.stimulusTrigger = [];
Data.ledTrigger = [];
Data.speaker = [];
Data.lever = [];
Data.baslerExposure = [];

handles = setfield(handles,'delay_time',5000);
handles = setfield(handles,'stim_freq',3.0);
handles = setfield(handles,'stim_width',1);
handles = setfield(handles,'stim_length',2);
handles = setfield(handles,'seq_length',30);
handles = setfield(handles,'num_trials',5);
handles = setfield(handles,'init_delay',5);
handles = setfield(handles,'outputData',[]);
handles = setfield(handles,'synchData',Data);
handles = setfield(handles,'output_file','synchfile.mat');
handles = setfield(handles,'daq_freq',1000); %sample rate is 10 kHz for output

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SynchGui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SynchGui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in update_button.
function update_button_Callback(hObject, eventdata, handles)
% hObject    handle to update_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% params.DAQFrequency = 1000;
% params.PreStimulusBaseline = 0;
% params.PostStimulusBaseline = 0;
% params.GENUSduration = 10;
% params.StimulusFrequency = 40;
% params.ModFrequency = 0.2;



% trigger = f_GENUSvis(params)*5;
% 
% figure
% plot(trigger)
% outputData = [trigger, trigger];

outputData = generateStimulus(handles.delay_time,handles.stim_freq,...
    handles.stim_width,handles.stim_length,handles.seq_length);
handles.outputData = outputData;
Total_Number_Points = size(outputData);
Output_DAQ_Frequency_Hz = 10000;
 t_OutputChannels = (1:Total_Number_Points)./Output_DAQ_Frequency_Hz;
    plot(t_OutputChannels,outputData(:,2),'.-g',t_OutputChannels,outputData(:,1),'.-b')
    hold on
     %plot(t_OutputChannels,zeros(size(t_OutputChannels,1)),'.-b')
    xlabel('Time (s)');
    ylabel('Voltage (V)');
    title('One stimulation sequence');
    legend('CCD trigger','Stimulation pattern');
    hold off
guidata(hObject,handles)



function Delay_Callback(hObject, eventdata, handles)
% hObject    handle to Delay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Delay as text
%        str2double(get(hObject,'String')) returns contents of Delay as a double
handles.delay_time = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function Delay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Delay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on Delay and none of its controls.
function Delay_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to Delay (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)



function stimfreq_Callback(hObject, eventdata, handles)
% hObject    handle to stimfreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stimfreq as text
%        str2double(get(hObject,'String')) returns contents of stimfreq as a double
handles.stim_freq = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function stimfreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stimfreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in start_button.
function start_button_Callback(hObject, eventdata, handles)
% hObject    handle to start_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Data
Data.frameReadOut = [];
Data.ts = [];
Data.stimulusTrigger = [];
Data.ledTrigger = [];
Data.speaker = [];
Data.lever = [];
Data.baslerExposure = [];
device = daq.getDevices; %Evren
s = daq.createSession('ni') ; %Evren
s.Rate = 1000;

ch1 = addAnalogOutputChannel(s,'Dev1','ao0','Voltage');%For Stimulation
ch2 = addAnalogOutputChannel(s,'Dev1','ao1','Voltage');%For CCD
ch3 = addAnalogInputChannel(s,'Dev1','ai1','Voltage'); %From Hamamatsu camera
ch4 = addAnalogInputChannel(s,'Dev1','ai2','Voltage'); %From air puff stimulation picospitzer
ch5 = addAnalogInputChannel(s,'Dev1','ai3','Voltage'); %From led to show led trigger time
ch6 = addAnalogInputChannel(s,'Dev1','ai4','Voltage'); %From Speaker to show signal time
ch7 = addAnalogInputChannel(s,'Dev1','ai5','Voltage'); %From lever to show lever press time
ch8 = addAnalogInputChannel(s,'Dev1','ai6','Voltage'); %From Basler output trigger
lh = addlistener(s,'DataAvailable', @(src,event) recordData(event.TimeStamps, event.Data));
%s.IsContinuous = true;
save(handles.output_file,'Data')
outputStim = [];
for i = 1:handles.num_trials
    outputStim = [outputStim;handles.outputData];
end
disp('Start')
tic
queueOutputData(s,outputStim);

trialTime = toc
tic
s.startForeground();
disp('End') 
trialTime = toc
clear ch1;
clear ch2;
clear ai;
clear ch3;
clear ch4;
clear ch5;
clear ch6;
clear ch7;
clear ch8;
delete(lh);
elapsedTime = toc

daqreset;
daqreset;
trials = organizeTrials(handles.daq_freq,handles.seq_length,handles.num_trials,Data);
try
    trials = findFrames(trials);
catch
end
save(handles.output_file,'trials')

% --- Executes on button press in restart_button.
function restart_button_Callback(hObject, eventdata, handles)
% hObject    handle to restart_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(exist('s') ~= 0)
    s.stop();
end



function stim_length_Callback(hObject, eventdata, handles)
% hObject    handle to stim_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stim_length as text
%        str2double(get(hObject,'String')) returns contents of stim_length as a double
handles.stim_length = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function stim_length_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stim_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function stim_width_Callback(hObject, eventdata, handles)
% hObject    handle to stim_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stim_width as text
%        str2double(get(hObject,'String')) returns contents of stim_width as a double
handles.stim_width = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function stim_width_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stim_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function seq_length_Callback(hObject, eventdata, handles)
% hObject    handle to seq_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of seq_length as text
%        str2double(get(hObject,'String')) returns contents of seq_length as a double

handles.seq_length = str2double(get(hObject,'String'));
guidata(hObject,handles);
% --- Executes during object creation, after setting all properties.
function seq_length_CreateFcn(hObject, eventdata, handles)
% hObject    handle to seq_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_seq_Callback(hObject, eventdata, handles)
% hObject    handle to num_seq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_seq as text
%        str2double(get(hObject,'String')) returns contents of num_seq as a double
handles.num_trials = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function num_seq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_seq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function init_delay_Callback(hObject, eventdata, handles)
% hObject    handle to init_delay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of init_delay as text
%        str2double(get(hObject,'String')) returns contents of init_delay as a double
handles.init_delay = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function init_delay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to init_delay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function output_file_Callback(hObject, eventdata, handles)
% hObject    handle to output_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of output_file as text
%        str2double(get(hObject,'String')) returns contents of output_file as a double
handles.output_file = get(hObject,'String');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function output_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to output_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
