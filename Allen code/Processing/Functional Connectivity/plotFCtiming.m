%% Description:
% - Retrieve the time tags from each superframe file of a fUS acquisition
% - For use with Functional Connectivity data (no triggers)

clearvars
%% Choose data path and load the parameters and timing data
datapath = uigetdir('J:\', 'Select the raw data path');
datapath = [datapath, '\'];

% load([datapath, ''])
load([datapath, 'params.mat'])
% load([datapath, 'startTimeTag'])
% load([datapath, 'daqStartTimetag.mat'])
% load([datapath, 'triggerData'])

parameterPrompt = {'Start file', 'End file'};
parameterDefaults = {'1', ''};
parameterUserInput = inputdlg(parameterPrompt, 'Input Parameters', 1, parameterDefaults);

startFile = str2double(parameterUserInput{1});
endFile = str2double(parameterUserInput{2});

clearvars parameterPrompt parameterDefaults parameterUserInput

RFfilenameStructure = ['RF-', num2str(round(P.maxAngle)), '-', num2str(P.na), '-', num2str(P.frameRate), '-', num2str(P.numFramesPerBuffer), '-1-'];

% Load the first timetag
load([datapath, RFfilenameStructure, num2str(1)], 'timetag');
startTimetag = timetag; % The timetag at the start of the Verasonics sequence (after the trigger starts it)
clearvars timetag

% Load the DAQ-read data (e.g., accelerometer output)
load([datapath, 'daqStartTimetag.mat']) % Timetag at the start of the DAQ run
load([datapath, 'daqData.mat']) % Timetag at the start of the DAQ run

%% Get the timestamp of the end of each superframe (relative to the acqStart)

sfTimeTags = zeros(endFile - startFile + 1, 1); % Superframe time tags relative to the start superframe timetag

for filenum = startFile:endFile
    load([datapath, RFfilenameStructure, num2str(filenum)], 'timetag');
    sfTimeTags(filenum) = seconds(timetag - startTimetag);
end
clearvars timetag

figure; plot(diff(sfTimeTags)); xlabel('Timetag pair #'); ylabel('Timetag difference [s]')
%% Move the superframe time tags back the width of each superframe (so the time corresponds to the start of each superframe acquisition)
% ** actually, this is not necessary for sfTimeTags since we are using the
% relative timing for those, but adjusting for the sfWidth would be
% necessary for the DAQ timing **

%%%%%%% check if this works for the uneven spacing!! %%%%%%%
sfWidth = P.numFramesPerBuffer / P.frameRate; % How long each superframe takes to acquire
% sfWidth_adj = round(sfWidth * P.daqrate); % Adjust for the DAQ sampling rate (so we can plot it)

% sfTimeTags_end = sfTimeTags;
% sfTimeTags_start = sfTimeTags - sfWidth;

%% Adjust the DAQ timestamps to match the start of the 1st superframe
daqTimeTags_daqstart = timeStamp; % DAQ time stamps relative to the DAQ start [s]
VSXVsDAQDelay = seconds(startTimetag - daqStartTimetag) - sfWidth; % Delay between the DAQ start and the (start of the) first superframe acquisition

daqTimeTags = daqTimeTags_daqstart - VSXVsDAQDelay; % DAQ time stamps relative to the first superframe start [s]

% figure; plot(daqTimeTags, inScanData); xlabel("Time [s]")

%% Manually input the behavioral camera start time...

% ------------- ............. -------------

%% ?

% sfTimeTagsDAQStart = sfTimeTags + seconds(startTimetag - daqStartTimetag); % Superframe timetags relative to the DAQ start
% sfTimeTagsDAQStart_adj = round(sfTimeTagsDAQStart * P.daqrate); % Round to the nearest time according to the DAQ rate (so we can plot it)
% % figure(tf)
% % 
% % 
% % maxValue = max(max(airPuffOutput), max(P.Mcr_fcp.vts.signal)); % max value across both timecourses to get a max value for shading
% % minValue = 0;
% % yshade = [maxValue, minValue, minValue, maxValue]; % Go top right and clockwise for the shading patch vertices
% % 
% % % Go through the superframe time tags and shade the corresponding regions
% % for tti = 1:length(sfTimeTagsDAQStart_adj) % time tag index
% %     tt = sfTimeTagsDAQStart_adj(tti); % time tag
% % %     xline(round(tt))
% %     xshade = [tt, tt, tt - sfWidth_adj, tt - sfWidth_adj];
% %     patch(xshade, yshade, 'g', 'FaceAlpha', .3) % Plot the shaded region
% % end
% % legend('Air puff output', 'Verasonics trigger', 'Superframe acquisition')

%% Add relevant variables to a structure 'TD' (Timing Data)
TD.sfStartTimetag = startTimetag;       % Timetag at the end of the first superframe acquired
TD.sfTimeTags = sfTimeTags;             % Timestamps of each superframe relative to the first superframe acquisition time
TD.sfWidth = sfWidth;                   % Width/duration of each superframe [s]
TD.inScanData = inScanData;             % Output data that the DAQ reads (e.g., accelerometer)
TD.daqStartTimetag = daqStartTimetag;   % Timetag at the start of the DAQ run
TD.VSXVsDAQDelay = VSXVsDAQDelay;       % Delay between the DAQ start and the (start of the) first superframe acquisition
TD.daqTimeTags_daqstart = daqTimeTags_daqstart; % DAQ time stamps relative to the DAQ start [s]
TD.daqTimeTags = daqTimeTags;           % DAQ time stamps relative to the first superframe start [s]
% TD.sfTimeTagsDAQStart = sfTimeTagsDAQStart;
% TD.sfTimeTagsDAQStart_adj = sfTimeTagsDAQStart_adj;

% TD.sfWidth_adj = sfWidth_adj;
% TD.timeStamp = timeStamp;

%% Save TD as a mat file
TDpath = uigetdir('J:\', 'Select the TD save path');
TDpath = [TDpath, '\'];

save([TDpath, 'TD.mat'], 'TD')