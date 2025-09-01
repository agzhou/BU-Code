%% Description:
% - Retrieve the time tags from each superframe file of a fUS acquisition
% - Plot the various triggers with each superframe's duration overlaid

%% Choose data path and load the parameters and timing data
datapath = uigetdir('J:\Allen\Data\', 'Select the raw data path');
datapath = [datapath, '\'];

% load([datapath, ''])
load([datapath, 'params.mat'])
load([datapath, 'startTimeTag'])
load([datapath, 'daqStartTimetag.mat'])
load([datapath, 'triggerData'])

parameterPrompt = {'Start file', 'End file'};
parameterDefaults = {'1', ''};
parameterUserInput = inputdlg(parameterPrompt, 'Input Parameters', 1, parameterDefaults);

startFile = str2double(parameterUserInput{1});
endFile = str2double(parameterUserInput{2});

clearvars parameterPrompt parameterDefaults parameterUserInput

RFfilenameStructure = ['RF-', num2str(P.maxAngle), '-', num2str(P.na), '-', num2str(P.frameRate), '-', num2str(P.numFramesPerBuffer), '-1-'];
acqStart = timetag; % The timetag at the start of the Verasonics sequence (after the trigger starts it)

%% Get the air puff output and Verasonics trigger timecourses

airPuffOutput = inScanData(:, 1); % air puff output
% clearvars inScanData
tf = figure; plot(airPuffOutput)
hold on
plot(P.Mcr_fcp.vts.signal) % Verasonics trigger
% hold off
legend('Air puff output', 'Verasonics trigger')
xlabel("Samples [" + num2str(P.daqrate) + " Hz sampling rate]")

sfTimeTags = zeros(endFile - startFile + 1, 1); % Superframe time tags relative to the acqStart timetag

%% Get the timestamp of the end of each superframe (relative to the acqStart)

for filenum = startFile:endFile
    load([datapath, RFfilenameStructure, num2str(filenum)], 'timetag');
    sfTimeTags(filenum) = seconds(timetag - acqStart);
end
clearvars timetag

%% 
sfTimeTagsDAQStart = sfTimeTags + seconds(acqStart - daqStartTimetag); % Superframe timetags relative to the DAQ start
sfTimeTagsDAQStart_adj = round(sfTimeTagsDAQStart * P.daqrate); % Round to the nearest time according to the DAQ rate (so we can plot it)
figure(tf)

%%%%%%% check if this works for the uneven spacing!! %%%%%%%
sfWidth = P.numFramesPerBuffer / P.frameRate; % How long each superframe takes to acquire
sfWidth_adj = round(sfWidth * P.daqrate); % Adjust for the DAQ sampling rate (so we can plot it)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxValue = max(max(airPuffOutput), max(P.Mcr_fcp.vts.signal)); % max value across both timecourses to get a max value for shading
minValue = 0;
yshade = [maxValue, minValue, minValue, maxValue]; % Go top right and clockwise for the shading patch vertices

% Go through the superframe time tags and shade the corresponding regions
for tti = 1:length(sfTimeTagsDAQStart_adj) % time tag index
    tt = sfTimeTagsDAQStart_adj(tti); % time tag
%     xline(round(tt))
    xshade = [tt, tt, tt - sfWidth_adj, tt - sfWidth_adj];
    patch(xshade, yshade, 'g', 'FaceAlpha', .3) % Plot the shaded region
end
legend('Air puff output', 'Verasonics trigger', 'Superframe acquisition')

%% Add relevant variables to a structure 'TD' (Timing Data)
TD.acqStart = acqStart;
TD.airPuffOutput = airPuffOutput;
TD.daqStartTimetag = daqStartTimetag;
TD.sfTimeTags = sfTimeTags;
TD.sfTimeTagsDAQStart = sfTimeTagsDAQStart;
TD.sfTimeTagsDAQStart_adj = sfTimeTagsDAQStart_adj;
TD.sfWidth = sfWidth;
TD.sfWidth_adj = sfWidth_adj;
TD.timeStamp = timeStamp;