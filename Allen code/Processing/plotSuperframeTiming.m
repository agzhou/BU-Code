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

RFfilenameStructure = ['RF-', num2str(P.maxAngle), '-', num2str(P.na), '-', num2str(P.frameRate), '-', num2str(P.numFramesPerBuffer), '-1-'];

load([datapath, RFfilenameStructure, num2str(1)], 'timetag');
acqStart = timetag; % The timetag at the start of the Verasonics sequence (after the trigger starts it)


%% Get the timestamp of the end of each superframe (relative to the acqStart)
sfTimeTags = zeros(endFile - startFile + 1, 1); % Superframe time tags relative to the acqStart timetag

for filenum = startFile:endFile
    load([datapath, RFfilenameStructure, num2str(filenum)], 'timetag');
    sfTimeTags(filenum) = seconds(timetag - acqStart);
end
clearvars timetag

figure; plot(diff(sfTimeTags), '-o', 'LineWidth', 2); xlabel('s'); title('Time per superframe')