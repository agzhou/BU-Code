
% Try to write my own air puff code
% Requires Data Acquisition Toolbox and the NI package

%% Set up the hardware and channels
d = daq('ni'); % Create the DAQ object
airPuffInputCh = addoutput(d, 'Dev1', 'ao0', 'Voltage');  % Trigger to the air puff input
airPuffOutputCh = addinput(d, 'Dev1', 'ai6', 'Voltage');  % Record the actual air puff actions
d.Rate = 1000; % DAQ rate, I guess this is samples per second

%% Create the air puff input signal
% Air puff input signal (apis)
apis.delay_time = 5000; % Delay before the start of the stimulation [ms]
apis.stim_freq = 3;     % Stimulation (square wave) frequency [pulses/sec]
apis.stim_width = 100;  % Width of each square wave [ms]
apis.stim_length = 5;   % Duration of the stimulation [s]
apis.seq_length = 20;   % Total duration of the trial [s]

% airPuffInputSignal = ones(1000, 1) .* 5;
% airPuffInputSignal = linspace(-1, 1, 2200)' .* 3.3;
airPuffInputSignal = generateStimulus(apis.delay_time, apis.stim_freq,...
    apis.stim_width, apis.stim_length, apis.seq_length);
airPuffInputSignal = airPuffInputSignal(:, 1); % don't need the camera triggers like they did
% numCycles = 3;
% t = linspace(0, 2*pi * numCycles)';
% airPuffInputSignal = square(1 .* t) .* 3.3;

%% Run the trial in the foreground

disp('Start')
readwrite(d, airPuffInputSignal);
disp('End')

%% Run the trial in the background
preload(d, airPuffInputSignal);
start(d);