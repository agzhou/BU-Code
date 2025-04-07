
% Trying to write my own air puff code
% Requires Data Acquisition Toolbox and the NI package
% Connect the air puffer (PicoSpritzer III) to the NI DAQ
function [Mcr_d, Mcr_fcp] = controlAirPuff_func
    %% Set up the hardware and channels
    Mcr_d = daq('ni'); % Create the DAQ object
    airPuffInputCh = addoutput(Mcr_d, 'Dev1', 'ao0', 'Voltage');  % Trigger to the air puff input
    airPuffOutputCh = addinput(Mcr_d, 'Dev1', 'ai6', 'Voltage');  % Record the actual air puff actions
    verasonicsTriggerCh = addoutput(Mcr_d, 'Dev1', 'ao1', 'Voltage');  % Trigger to start the Verasonics acquisition

    Mcr_d.Rate = 1000; % DAQ rate [samples per second]
    
    %% Create the air puff input signal
    % Functional control parameters (fcp)
    % Air puff input signal (apis)
    Mcr_fcp.apis.delay_time = 10000; % Delay before the start of the stimulation [ms]
    Mcr_fcp.apis.stim_freq = 3;     % Stimulation (square wave) frequency [pulses/sec]
    Mcr_fcp.apis.stim_width = 100;  % Width of each square wave [ms]
    % apis.stim_length = 5;   % Duration of the stimulation [s]
    Mcr_fcp.apis.stim_length = 1;   % Duration of the stimulation [s]
    % apis.seq_length = 20;   % Total duration of the trial [s]
    Mcr_fcp.apis.seq_length = 20;   % Total duration of the trial [s]
    
    % airPuffInputSignal = ones(1000, 1) .* 5;
    % airPuffInputSignal = linspace(-1, 1, 2200)' .* 3.3;
    Mcr_fcp.apis.signal = generateStimulus(Mcr_fcp.apis.delay_time, Mcr_fcp.apis.stim_freq,...
        Mcr_fcp.apis.stim_width, Mcr_fcp.apis.stim_length, Mcr_fcp.apis.seq_length);
    Mcr_fcp.apis.signal = Mcr_fcp.apis.signal(:, 1); % don't need the camera triggers like they did
    
    % Verasonics trigger signal (vts)
    Mcr_fcp.vts.delay_s = 4;
    Mcr_fcp.vts.pulse_width_s = 0.5;
    Mcr_fcp.vts.total_duration_s = Mcr_fcp.apis.seq_length;
    Mcr_fcp.vts.signal = generateSingleTriggerSignal(Mcr_d.Rate, Mcr_fcp.vts.delay_s, Mcr_fcp.vts.pulse_width_s, Mcr_fcp.vts.total_duration_s);
    
    % in case the lengths of each signal are somehow different - they need
    % to be the same to be concatenated later into a matrix for 'preload'
    Mcr_fcp.vts.signal = padarray(Mcr_fcp.vts.signal, length(Mcr_fcp.apis.signal) - length(Mcr_fcp.vts.signal), 'post');

    % numCycles = 3;
    % t = linspace(0, 2*pi * numCycles)';
    % airPuffInputSignal = square(1 .* t) .* 3.3;
    
    %% Run the trial in the foreground
    
    % disp('Start')
    % readwrite(d, airPuffInputSignal);
    % disp('End')
    
    %% Run the trial in the background
    % global data;
    stop(Mcr_d);
    flush(Mcr_d);
    
    % d.ScansAvailableFcn = @plotMyData; % add a callback function to use
    
%     preload(Mcr_d, Mcr_fcp.apis.signal);
    preload(Mcr_d, [Mcr_fcp.apis.signal, Mcr_fcp.vts.signal]);
    numTrials = 1;
    % start(d, 'NumScans', numTrials);
    disp('Start')
    start(Mcr_d)
    
    %% Read the output data (what the air puff actually did)
%     [inScanData, timeStamp, triggerTime] = read(d, seconds(apis.seq_length), "OutputFormat", "Matrix");
%     disp('End')
end