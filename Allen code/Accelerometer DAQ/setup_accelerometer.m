
% Description: Get accelerometer

% Requires Data Acquisition Toolbox and the NI package
% Connect the air puffer (PicoSpritzer III) to the NI DAQ
function [Mcr_d, Mcr_fcp] = setup_accelerometer()

    %% User input: parameters for the DAQ and its Verasonics start trigger
    parameterPrompt = {'Verasonics trigger delay [s]', 'Verasonics trigger pulse width [ms]', 'DAQ sample rate [Hz]', 'Experiment length [s]'};
    parameterDefaults = {'10', '500', '1000', '360'};
    parameterUserInput = inputdlg(parameterPrompt, 'Input DAQ Parameters', 1, parameterDefaults);
    
    % Verasonics trigger signal (vts)
    vts.delay_s = str2double(parameterUserInput{1}); % Delay between starting VSX and sending the start trigger [s]
    vts.pulse_width_ms = str2double(parameterUserInput{2}); % Pulse width of the Verasonics start trigger (square pulse) [ms]

    daqrate = str2double(parameterUserInput{3}); % DAQ sample rate [Hz]
    seq_length = str2double(parameterUserInput{4}); % Experiment/sequence duration [s]
    
    Mcr_fcp.daqrate = daqrate;

    %% Set up the hardware and channels
    Mcr_d = daq('ni'); % Create the DAQ object

    % User prompt for which DAQ channels to use
    channelPrompt = {'Accelerometer X channel (analog in)', 'Accelerometer Y channel (analog in)', 'Accelerometer Z channel (analog in)', 'Verasonics trigger channel (analog out)'};
    channelDefaults = {'3', '4', '5', '1'};
    channelUserInput = inputdlg(channelPrompt, 'Input DAQ channel numbers', 1, channelDefaults);
    
%     accelZOutputCh = addinput(Mcr_d, 'Dev1', 'ai5', 'Voltage');  % Record the Z acceleration
    channelIndTemp = 1;
    accelXOutputCh = addinput(Mcr_d, 'Dev1', "ai" + channelUserInput{channelIndTemp}, 'Voltage');  % Record the X acceleration
    
    channelIndTemp = channelIndTemp + 1;
    accelYOutputCh = addinput(Mcr_d, 'Dev1', "ai" + channelUserInput{channelIndTemp}, 'Voltage');  % Record the Y acceleration

    channelIndTemp = channelIndTemp + 1;
    accelZOutputCh = addinput(Mcr_d, 'Dev1', "ai" + channelUserInput{channelIndTemp}, 'Voltage');  % Record the Z acceleration
    
    channelIndTemp = channelIndTemp + 1;
%     verasonicsTriggerCh = addoutput(Mcr_d, 'Dev1', 'ao1', 'Voltage');  % Trigger to start the Verasonics acquisition
    verasonicsTriggerCh = addoutput(Mcr_d, 'Dev1', "ao" + channelUserInput{channelIndTemp}, 'Voltage');  % Trigger to start the Verasonics acquisition
%     verasonicsOutputCh = addinput(Mcr_d, 'Dev1', 'ai7', 'Voltage');  % Record the start of (triggers from) each acquired superframe/buffer

%     Mcr_d.Rate = 1000; % DAQ rate [samples per second]
    Mcr_d.Rate = daqrate; % DAQ rate [samples per second]
    
    %% Create the air puff input signal
    % Functional control parameters (fcp)

%     % Air puff input signal (apis)
%     Mcr_fcp.apis = apis;
%     
%     % Use the previous lab code to generate the signal
%     Mcr_fcp.apis.signal = generateStimulus_variablerate(Mcr_fcp.apis.delay_time_ms, Mcr_fcp.apis.stim_freq_Hz, ...
%         Mcr_fcp.apis.stim_width_ms, Mcr_fcp.apis.stim_length_s, Mcr_fcp.apis.seq_length_s, Mcr_d.Rate);
%     Mcr_fcp.apis.signal = Mcr_fcp.apis.signal(:, 1); % don't need the camera triggers like they did
%     
    % Verasonics trigger signal (vts)
    Mcr_fcp.vts = vts;
    Mcr_fcp.vts.pulse_width_s = vts.pulse_width_ms / 1000;
    Mcr_fcp.vts.signal = generateSingleTriggerSignal(Mcr_d.Rate, Mcr_fcp.vts.delay_s, Mcr_fcp.vts.pulse_width_s);
    
%     % Repeat the stimulus signal for each trial
%     Mcr_fcp.apis.signal = repmat(Mcr_fcp.apis.signal, numTrials, 1);

%     % Pad the stimulus signal with zeros while the Verasonics trigger delay
%     % is going
%     Mcr_fcp.apis.signal = [zeros(length(Mcr_fcp.vts.signal), 1); Mcr_fcp.apis.signal];

    % Pad the Verasonics trigger delay with zeros according to the length
    % of the experiment. (Also, in case the lengths of each signal are 
    % somehow different - they need to be the same to be concatenated 
    % later into a matrix for 'preload'
%     Mcr_fcp.vts.signal = padarray(Mcr_fcp.vts.signal, length(Mcr_fcp.apis.signal) - length(Mcr_fcp.vts.signal), 'post');
    Mcr_fcp.vts.signal = padarray(Mcr_fcp.vts.signal, seq_length * daqrate, 'post'); % Start counting the experiment length after the trigger pulse occurs. Note that my code adds 2x the trigger pulse width to the trigger signal originally to give some buffer room.

    %% Run the trial in the background
    % global data;
    stop(Mcr_d);
    flush(Mcr_d);
    
    % d.ScansAvailableFcn = @plotMyData; % add a callback function to use
    
%     preload(Mcr_d, [Mcr_fcp.apis.signal, Mcr_fcp.vts.signal]);
    preload(Mcr_d, [Mcr_fcp.vts.signal]);
    % numTrials = 1;
    % start(d, 'NumScans', numTrials);
    disp('==== Starting the experiment ====')
    start(Mcr_d)             % Finite-length scan
%     start(Mcr_d, 'continuous') % Continuous scan

    %% Read the output data (what the DAQ actually saw)
%     [inScanData, timeStamp, triggerTime] = read(d, seconds(apis.seq_length), "OutputFormat", "Matrix");
%     disp('End')
end