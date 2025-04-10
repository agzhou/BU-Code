function [apis, vts, daqrate, numTrials] = functionalParameterInputPrompt
%     functionalParameterPrompt = {'Stim delay time [ms]', 'Stim frequency [Hz]', 'Stim pulse width [ms]', 'Stim duration within each trial [s]', 'Stim trial duration [s]', 'Verasonics trigger delay [s]', 'Verasonics trigger delay total duration [s]', 'Verasonics trigger pulse width [ms]'};
    functionalParameterPrompt = {'Stim delay time [ms]', 'Stim frequency [Hz]', 'Stim pulse width [ms]', 'Stim duration within each trial [s]', 'Stim trial duration [s]', 'Verasonics trigger delay [s]', 'Verasonics trigger pulse width [ms]', 'DAQ sample rate [Hz]', 'Number of trials per run'};
%     functionalParameterDefaults = {'5000', '3', '100', '1', '20', '6', '500', '1000', '1'};
    functionalParameterDefaults = {'5000', '3', '100', '5', '30', '10', '500', '1000', '1'};
%     functionalParameterDefaults = {'5000', '3', '100', '1', '20', '10', '500', '1000', '1'};
    functionalParameterUserInput = inputdlg(functionalParameterPrompt, 'Input Functional Stimulus Parameters', 1, functionalParameterDefaults);
    
    apis.delay_time_ms = str2double(functionalParameterUserInput{1}); % Delay before the start of the stimulation [ms]
    apis.stim_freq_Hz = str2double(functionalParameterUserInput{2});     % Stimulation (square wave) frequency [pulses/sec]
    apis.stim_width_ms = str2double(functionalParameterUserInput{3});  % Width of each square wave [ms]
    apis.stim_length_s = str2double(functionalParameterUserInput{4});   % Duration of the stimulation [s]
    apis.seq_length_s = str2double(functionalParameterUserInput{5});   % Total duration of the trial [s]
        
    % Verasonics trigger signal (vts)
    vts.delay_s = str2double(functionalParameterUserInput{6});
%     vts.total_duration_s = str2double(functionalParameterUserInput{7});
%     vts.pulse_width_ms = str2double(functionalParameterUserInput{8});
    vts.pulse_width_ms = str2double(functionalParameterUserInput{7});
%     vts.total_duration_s = apis.seq_length_s; % same as the air puff total duration

    daqrate = str2double(functionalParameterUserInput{8});
    numTrials = str2double(functionalParameterUserInput{9});
end