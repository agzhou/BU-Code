% Generate a single square wave trigger according to the DAQ sampling rate
% [samples/second], the delay between the start of the signal and the rising 
% edge [s], the pulse width [s], and the total duration of the signal [s].

function [Mcr_triggerSignal] = generateSingleTriggerSignal(daqrate, delay_s, pulse_width_s, total_duration_s)
%     % testing
%     Mcr_daqrate = Mcr_d.Rate;
%     delay_s = 1;
%     pulse_width_s = 0.5;
%     total_duration_s = 5;
    

    voltageHigh = 5; % 5 V default for a high signal
    Mcr_triggerSignal = zeros(total_duration_s * daqrate, 1);
    Mcr_triggerSignal(round(delay_s * daqrate):round((delay_s + pulse_width_s) * daqrate)) = voltageHigh;

end