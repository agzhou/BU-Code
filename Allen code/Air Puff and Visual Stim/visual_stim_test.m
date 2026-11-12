Mcr_d = daq('ni'); % Create the DAQ object
Mcr_d.Rate = 1000; % DAQ rate [samples per second]

stimInputCh = addoutput(Mcr_d, 'Dev1', 'ao0', 'Voltage');  % Trigger to the visual stim fiber

% 4. Generate a simple waveform (e.g., 1 second of a 5 Hz sine wave from 0 to 5V)
t = (0:1/Mcr_d.Rate:1-1/Mcr_d.Rate)';
signal = 2.5 * sin(2 * pi * 5 * t) + 2.5; 

% 5. Write data to the NI device
write(Mcr_d, signal);