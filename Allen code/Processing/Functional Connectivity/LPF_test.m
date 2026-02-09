%% Set up a Low Pass Filter
fc = 0.2; % Cutoff frequency [Hz]
fs = 1/mean(diff(t)); % Sampling frequency [Hz]
LPF.order = 3; % Butterworth filter order

[LPF.b, LPF.a] = butter(LPF.order, fc/(fs/2), 'low');

%% Apply low pass filter
    % dim = length(size(PDIallSF_reg)); % Operate on the time dimension
    % PDIallSF_reg_LPF = filter(LPF.b, LPF.a, PDIallSF_reg, [], dim);

dim = length(size(PDIallBlocks_reg)); % Operate on the time dimension
PDIallBlocks_reg_LPF = filter(LPF.b, LPF.a, PDIallBlocks_reg, [], dim);

%% Calculate GVTD after LPF
size_PDIallSF_reg = size(PDIallBlocks_reg);
numVoxelsInVolume = size_PDIallSF_reg(1) * size_PDIallSF_reg(2) * size_PDIallSF_reg(3); % # of voxels in each volume
GVTD_LPF = squeeze( sum( diff(PDIallBlocks_reg_LPF, 1, length(size(PDIallBlocks_reg_LPF))) .^ 2, [1, 2, 3] ) ./ numVoxelsInVolume ) .^ 0.5;
GVTD_LPF(end + 1) = NaN; % pad the end with a NaN, since there is no forward point past the last time point

%%
figure
yyaxis left
plot(GVTD)
yyaxis right
plot(GVTD_LPF)