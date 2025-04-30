%% load params and stuff
IQpath = uigetdir('D:\Allen\Data\', 'Select the IQ data path');
IQpath = [IQpath, '\'];

% Load parameters
% if ~exist('P', 'var')
%     load([IQpath, '..\params.mat'])
% end
% Load acquisition parameters: params.mat
if ~exist('P', 'var')
    % Choose and load the params.mat file (from the acquisition)
    [params_filename, params_pathname, ~] = uigetfile('*.mat', 'Select the params file', [IQpath, '..\params.mat']);
    load([params_pathname, params_filename])
end

% Load Verasonics reconstruction parameters: datapath\PData.mat
if ~exist('PData', 'var')
    load([IQpath, 'PData.mat'])
end

IQfilenameStructure = ['IQ-', num2str(P.maxAngle), '-', num2str(P.na), '-', num2str(P.frameRate), '-', num2str(P.numFramesPerBuffer), '-1-'];

savepath = uigetdir('D:\Allen\Data\', 'Select the save path');
savepath = [savepath, '\'];

addpath([cd, '\Speckle tracking']) % add path for the g1 calculation functions

% Load the timing data
[timingFilePathFN, timingFilePath] = uigetfile([IQpath, '..\Timing data'], 'Select the timing data');
timingFilePath = [timingFilePath, timingFilePathFN];
load(timingFilePath)
% load(timingFilePath, 'acqStart', 'airPuffOutput', 'daqStartTimetag', 'sfTimeTags', 'sfTimeTagsDAQStart', 'sfTimeTagsDAQStart_adj', 'sfWidth', 'sfWidth_adj', 'timeStamp')
%% Define some parameters (add this to a prompt later)
sv_threshold_lower = 10;
sv_threshold_upper = 150;

startFile = 1;
endFile = 148;

taustep = 1/P.frameRate;
% tau = taustep:taustep:(P.numFramesPerBuffer * taustep);
tau = 0:taustep:((P.numFramesPerBuffer - 1) * taustep);
tau_ms = tau .* 1000; % Assuming even time spacing between frames

tau1_index_CBF = 2;
tau2_index_CBF = 6;
tau1_index_CBV = 2;

%% Main loop
for filenum = startFile:endFile
% for filenum = [37, 110, 111, 123:endFile]
% for filenum = 7
    tic
    load([IQpath, IQfilenameStructure, num2str(filenum)])
    
    IQ = squeeze(IData + 1i .* QData);
    clearvars IData QData
    
    % SVD decluttering
    [xp, yp, zp, nf] = size(IQ);
    
    [PP, EVs, V_sort] = getSVs2D(IQ);
    disp('SVs decomposed')
    [IQf] = applySVs2D(IQ, PP, EVs, V_sort, sv_threshold_lower, sv_threshold_upper);
    disp('SVD filtered images put together')

    clearvars IQ

    % Use the IQf with separated negative and positive frequency components
    [IQf_separated, IQf_FT_separated] = separatePosNegFreqs(IQf);
    
%     g1_n = g1T(IQf_separated{1}, 10);
%     [CBFi_n, CBVi_n] = g1_to_CBi(g1_n, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
%     g1_p = g1T(IQf_separated{2}, 10);
%     [CBFi_p, CBVi_p] = g1_to_CBi(g1_p, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
% 
%     % g1 = g1T(IQf, 10); % Only get the first 10 points
%     g1 = g1T(IQf);
%     [CBFi, CBVi] = g1_to_CBi(g1, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
% 
% %     savefast([savepath, 'fUSdata-', num2str(filenum), '.mat'], g1, CBFi, CBVi);
%     save([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'g1', 'CBFi', 'CBVi', '-v7.3', '-nocompression');
%     disp("fUS result for file " + num2str(filenum) + " saved" )

    [PDI] = calcPowerDoppler(IQf_separated);
    [CDI] = calcColorDoppler(IQf_FT_separated, P);

    save([savepath, 'PDI_CDI-', num2str(filenum), '.mat'], 'PDI', 'CDI', '-v7.3', '-nocompression');
    disp("PDI and CDI for file " + num2str(filenum) + " saved" )

    toc
    
end
% savefast([savepath, 'fUS_proc_params.mat'], 'sv_threshold_lower', 'sv_threshold_upper', 'tau', 'tau_ms', 'tau1_index_CBF', 'tau2_index_CBF', 'tau1_index_CBV');
savefast([savepath, 'PDI_CDI_proc_params.mat'], 'sv_threshold_lower', 'sv_threshold_upper');

%% g1 adjustment test (see MAIN_g1fUS_invivo_annotated.m)
g1_shift = g1(:, :, :, 2:end);
g1Temp = reshape(g1_shift, [size(g1_shift, 1) * size(g1_shift, 2) * size(g1_shift, 3), size(g1_shift, 4)]); % ** stack the spatial dimensions of the g1 **
                
% ** This seems like noise removal (CR = Clutter Removal?) **

% ** The difference between the first 2 g1 values is greater than 2x the
%    difference between the 2nd and 3rd g1 values **
% ggCR1=(abs( g1Temp(:,1) - g1Temp(:,2) ) > 2*abs(( g1Temp(:,2) - g1Temp(:,3) )));

% ** This is a series of masks that only keep the g1 at voxels where the
%    1st, 2nd, and 3rd g1 values meet some threshold criteria:
%    - |g1| at tau1 > 0.55 (I assume this means just if there is any blood
%                           signal at this voxel)
%    - |g1| at tau2 < 0.25 (If the g1 decays quickly enough)
%    - |g1| at tau2 < g1 at tau3 **
% ggCR2=( abs(g1Temp(:,1)) > 0.55) .* ( abs(g1Temp(:,2)) < 0.25 ) .* ( abs(g1Temp(:,2)) < abs(g1Temp(:,3)) );

% ********* trying my own thresholds *********
ggCR2 = ( abs(g1Temp(:,1)) > 0.5) .* ( abs(g1Temp(:,2)) < abs(g1Temp(:,1)) );
ggCR1 = (abs( abs(g1Temp(:,1)) - abs(g1Temp(:,2)) ) > 2 * abs(( abs(g1Temp(:,2)) - abs(g1Temp(:,3)) )));

ggCR0 = ggCR2;
% ggCR0=((ggCR1+ggCR2)>0); % modified by Bingxue Liu; ** logical 'OR' operation between the two conditions **
%GG2(:,1)=(1-ggCR0).*GG2(:,1)+ggCR0.*(GG2(:,2)+(abs(real(GG2(:,2)-GG2(:,3)))+1i*abs(imag(GG2(:,2)-GG2(:,3))))*1.5);
GG2temp(:,1)=(1-ggCR0).*g1Temp(:,1)+ggCR0.*(g1Temp(:,2)+(abs(real(g1Temp(:,1)-g1Temp(:,2)))+1i*(imag(g1Temp(:,2)-g1Temp(:,3))))*1);% abs; real GG1-GG2
GG2temp(find(abs(GG2temp)>1)) = g1Temp(find(abs(GG2temp)>1),1); %
g1Temp(:,1) = GG2temp; % modified by Bingxue Liu
g1Adj = reshape(g1Temp, [size(g1_shift, 1), size(g1_shift, 2), size(g1_shift, 3), size(g1_shift, 4)]); % ** unstack the spatial dimensions **
      
%% Plot CBVi and CBFspeedi from the smoothed g1
tau1_index_CBF = 2;
tau2_index_CBF = 6;
tau1_index_CBV = 2;

[CBFi_adj, CBVi_adj] = g1_to_CBi(g1Adj, tau_ms(2:end) ./ 1000, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
%
plotMIPs(CBFi_adj, 1)
plotMIPs(CBVi_adj, 1)
plotMIPs(CBVi, 1)

%% Look at the g1 adjustment criteria
ggCR1_rs = reshape(ggCR1, [size(g1_shift, 1), size(g1_shift, 2), size(g1_shift, 3)]); % ** unstack the spatial dimensions **
figure; imagesc(squeeze(max(ggCR1_rs, [], 1))')

ggCR2_rs = reshape(ggCR2, [size(g1_shift, 1), size(g1_shift, 2), size(g1_shift, 3)]); % ** unstack the spatial dimensions **
figure; imagesc(squeeze(max(ggCR2_rs, [], 1))') %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% More smoothing from the previous code
covB=ones(3,3); covB(3,3)=9; covB=covB/sum(covB(:)); % no clue wtf this is

g1Smoothed = g1;
for itau=1:size(g1, 4)
    g1Smoothed(:, :, :, itau) = convn(g1(:, :, :, itau), covB, 'same');
end
g1Smoothed = smoothdata(g1Smoothed, 4, 'sgolay', 9);

%% Plot CBVi and CBFspeedi from the smoothed g1
tau1_index_CBF = 2;
tau2_index_CBF = 6;
tau1_index_CBV = 2;

[CBFi_smoothed, CBVi_smoothed] = g1_to_CBi(g1Smoothed, tau_ms ./ 1000, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
%
% plotMIPs(CBFi_smoothed, 1)
plotMIPs(CBVi_smoothed, 1)
plotMIPs(CBVi, 1)

%% Phase testing
% plotMIPs(squeeze(IData(:, :, :, 1, 1)), 1)
% plotMIPs(squeeze(QData(:, :, :, 1, 1)), 1)
phasetest = squeeze(QData(:, :, :, 1, 1) ./ IData(:, :, :, 1, 1));
phasetest = phasetest(abs(phasetest) < 10);
plotMIPs(phasetest(:, :, :, 1), 1)
findfigs
%% Power Doppler
[PDI] = calcPowerDoppler(IQf_separated);
% plotMIPs(PDI{1}, 0.8)
% plotMIPs(PDI{2}, 0.8)
plotMIPs(PDI{3}, 0.8)

% volumeViewer(PDI{3})
% volumeSegmenter(PDI{3})
%% Plot Power Doppler results
PDIn = PDI{3}; % normalized PDI of all frequencies
PDIn = PDIn ./ max(PDIn, [], 'all');
plotMIPs(PDIn, 1)

%% trying to denoise the Power Doppler
% PDInt = PDIn; % normalized, thresholded
% PDInt(PDInt < 0.07) = 0;
% plotMIPs(PDInt, 1)
% %%
% test = NLMF(PDIn);
% %%
% plotMIPs(test, 1)
%% Color Doppler
[CDI] = calcColorDoppler(IQf_FT_separated, P);
%% Plot Color Doppler
plotMIPs(CDI{1}, 1)
plotMIPs(CDI{2}, 1)
plotMIPs(CDI{3}, 1)

% volumeViewer(CDI{})
% volumeSegmenter(CDI{1})

%% CBVi and CBFi MIP over the whole dimension with negative and positive components
% plotMIPs(CBVi_n, 1)
% plotMIPs(CBFi_n, 1)

% plotMIPs(CBVi_p, 1)
% plotMIPs(CBFi_p, 1)

plotMIPs(CBVi, 1)
plotMIPs(CBFi, 1)

%% Plot the magnitude of g1 at some point, of the adjusted g1
figure; plot(tau_ms(2:size(g1, 4)), abs(squeeze(g1_shift(40, 45, 61, :))), '-o')
title('|g1| at 40, 45, 61')
xlabel('Tau [ms]')
ylabel('|g1|')

% hold on
figure
plot(tau_ms(2:size(g1, 4)), abs(squeeze(g1Adj(40, 45, 61, :))), '-o')
% hold off
% legend('Original g1', 'Adjusted g1')
%% Plot the base vs. smoothed g1
figure; plot(tau_ms(1:size(g1, 4)), abs(squeeze(g1(40, 45, 61, :))), '-o')
title('|g1| at 40, 45, 61')
xlabel('Tau [ms]')
ylabel('|g1|')

figure
plot(tau_ms(1:size(g1Smoothed, 4)), abs(squeeze(g1Smoothed(40, 45, 61, :))), '-o')

%%
figure; plot(tau_ms(2:size(g1, 4)), abs(squeeze(g1_shift(30, 60, 71, :))), '-o')
title('|g1| at 30, 60, 71')
xlabel('Tau [ms]')
ylabel('|g1|')

% hold on
figure
plot(tau_ms(2:size(g1, 4)), abs(squeeze(g1Adj(30, 60, 71, :))), '-o')
% hold off
% legend('Original g1', 'Adjusted g1')
%%
% [CBFi, CBVi] = g1_to_CBi(g1, tau_ms, 2, 3, 2); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
%%
figure; imagesc(squeeze(CBF(40, :, :))')
title('CBFi - xz plane')
xlabel('x pixels')
ylabel('z pixels')
figure; imagesc(squeeze(CBF(:, 40, :))')
title('CBFi - yz plane')
xlabel('x pixels')
ylabel('z pixels')

%% CBF MIP over the whole dimension
figure; imagesc(squeeze(max(CBFi_smoothed, [], 1))' .^ 1); colormap hot
title('CBFi - xz MIP')
xlabel('x pixels')
ylabel('z pixels')
figure; imagesc(squeeze(max(CBFi, [], 2))' .^ 1); colormap hot
title('CBFi - yz MIP')
xlabel('x pixels')
ylabel('z pixels')

% figure; imagesc(squeeze(sum(CBF, 1))' .^ 1); colormap hot
% title('CBFi - sum over y')
% xlabel('x pixels')
% ylabel('z pixels')
%%
figure; imagesc(squeeze(CBV(40, :, :))'); colormap hot
title('CBVi - xz plane')
xlabel('x pixels')
ylabel('z pixels')
figure; imagesc(squeeze(CBV(:, 40, :))'); colormap hot
title('CBVi - yz plane')
xlabel('x pixels')
ylabel('z pixels')
%% CBV MIP over the whole dimension
gamcp = 1; % gamma compression power
figure; imagesc(squeeze(max(CBVi, [], 1))' .^ gamcp); colormap hot; colorbar
title('CBVi - xz MIP')
xlabel('y pixels')
ylabel('z pixels')
figure; imagesc(squeeze(max(CBVi, [], 2))' .^ gamcp); colormap hot; colorbar
title('CBVi - yz MIP')
xlabel('x pixels')
ylabel('z pixels')

%% 

% volumeViewer(CBV .^ gamcp)


%% Store all the CBVi across the experiment into one matrix
load([savepath, 'fUSdata-', num2str(1), '.mat'], 'CBVi')
CBViallSF = zeros([size(CBVi), endFile - startFile + 1]); % Matrix with the CBVi for every superframe
CBViallSF(:, :, :, 1) = CBVi;
for filenum = startFile + 1:endFile
    load([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'CBVi')
    CBViallSF(:, :, :, filenum) = CBVi;
end

%% Store all the PDI across the experiment into one matrix
load([savepath, 'PDI_CDI-', num2str(1), '.mat'], 'PDI', 'CDI')
% PDIallSF = cell([length(PDI), endFile - startFile + 1]); % Matrix with the CBVi for every superframe
PDIallSF = cell([size(PDI)]); % Matrix with the CBVi for every superframe
% PDIallSF(:,  1) = PDI;
CDIallSF = cell([size(CDI)]); % Matrix with the CBVi for every superframe
% CDIallSF(:,  1) = CDI;


% for filenum = startFile + 1:endFile
for filenum = startFile:endFile
    load([savepath, 'PDI_CDI-', num2str(filenum), '.mat'], 'PDI', 'CDI')
    for i = 1:3
        PDIallSF{i} = cat(4, PDIallSF{i}, PDI{i});
        CDIallSF{i} = cat(4, CDIallSF{i}, CDI{i});
    end
end

%% Visualize the PDI and CDI across the experiment
% mwr = 30:50; % MIP window range
% mdim = 1; % MIP dimension
% newsize = size(CBViallSF);
% newsize(mdim) = 1; % Set the size of the new variable to 1
% CBViMIPStack = zeros(newsize);
PDIMIPStack = squeeze(max(PDIallSF{1}(30:50, :, :, :), [], 1));
CDIMIPStack = squeeze(max(CDIallSF{1}(30:50, :, :, :), [], 1));

%%
generateTiffStack_acrossframes(PDIallSF{3} .^ 0.7, [8.8, 8.8, 8], 'hot', 1:80)
%% Calculate rCBV and rCBF
% rCBV = CBViallSF(:, :, :, 2:end) ./ CBViallSF(:, :, :, 1:end-1);
rCBV = CBViallSF ./ CBViallSF(:, :, :, 1); % Measure relative to the "baseline", which I'm choosing as superframe 1
% Need to add the timetags

%% Visualize the CBVi across the experiment
% mwr = 30:50; % MIP window range
% mdim = 1; % MIP dimension
% newsize = size(CBViallSF);
% newsize(mdim) = 1; % Set the size of the new variable to 1
% CBViMIPStack = zeros(newsize);
CBViMIPStack = squeeze(max(CBViallSF(30:50, :, :, :), [], 1));
%% Check different MIPs across superframes
yr = 30:40;
generateTiffStack_acrossframes(CBViallSF .^ 0.7, [8.8, 8.8, 8], 'hot', yr)
% generateTiffStack_acrossframes(CBViallSF .^ 1, [8.8, 8.8, 8], 'hot', yr)
%% Plot the rCBV at some point
% Increasing y is going towards the back of the brain
% Increasing x is going from the right to the left of the brain if we align
% with -y (look towards the front)

% pt = [40, 45, 61];
pt = [40, 56, 26];
high_values_risingedges = squeeze(rCBV(pt(1), pt(2), pt(3), :));
test_ma = movmean(high_values_risingedges, 1);
% figure; plot(test, '-o')
figure; plot(test_ma, '-o')
title("rCBV at " + num2str(pt(1)) + ", " +  num2str(pt(2)) + ", " +num2str(pt(3)))
xlabel('')
ylabel('rCBV')

%% Separate each trial
ah = 3; % Approximate a cutoff value for analog high
figure; plot(TD.airPuffOutput)
ind_above_ah = find(TD.airPuffOutput > ah); % Get indices of the air puff output above analog high
ind_shift_below_ah = find(TD.airPuffOutput(ind_above_ah - 1) < ah); % See which indices above analog high have an analog low when shifted by -1 (rising edge)
ind_rising_edge = ind_above_ah(ind_shift_below_ah); % Store the original indices for the rising edges
% hold on
% plot(ind_rising_edge, ones(size(ind_rising_edge)) .* 5, 'o')
% hold off

stim_starts_gap = (P.Mcr_fcp.apis.seq_length_s - P.Mcr_fcp.apis.stim_length_s) * P.daqrate; % How long we expect the stim gap to be between the end of one stim to the start of the next
stim_prestart_baseline = (P.Mcr_fcp.apis.delay_time_ms / 1e3) * P.daqrate; % The duration between the baseline period and the corresponding stim start
stim_starts = ind_rising_edge([true; diff(ind_rising_edge) > stim_starts_gap]); % Add a 1/true at the beginning index for the first stim
hold on
plot(stim_starts, ones(size(stim_starts)) .* 5, 'o') % Plot the calculated start points of each stim period
hold off

clearvars ind_above_ah ind_shift_below_ah ind_rising_edge
% figure; plot(TD.sfTimeTagsDAQStart_adj) % plot the time tags for each superframe, adjusted to match the DAQ sampling rate

trial_windows = cell(size(stim_starts)); % Cell array of size (# trials, 1). Each cell contains the time points (according to the DAQ rate) that correspond to that trial.
trial_sf = cell(size(trial_windows));    % Cell array of size (# trials, 1). Each cell contains the superframe indices that started within that trial.

sfStarts = (TD.sfTimeTagsDAQStart_adj - TD.sfWidth_adj); % Adjust the superframe time tags so each index is at the start of the superframe acquisition

% Go through each trial within the run and assign the trial timepoints and the corresponding superframe indices
for trial = 1:length(trial_windows)
    trial_windows{trial} = stim_starts(trial) - stim_prestart_baseline : stim_starts(trial) + stim_starts_gap;

    trial_sf{trial} = find(sfStarts >= trial_windows{trial}(1) & sfStarts <= trial_windows{trial}(end));
end
clearvars trial

%% Assign the superframe trial binning to CBVi and PDI
trial_CBVi = cell(size(trial_sf));
trial_PDI = cell(size(trial_sf)); % use the all frequency PDI
minNumPts = Inf;
for trial = 1:length(trial_sf)
    trial_CBVi{trial} = CBViallSF(:, :, :, trial_sf{trial});
    trial_PDI{trial} = PDIallSF{3}(:, :, :, trial_sf{trial});
    minNumPts = min(minNumPts, length(trial_sf{trial})); % Get the minimum number of measurement points across all trials
end


%% Trial averaging
temp_size = size(trial_CBVi{1}); temp_size(4) = minNumPts; % NEED TO THINK ABOUT THE ALIGNMENT
trialAvg_CBVi = zeros(temp_size);
clearvars temp_size

for trial = 1:length(trial_sf)
    trialAvg_CBVi = trialAvg_CBVi + trial_CBVi{trial}(:, :, :, 1:minNumPts);
end
trialAvg_CBVi = trialAvg_CBVi ./ length(trial_sf);

% HRF_analytical = ;
% generateTiffStack_multi([{trialAvg_CBVi(:, :, :, 10) .^ 0.7}], [8.8, 8.8, 8], 'hot', 5)
% yr = 30:40;
yr = 1:80;
generateTiffStack_acrossframes(trialAvg_CBVi .^ 0.7, [8.8, 8.8, 8], 'hot', yr)
figure; imagesc(squeeze(max(trialAvg_CBVi(1:80, :, :, 1), [], 1))'); colormap hot

%% Plot the CBVi trial average at some point
% Increasing y is going towards the back of the brain
% Increasing x is going from the right to the left of the brain if we align
% with -y (look towards the front)

pt = [35, 26, 33];
high_values_risingedges = squeeze(trialAvg_CBVi(pt(1), pt(2), pt(3), :));
test_ma = movmean(high_values_risingedges, 1);
% figure; plot(test, '-o')
figure; plot(test_ma, '-o')
title("CBVi at " + num2str(pt(1)) + ", " +  num2str(pt(2)) + ", " +num2str(pt(3)))
xlabel('')
ylabel('CBVi')
%% Helper functions
