%% Description:
%       3D (tl-)fUS processing
%       Timing data should be processed with plotfUStiming.m first

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
[timingFilePathFN, timingFilePath] = uigetfile([IQpath, '..\Timing data\TD.mat'], 'Select the timing data');
timingFilePath = [timingFilePath, timingFilePathFN];
load(timingFilePath)
% load(timingFilePath, 'acqStart', 'airPuffOutput', 'daqStartTimetag', 'sfTimeTags', 'sfTimeTagsDAQStart', 'sfTimeTagsDAQStart_adj', 'sfWidth', 'sfWidth_adj', 'timeStamp')
%% Define some parameters

parameterPrompt = {'Start file number', 'End file number', 'SVD lower bound', 'SVD upper bound', 'Tau 1 index for CBFspeed', 'Tau 2 index for CBFspeed', 'Tau 1 index for CBV'};
parameterDefaults = {'1', '', '20', '', '2', '5', '2'};
parameterUserInput = inputdlg(parameterPrompt, 'Input Parameters', 1, parameterDefaults);

% define # of files manually for now
% str2double(parameterUserInput{});
startFile = str2double(parameterUserInput{1});
endFile = str2double(parameterUserInput{2});
numFiles = endFile - startFile + 1;
sv_threshold_lower = str2double(parameterUserInput{3});
sv_threshold_upper = str2double(parameterUserInput{4});
tau1_index_CBF = str2double(parameterUserInput{5});
tau2_index_CBF = str2double(parameterUserInput{6});
tau1_index_CBV = str2double(parameterUserInput{7});

clearvars parameterPrompt parameterDefaults parameterUserInput

taustep = 1/P.frameRate;
% tau = taustep:taustep:(P.numFramesPerBuffer * taustep);
tau = 0:taustep:((P.numFramesPerBuffer - 1) * taustep);
tau_ms = tau .* 1000; % Assuming even time spacing between frames

% tau1_index_CBF = 2;
% tau2_index_CBF = 6;
% tau1_index_CBV = 2;

%% Main loop
% for filenum = startFile:endFile
% for filenum = [2:endFile]
% for filenum = [285:-1:189]
for filenum = 100:223
% for filenum = 1

    % Load the IQ data
    tic
    load([IQpath, IQfilenameStructure, num2str(filenum)])
    
    IQ = single(squeeze(IData + 1i .* QData));
    clearvars IData QData

    % figure; imagesc(squeeze(max(abs(IQ(:, :, :, 2)), [], 1))')
    
    % Crop the IQ first 
    zstart = 40;
    zend = size(IQ, 3);
    IQm = IQ(:, :, zstart:zend, :);
%     figure; imagesc(squeeze(max(abs(IQm(:, :, :, 2)), [], 1))')

    %%%%%%%%%%%%%% IF USING THE PREDEFINED MASK %%%%%%%%%%%%
%     IQm(coronal_mask_rep) = 0; % Apply the brain mask to the IQ: set the non-brain voxels equal to 0
    %     [xp, yp, zp, nf] = size(IQm);
    
    % SVD decluttering
    [PP, EVs, V_sort] = getSVs2D(IQm);
    disp('SVs decomposed')
    [IQf, noise] = applySVs2D(IQm, PP, EVs, V_sort, sv_threshold_lower, sv_threshold_upper);
    disp('SVD filtered images put together')

%     volumeViewer(abs(IQf(:, :, :, 1)))
%     figure; imagesc(squeeze(abs(max(IQf(:, :, :, 1), [], 1)))')
    % clearvars IQ

    % Use the IQf with separated negative and positive frequency components
%     [IQf_separated, IQf_FT_separated] = separatePosNegFreqs(IQf);
    
    numg1pts = 20; % Only calculate the first N points
%     g1_n = g1T(IQf_separated{1}, numg1pts);
% %     [CBFsi_n, CBVi_n] = g1_to_CBi(g1_n, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
%     g1_p = g1T(IQf_separated{2}, numg1pts);
%     [CBFsi_p, CBVi_p] = g1_to_CBi(g1_p, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
% 
    g1 = g1T(IQf, numg1pts);
%     g1 = g1T(IQf);
%     [CBFsi, CBVi] = g1_to_CBi(g1, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
% 
% %     savefast([savepath, 'fUSdata-', num2str(filenum), '.mat'], g1, CBFi, CBVi);

%     [PDI] = calcPowerDoppler(IQf_separated);
    PDI = sum(abs(IQf) .^ 2, 4) ./ size(IQf, 4);
%     [CDI] = calcColorDoppler(IQf_FT_separated, P);

%     figure; imagesc(squeeze(max(PDI, [], 1))' .^ 0.5); colormap hot
%     figure; imagesc(squeeze(max(PDI ./ noise, [], 1))' .^ 0.5); colormap hot
%     volumeViewer(PDI)

%     save([savepath, 'PDI_CDI-', num2str(filenum), '.mat'], 'PDI', 'CDI', '-v7.3', '-nocompression');
%     disp("PDI and CDI for file " + num2str(filenum) + " saved" )
%     save([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'g1', 'CBFsi', 'CBVi', 'PDI', 'CDI', '-v7.3', '-nocompression');
%     save([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'g1', 'CBFsi', 'CBVi', 'PDI', 'CDI', 'g1_n', 'g1_p', 'CBFsi_n', 'CBVi_n', 'CBFsi_p', 'CBVi_p',  '-v7.3', '-nocompression');
%     save([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'g1', 'g1_n', 'g1_p', 'PDI', 'CDI', '-v7.3', '-nocompression');
    save([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'g1', 'PDI', 'noise', '-v7.3', '-nocompression');
%     save([savepath, 'g1-', num2str(filenum), '.mat'], 'g1', 'g1_n', 'g1_p', '-v7.3', '-nocompression');
%     save([savepath, 'g1-', num2str(filenum), '.mat'], 'g1', '-v7.3', '-nocompression');

    disp("fUS result for file " + num2str(filenum) + " saved" )
%     disp("g1 result for file " + num2str(filenum) + " saved" )

    toc
    
end
% savefast([savepath, 'fUS_proc_params.mat'], 'sv_threshold_lower', 'sv_threshold_upper', 'tau', 'tau_ms', 'tau1_index_CBF', 'tau2_index_CBF', 'tau1_index_CBV');
savefast([savepath, 'fUS_proc_params.mat'], 'sv_threshold_lower', 'sv_threshold_upper', 'tau', 'tau_ms', 'numg1pts');
% savefast([savepath, 'PDI_CDI_proc_params.mat'], 'sv_threshold_lower', 'sv_threshold_upper');

%% Testing
testIQf = applySVs2D(IQ, PP, EVs, V_sort, 35, sv_threshold_upper);
figure; imagesc(squeeze(max(abs(testIQf(30:50, :, :, 1)), [], 1) .^ 0.5)'); colormap hot
% volumeViewer(abs(testIQf(:, :, :, 1)))
testPDI = mean(abs(IQf) .^ 2, 4);
figure; imagesc(squeeze(max(testPDI(30:50, :, :, 1), [], 1) .^ 0.5)'); colormap hot
figure; imagesc(squeeze(max(testPDI(:, :, :, 1), [], 3) .^ 0.5)'); colormap hot

%% Convert g1 into CBV, CBFspeed, etc.

g1_tau1_cutoff = 0.3;
% g1_tau1_cutoff = 0.0;
% tau_difference_cutoff = 0.2;

for filenum = startFile:endFile
% for filenum = [1]
%     load([savepath, 'g1-', num2str(filenum)], 'g1') % Load the saved g1 mat files
    load([savepath, 'fUSdata-', num2str(filenum)], 'g1') % Load the saved g1 mat files

    [g1A_mask] = createg1mask(g1, g1_tau1_cutoff, tau1_index_CBF, tau2_index_CBF);
%     [g1A_mask] = createg1mask(g1Avg, g1_tau1_cutoff, tau1_index_CBF, tau2_index_CBF);

    [CBFsi, CBVi] = g1_to_CBi(g1, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)

%     CBFsi(~g1A_mask) = -Inf; % Remove noisy points from the CBFspeed index (in theory)
    CBFsi(~g1A_mask) = 0; % Remove noisy points from the CBFspeed index (in theory)

    save([savepath, 'tlfUSdata-', num2str(filenum), '.mat'], 'CBFsi', 'CBVi', '-v7.3', '-nocompression');
%     save([savepath, 'tlfUSdatatest-', num2str(filenum), '.mat'], 'CBFsi', 'CBVi', '-v7.3', '-nocompression');
    disp("tl-fUS result for file " + num2str(filenum) + " saved" )

end
save([savepath, 'tlfUS_proc_params.mat'], 'tau1_index_CBV', 'tau1_index_CBF', 'tau2_index_CBF', 'g1_tau1_cutoff', 'g1A_mask');
% save([savepath, 'tlfUStest_proc_params.mat'], 'tau1_index_CBV', 'tau1_index_CBF', 'tau2_index_CBF', 'g1_tau1_cutoff');
figure; imagesc(squeeze(max(CBVi(30:50, :, :), [], 1) .^ 0.3)'); colormap hot
figure; imagesc(squeeze(max(CBVi(:, :, :), [], 3) .^ 0.5)'); colormap hot
vcmap = colormap_ULM;
figure; imagesc(squeeze(mean(CBFsi(30:50, :, :), 1))'); colormap(vcmap)

% generateTiffStack_multi({CBVi .^ 0.7}, [8.8, 8.8, 8], 'hot', 5)
%% Get and save PDI, CDI only
for filenum = startFile:endFile
% for filenum = 1
% for filenum = [285:-1:189]
    tic
    load([IQpath, IQfilenameStructure, num2str(filenum)])
%     load(['E:\Allen BME-BOAS-27 Data Backup\AZ03 Stroke RC15gV\fUS\05-06-2025 pre-stroke\IQ Data - Verasonics recon\', IQfilenameStructure, num2str(filenum)])
    
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
%     [IQf_separated, IQf_FT_separated] = separatePosNegFreqs(IQf);
    
    [PDI] = calcPowerDoppler(IQf);
%     [CDI] = calcColorDoppler(IQf_FT_separated, P);

%     save([savepath, 'PDI_CDI-', num2str(filenum), '.mat'], 'PDI', 'CDI', '-v7.3', '-nocompression');
    save([savepath, 'PDI-', num2str(filenum), '.mat'], 'PDI', '-v7.3', '-nocompression');
    disp("PDI and CDI for file " + num2str(filenum) + " saved" )

    toc
    
end
savefast([savepath, 'PDI_CDI_proc_params.mat'], 'sv_threshold_lower', 'sv_threshold_upper');


%% CBFspeed index with the Derivative Method
% vtest = sqrt( (abs(g1(:, :, :, tau1_index_CBF)) - abs(g1(:, :, :, tau2_index_CBF))) ./ (abs(g1(:, :, :, tau1_index_CBF)) .* ( tau(tau2_index_CBF)^2 - tau(tau1_index_CBF)^2 )) );
CBFspeed_DM = (abs(g1(:, :, :, tau1_index_CBF)) - abs(g1(:, :, :, tau2_index_CBF))) ./ (abs(g1(:, :, :, tau1_index_CBF)) .* ( tau(tau2_index_CBF)^2 - tau(tau1_index_CBF)^2 ));
CBFspeed_DM(CBFspeed_DM < 0) = 0;
CBFspeed_DM = sqrt(CBFspeed_DM);

CBFspeed_DM_masked = CBFspeed_DM;
CBFspeed_DM_masked(~g1A_mask) = 0;

% vtest_NLMF = NLMF(vtest ./ max(vtest, [], 'all')); % Try using nonlinear means filtering to get rid of those high-noise points
% vtest_gaussLPF = imgaussfilt3(vtest);

%% Random phase test
phasetest = atan2(imag(IQf), real(IQf));

figure; imagesc(squeeze(phasetest(40, :, :, 1))'); % colormap(vcmap)
figure; imagesc(squeeze(abs(IQf(40, :, :, 1)))'); % colormap(vcmap)

% figure; imagesc(squeeze(mean(phasetest(30:50, :, :), 1))'); % colormap(vcmap)

%%
figure; imagesc(squeeze(max(CBFspeed_DM_masked(30:50, :, :), [], 1))'); colormap(vcmap)
figure; imagesc(squeeze(mean(CBFspeed_DM_masked(30:50, :, :), 1))'); colormap(vcmap)
% figure; imagesc(squeeze(median(CBFspeed_DM_masked(30:50, :, :), 1))'); colormap(vcmap)
% volumeViewer(vtest_masked)

% generateTiffStack_MeanIPs_multi({CBFspeed_DM_masked .^ 1}, [8.8, 8.8, 8], vcmap, 10)
% generateTiffStack_multi({vtest_masked .^ 1}, [8.8, 8.8, 8], vcmap, 10)
% figure; imagesc(squeeze(max(vtest_NLMF(30:50, :, :), [], 1))'); colormap(vcmap)
% figure; imagesc(squeeze(mean(vtest_NLMF(30:50, :, :), 1))'); colormap(vcmap)
% volumeViewer(vtest_NLMF)

% figure; imagesc(squeeze(max(vtest_gaussLPF(30:50, :, :), [], 1))'); colormap(vcmap)
% figure; imagesc(squeeze(mean(vtest_gaussLPF(30:50, :, :), 1))'); colormap(vcmap)
% figure; imagesc(squeeze(median(vtest_gaussLPF(30:50, :, :), 1))'); colormap(vcmap)

% volumeViewer(vtest_gaussLPF)

%% Compare the derivative CBFspeed index method to the old log one
CBFsi_all_masked = CBFsi_all;
CBFsi_all_masked(~g1A_mask) = 0;

figure; imagesc(squeeze(max(CBFsi_all_masked(30:50, :, :), [], 1))'); colormap(vcmap)
figure; imagesc(squeeze(mean(CBFsi_all_masked(30:50, :, :), 1))'); colormap(vcmap)

% generateTiffStack_MeanIPs_multi({CBFsi_all_masked .^ 1}, [8.8, 8.8, 8], vcmap, 10)
% generateTiffStack_multi({CBFsi_all_masked .^ 1}, [8.8, 8.8, 8], vcmap, 10)

generateTiffStack_MeanIPs_multi({CBFsi_all_masked .^ 1, CBFspeed_DM_masked}, [8.8, 8.8, 8], vcmap, 20)

%% plot the base g1 at some point
% pt = [40, 58, 14];
% pt = [32, 34, 17];
pt = [32, 31, 35];
% pt = [40, 36, 128];
% figure; plot(tau_ms(1:size(g1, 4)), squeeze(abs(g1(pt(1), pt(2), pt(3), :))), '-o');
figure; plot(tau_ms(1:size(g1s{1}, 4)), squeeze(abs(g1s{1}(pt(1), pt(2), pt(3), :))), '-o');
hold on
plot(tau_ms(1:size(g1Avg, 4)), squeeze(abs(g1Avg(pt(1), pt(2), pt(3), :))), '-o');
% plot(tau_ms(1:size(g1s{2}, 4)), squeeze(abs(g1s{2}(pt(1), pt(2), pt(3), :))), '-o');
hold off
xlabel('tau [ms]')
ylabel('|g1|')

%% Use smoothdata
g1_smoothdata = smoothdata(g1, 4, 'sgolay'); % Smooth along the time dimension
%
figure; plot(tau_ms(1:size(g1, 4)), squeeze(abs(g1_smoothdata(pt(1), pt(2), pt(3), :))), '-o');
xlabel('tau [ms]')
ylabel('|g1| with smoothdata')

%% Trying some alternate CBFspeed index calculations
% CBFsi_all = squeeze(abs(g1(:, :, :, tau1_index_CBF)) - abs(g1(:, :, :, tau2_index_CBF)));
% CBFsi_p = squeeze(abs(g1_p(:, :, :, tau1_index_CBF)) - abs(g1_p(:, :, :, tau2_index_CBF)));
% CBFsi_n = squeeze(abs(g1_n(:, :, :, tau1_index_CBF)) - abs(g1_n(:, :, :, tau2_index_CBF)));
% volumeViewer(CBFsi_all)
% figure; imagesc(squeeze(max(CBFsi_all(30:50, :, :), [], 1))')

vcmap = colormap_ULM; % velocity colormap
% generateTiffStack_multi({CBFsi_all}, [8.8, 8.8, 8], vcmap, 5)
% generateTiffStack_multi({CBFsi_p}, [8.8, 8.8, 8], vcmap, 5)
% generateTiffStack_multi({CBFsi_n}, [8.8, 8.8, 8], vcmap, 5)

% CBFsi_all_test = squeeze(abs(g1(:, :, :, tau1_index_CBF)) ./ abs(g1(:, :, :, tau2_index_CBF))); %./ (tau2_index_CBF - tau1_index_CBF);
% CBFsi_all_test = squeeze(log(abs(g1(:, :, :, tau1_index_CBF)) ./ abs(g1(:, :, :, tau2_index_CBF)))); %./ (tau2_index_CBF - tau1_index_CBF);
% CBFsi_all_test = squeeze(abs(g1(:, :, :, tau1_index_CBF)) ./ abs(g1(:, :, :, tau2_index_CBF))); %./ (tau2_index_CBF - tau1_index_CBF);
CBFsi_all_test = squeeze(log(abs(g1(:, :, :, tau1_index_CBF))) - log(abs(g1(:, :, :, tau2_index_CBF)))); %./ (tau2_index_CBF - tau1_index_CBF);
figure; imagesc(squeeze(max(CBFsi_all_test(30:50, :, :), [], 1))'); colormap(vcmap)

figure; imagesc(squeeze(mean(CBFsi_all_test(30:50, :, :), 1))'); colormap(vcmap)

volumeViewer(CBFsi_all_test)
%% g1 adjustment test (see MAIN_g1fUS_invivo_annotated.m)
g1_shift = g1(:, :, :, 2:end);
g1Temp = reshape(g1_shift, [size(g1_shift, 1) * size(g1_shift, 2) * size(g1_shift, 3), size(g1_shift, 4)]); % ** stack the spatial dimensions of the g1 **
                
% ** This seems like noise removal (CR = Clutter Removal? CRiteria?) **

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
      
%%
figure; plot(tau_ms(1:size(g1Adj, 4)), squeeze(abs(g1Adj(pt(1), pt(2), pt(3), :))), '-o');
xlabel('tau [ms]')
ylabel('|g1|')
%% Plot CBVi and CBFspeedi from the smoothed g1
tau1_index_CBF = 2;
tau2_index_CBF = 10;
tau1_index_CBV = 2;

[CBFi_adj, CBVi_adj] = g1_to_CBi(g1Adj, tau_ms(2:end) ./ 1000, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
%
figure; imagesc(squeeze(mean(CBFi_adj(30:50, :, :), 1))')
% plotMIPs(CBFi_adj, 1)
% plotMIPs(CBVi_adj, 1)
% plotMIPs(CBVi, 1)

%% Look at the g1 adjustment criteria
ggCR1_rs = reshape(ggCR1, [size(g1_shift, 1), size(g1_shift, 2), size(g1_shift, 3)]); % ** unstack the spatial dimensions **
figure; imagesc(squeeze(max(ggCR1_rs, [], 1))')

ggCR2_rs = reshape(ggCR2, [size(g1_shift, 1), size(g1_shift, 2), size(g1_shift, 3)]); % ** unstack the spatial dimensions **
figure; imagesc(squeeze(max(ggCR2_rs, [], 1))') %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% More smoothing from the previous code
covB=ones(3,3); covB(3,3)=9; covB=covB/sum(covB(:));

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
figure; imagesc(squeeze(max(PDI{3}(30:50, :, :), [], 1) .^ 0.7)'); colormap hot
% generateTiffStack_multi({PDI{3} .^ 0.9}, [8.8, 8.8, 8], 'hot', 5)
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
plotMIPs(CBFsi, 1)

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
figure; imagesc(squeeze(max(CBFsi, [], 2))' .^ 1); colormap hot
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

%% Store all the updated CBVi and CBFsi across the experiment into one matrix
load([savepath, 'tlfUSdata-', num2str(1), '.mat'], 'CBFsi', 'CBVi')
CBViallSF = zeros([size(CBVi), endFile - startFile + 1]); % Matrix with the CBVi for every superframe
CBViallSF(:, :, :, 1) = CBVi;

CBFsiallSF = zeros([size(CBFsi), endFile - startFile + 1]); % Matrix with the CBFsi for every superframe
CBFsiallSF(:, :, :, 1) = CBFsi;

for filenum = startFile + 1:endFile
    load([savepath, 'tlfUSdata-', num2str(filenum), '.mat'], 'CBFsi', 'CBVi')
    CBViallSF(:, :, :, filenum) = CBVi;
    CBFsiallSF(:, :, :, filenum) = CBFsi;
end

%% Store all the PDI across the experiment into one matrix (with separated frequencies)
% % load([savepath, 'PDI_CDI-', num2str(1), '.mat'], 'PDI', 'CDI')
% load([savepath, 'fUSdata-', num2str(1), '.mat'], 'PDI', 'CDI')
% % PDIallSF = cell([length(PDI), endFile - startFile + 1]); % Matrix with the CBVi for every superframe
% PDIallSF = cell([size(PDI)]); % Matrix with the CBVi for every superframe
% % PDIallSF(:,  1) = PDI;
% CDIallSF = cell([size(CDI)]); % Matrix with the CBVi for every superframe
% % CDIallSF(:,  1) = CDI;
% 
% % for filenum = startFile + 1:endFile
% for filenum = startFile:endFile
% %     load([savepath, 'PDI_CDI-', num2str(filenum), '.mat'], 'PDI', 'CDI')
%     load([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'PDI', 'CDI')
% %     PDI = load([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'PDI')
% %     CDI = load([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'CDI')
% 
%     for i = 1:3
% %         PDIallSF1 = cat(4, PDIallSF1, PDI{1});
% %         CDIallSF1 = cat(4, CDIallSF1, CDI{1});
% %         PDIallSF2 = cat(4, PDIallSF2, PDI{2});
% %         CDIallSF2 = cat(4, CDIallSF2, CDI{2});
% %         PDIallSF3 = cat(4, PDIallSF3, PDI{3});
% %         CDIallSF3 = cat(4, CDIallSF3, CDI{3});
%         PDIallSF{i} = cat(4, PDIallSF{i}, PDI{i});
%         CDIallSF{i} = cat(4, CDIallSF{i}, CDI{i});
%     end
% end

%% Store all the PDI across the experiment into one matrix
% load([savepath, 'PDI_CDI-', num2str(1), '.mat'], 'PDI', 'CDI')
% load([savepath, 'fUSdata-', num2str(1), '.mat'], 'PDI', 'CDI')
load([savepath, 'fUSdata-', num2str(1), '.mat'], 'PDI')
PDIallSF = zeros([size(PDI), endFile - startFile + 1]); % Matrix with the CBVi for every superframe
PDIallSF(:, :, :, 1) = PDI;

% for filenum = startFile + 1:endFile
for filenum = startFile:endFile
%     load([savepath, 'PDI_CDI-', num2str(filenum), '.mat'], 'PDI', 'CDI')
    load([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'PDI')

    PDIallSF(:, :, :, filenum) = PDI;
end

%% Visualize the PDI and CDI across the experiment
% mwr = 30:50; % MIP window range
% mdim = 1; % MIP dimension
% newsize = size(CBViallSF);
% newsize(mdim) = 1; % Set the size of the new variable to 1
% CBViMIPStack = zeros(newsize);
PDIMIPStack = squeeze(max(PDIallSF{1}(30:50, :, :, :), [], 1));
CDIMIPStack = squeeze(max(CDIallSF{1}(30:50, :, :, :), [], 1));

%
generateTiffStack_acrossframes(PDIallSF{3} .^ 0.7, [8.8, 8.8, 8], 'hot', 1:80)

%% Visualize the CBVi across the experiment
% mwr = 30:50; % MIP window range
% mdim = 1; % MIP dimension
% newsize = size(CBViallSF);
% newsize(mdim) = 1; % Set the size of the new variable to 1
% CBViMIPStack = zeros(newsize);
CBViMIPStack = squeeze(max(CBViallSF(30:50, :, :, :), [], 1));
%% Check different MIPs across superframes
yr = 20:40;
generateTiffStack_acrossframes(CBViallSF .^ 0.7, [8.8, 8.8, 8], 'hot', yr)
% generateTiffStack_acrossframes(CBViallSF .^ 1, [8.8, 8.8, 8], 'hot', yr)


%% Separate each trial
ah = 3; % Approximate a cutoff value for analog high

ind_above_ah = find(TD.airPuffOutput > ah); % Get indices of the air puff output above analog high
ind_shift_below_ah = find(TD.airPuffOutput(ind_above_ah - 1) < ah); % See which indices above analog high have an analog low when shifted by -1 (rising edge)
ind_rising_edge = ind_above_ah(ind_shift_below_ah); % Store the original indices for the rising edges
% hold on
% plot(ind_rising_edge, ones(size(ind_rising_edge)) .* 5, 'o')
% hold off

stim_starts_gap = (P.Mcr_fcp.apis.seq_length_s - P.Mcr_fcp.apis.stim_length_s) * P.daqrate; % How long we expect the stim gap to be between the end of one stim to the start of the next
stim_prestart_baseline = (P.Mcr_fcp.apis.delay_time_ms / 1e3) * P.daqrate; % The duration between the baseline period and the corresponding stim start
stim_starts = ind_rising_edge([true; diff(ind_rising_edge) > stim_starts_gap]); % Add a 1/true at the beginning index for the first stim

% Plot the air puff signal and the calculated start points of each stim period
figure; plot(TD.airPuffOutput)
hold on
plot(stim_starts, ones(size(stim_starts)) .* 5, 'o')
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

%% Resample the trials for the hemodynamic parameters
interp_factor = 100;
[trial_CBVi_usi] = resampleTrials(CBViallSF, trial_sf, trial_windows, sfStarts, P, interp_factor);
[trial_PDI_usi] = resampleTrials(PDIallSF, trial_sf, trial_windows, sfStarts, P, interp_factor);

%% Calculate the relative hemodynamic changes for each trial

[trial_CBVi_usi_baseline, trial_rCBV_usi] = fUS_calc_rHP(trial_CBVi_usi, P, interp_factor);
[trial_PDI_usi_baseline, trial_rPDI_usi] = fUS_calc_rHP(trial_PDI_usi, P, interp_factor);

%% Trial average the relative hemodynamic changes

rCBV_TA = fUS_trialAverage(trial_rCBV_usi);
rPDI_TA = fUS_trialAverage(trial_rPDI_usi);

%% Correlation on the trial average

% Resample the stim pattern/predicted HRF
trial_stim_pattern = zeros(P.Mcr_fcp.apis.seq_length_s * P.daqrate / interp_factor, 1);
trial_stim_pattern(P.Mcr_fcp.apis.delay_time_ms/1000 * P.daqrate / interp_factor : ...
    P.Mcr_fcp.apis.delay_time_ms/1000 * P.daqrate / interp_factor + ...
    P.Mcr_fcp.apis.stim_length_s * P.daqrate / interp_factor) = 1;
figure; plot(trial_stim_pattern); title('Trial stim pattern')

zt = 1;
[r_rCBV, z_rCBV, am_rCBV] = activationMap3D(rCBV_TA, trial_stim_pattern, zt);

volumeViewer(r_rCBV)
volumeViewer(z_rCBV)
volumeViewer(am_rCBV)
figure; imagesc(squeeze(max(r_rCBV(:, :, :), [], 1))'); colormap jet; clim([-1, 1])
% figure; imagesc(z_rCBV)
% figure; imagesc(am_rCBV); colormap jet; title("Activation Map (rCBV) with z threshold = " + num2str(zt))
figure; imagesc(squeeze(max(am_rCBV(:, :, :), [], 1))'); colormap jet; title("Activation Map (rCBV) with z threshold = " + num2str(zt))

%% Plot activation at each slice
for slice = 1:10
    my_inds = (slice-1)*5:slice*5;
    my_inds = my_inds+1;
    figure; imagesc(squeeze(mean(test(my_inds, :, :), 1))')
    title(num2str(mean(my_inds)))
end

kernel = ones(3, 3, 3);
kernel(2, 2, 2) = 3;%sum(kernel, 'all');

test_r_CBVi_relative_change_conv = convn(test_r_CBVi_relative_change, kernel, 'same');

%% Helper functions

function [g1A_mask] = createg1mask(g1, g1_tau1_cutoff, tau1_index_CBF, tau2_index_CBF)

    g1A_T = {};
    
    g1A_T{1} = abs(g1(:, :, :, 2)) > g1_tau1_cutoff; % First treatment: tau1 is above some cutoff (make sure there is some actual blood signal there)
    g1A_T{2} = abs(g1(:, :, :, tau1_index_CBF)) > abs(g1(:, :, :, tau2_index_CBF)); % Keep the voxels where |g1(tau1)| > |g1(tau2)| (noise might have the g1 randomly increase with tau, but it should not happen with a voxel where there is a real blood signal)
    g1A_T{3} = abs(g1(:, :, :, tau1_index_CBF)) > 2 .* abs(g1(:, :, :, tau2_index_CBF)); % Keep the voxels where |g1(tau1)| > 2 * |g1(tau2)| (same as #2, but more severe)
    % g1A_T{4} = abs(g1(:, :, :, tau1_index_CBF)) - 1 .* abs(g1(:, :, :, tau2_index_CBF)) > tau_difference_cutoff; % Keep the voxels where |g1(tau1)| > 2 * |g1(tau2)| (same as #2, but more severe)
    
    g1A_mask = true(size(g1A_T{1})); % Mask of voxels to keep for the g1 treatments
    for i = 1:length(g1A_T)
        g1A_mask = and(g1A_mask, g1A_T{i});
    end

end

