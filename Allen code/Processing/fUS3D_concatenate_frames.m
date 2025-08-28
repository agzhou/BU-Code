
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

%% Choose which superframes to concatenate
sframes = [40:42]; % Superframes to use
ncf = length(sframes); % # of concatenated superframes

%%
IQcat = [];
for sfi = sframes
    load([IQpath, IQfilenameStructure, num2str(sfi)])
    IQ = single(squeeze(IData + 1i .* QData));
    clearvars IData QData

    IQcat = cat(4, IQcat, IQ);
end
%%
%     IQm = IQ(:, :, 40:end, :);
%     IQm = IQ(:, :, 50:end, :);
    % IQm = IQ(:, :, 1:100, :);
%     figure; imagesc(squeeze(max(abs(IQm(:, :, :, 2)), [], 1))')

    %%%%%%%%%%%%%% IF USING THE MASK %%%%%%%%%%%%
    IQmcat = IQcat;
    coronal_mask_rep_cat = repmat(coronal_mask_rep, 1, 1, 1, ncf);
    IQmcat(coronal_mask_rep_cat) = 0; % Apply the brain mask to the IQ: set the non-brain voxels equal to 0
    IQmcat = IQmcat(:, :, 40:end, :); % Crop to save memory

    % Apply the HPF
    dimcat = length(size(IQmcat)); % Operate on the time dimension
    IQmcat_HPF = filter(HPF_b, HPF_a, IQmcat, [], dimcat);

    % Determine the optimal SV thresholds with the spatial similarity matrix
    [xpcat, ypcat, zpcat, nfcat] = size(IQmcat);
    PPcat = reshape(IQmcat, [xpcat*ypcat*zpcat, nfcat]);
    tic
%     [U, S, V] = svd(PP); % Already sorted in decreasing order
    [Ucat, Scat, Vcat] = svd(PPcat, 'econ'); % Already sorted in decreasing order
    SVscat = diag(Scat);
%     disp('Full SVD done')
%     toc
% 
%     % Plot one SVD subspace as an image
%     subspace = 180;
%     subspace_img = reshape(Ucat(:, subspace) * SVscat(subspace) * Vcat(:, subspace)', [xpcat, ypcat, zpcat, nfcat]);
%     figure; imagesc(squeeze(max(abs(subspace_img(:, :, :, 2)), [], 1))')
%     volumeViewer(abs(subspace_img(:, :, :, 2)))
% 
%     SSMcat = plotSSM(Ucat, false);
    SSMcat = plotSSM(Ucat, true);
%     [~, a_opt, b_opt] = fitSSM(SSM, false); % Get the optimal singular value thresholds
% %     [~, a_opt, b_opt] = fitSSM(SSM, true); % Get the optimal singular value thresholds
    

    [PPcat, EVscat, V_sortcat] = getSVs2D(IQmcat);
%     disp('SVs decomposed')
%     [IQf_HPF, noise] = applySVs2D(IQmcat_HPF, PP, SVscat, Vcat, a_opt, b_opt);
    [IQfcat, noisecat] = applySVs2D(IQmcat, PPcat, SVscat, Vcat, sv_threshold_lower, sv_threshold_upper);
%     [IQfcat, noisecat] = applySVs2D(IQmcat, PPcat, EVscat, V_sortcat, sv_threshold_lower_cat, sv_threshold_upper_cat);
%     disp('SVD filtered images put together')

%     volumeViewer(abs(IQfcat(:, :, :, 1)))
%     figure; imagesc(squeeze(abs(max(IQfcat(:, :, :, 1), [], 1)))')
    % clearvars IQ

    % Use the IQf with separated negative and positive frequency components
%     [IQf_separated, IQf_FT_separated] = separatePosNegFreqs(IQf);
    
    numg1pts = 20; % Only calculate the first N points
%     g1_n = g1T(IQf_separated{1}, numg1pts);
% %     [CBFsi_n, CBVi_n] = g1_to_CBi(g1_n, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
%     g1_p = g1T(IQf_separated{2}, numg1pts);
%     [CBFsi_p, CBVi_p] = g1_to_CBi(g1_p, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
% 
    g1 = g1T(IQfcat, numg1pts);

%     [PDI] = calcPowerDoppler(IQf_separated);
    PDIcat = sum(abs(IQfcat) .^ 2, 4) ./ size(IQfcat, 4);
%     [CDI] = calcColorDoppler(IQf_FT_separated, P);

%     figure; imagesc(squeeze(max(PDIcat, [], 1))' .^ 0.5); colormap hot
    figure; imagesc(squeeze(max(PDIcat ./ noisecat, [], 1))' .^ 0.5); colormap hot
%     volumeViewer(PDIcat)

%     save([savepath, 'PDI_CDI-', num2str(filenum), '.mat'], 'PDI', 'CDI', '-v7.3', '-nocompression');
%     disp("PDI and CDI for file " + num2str(filenum) + " saved" )
%     save([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'g1', 'CBFsi', 'CBVi', 'PDI', 'CDI', '-v7.3', '-nocompression');
%     save([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'g1', 'CBFsi', 'CBVi', 'PDI', 'CDI', 'g1_n', 'g1_p', 'CBFsi_n', 'CBVi_n', 'CBFsi_p', 'CBVi_p',  '-v7.3', '-nocompression');
%     save([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'g1', 'g1_n', 'g1_p', 'PDI', 'CDI', '-v7.3', '-nocompression');
    save([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'g1', 'PDIcat', 'noisecat', 'SVscat', 'SSM', 'a_opt', 'b_opt', '-v7.3', '-nocompression');
%     save([savepath, 'g1-', num2str(filenum), '.mat'], 'g1', 'g1_n', 'g1_p', '-v7.3', '-nocompression');
%     save([savepath, 'g1-', num2str(filenum), '.mat'], 'g1', '-v7.3', '-nocompression');

    disp("fUS result for file " + num2str(filenum) + " saved" )
%     disp("g1 result for file " + num2str(filenum) + " saved" )

%%
g1cat = {};
g1cat{1} = g1T(IQfcat(:, :, :, 1:496), numg1pts);
g1cat{2} = g1T(IQfcat(:, :, :, 497:992), numg1pts);
g1cat{3} = g1T(IQfcat(:, :, :, 993:end), numg1pts);

CBVicat = {};
CBFsicat = {};
[CBFsicat{1}, CBVicat{1}] = g1_to_CBi(g1cat{1}, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
[CBFsicat{2}, CBVicat{2}] = g1_to_CBi(g1cat{2}, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
[CBFsicat{3}, CBVicat{3}] = g1_to_CBi(g1cat{3}, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)

figure; imagesc(squeeze(max(CBVicat{1} .^ 0.5, [], 1))'); colormap hot
figure; imagesc(squeeze(max(CBVicat{2} .^ 0.5, [], 1))'); colormap hot
figure; imagesc(squeeze(max(CBVicat{3} .^ 0.5, [], 1))'); colormap hot