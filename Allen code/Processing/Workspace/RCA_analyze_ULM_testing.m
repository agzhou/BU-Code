%% Description
% ULM data analysis:
% Take reconstructed data from individual buffers/batches/superframes, each with
% some number of frames/subframes, each containing N acquisitions/angles.

% Re-sample with a rolling method to get more effective frames.

% SVD to separate the bubble signals from the tissue signal and other
% clutter

% ...
%% Use parallel processing for speed
% https://www.mathworks.com/matlabcentral/answers/91744-how-can-i-check-if-matlabpool-is-running-when-using-parallel-computing-toolbox

pp = gcp('nocreate');
if isempty(pp)
    % There is no parallel pool
    parpool LocalProfile1

end

%% Load parameters and make folder for saving the processed data
datapath = 'G:\Allen\Data\01-29-2025 AZ001 ULM\RC15gV\run 1 left eye\';
if ~exist('P', 'var')
    load([datapath, 'params.mat'])
end

numFiles = 96; % define # of files manually for now

IQfolderName = 'IQ Data - Verasonics Recon\'; % 'IQ data\'
saveFolderName = 'Processed Data 03-14-2025 2 refinement inside maskhole\';
% savepath = [datapath, saveFolderName];
% mkdir([datapath, saveFolderName])
savepath = ['F:\Allen\Data\01-29-2025 AZ001 ULM\RC15gV\run 1 left eye\', saveFolderName];
mkdir(savepath)
% savepath = 'F:\Allen\Data\01-29-2025 AZ001 ULM\RC15gV\run 1 left eye\FMAS Processed Data\';

filename_structure = ['IQ-', num2str(P.maxAngle), '-', num2str(P.na), '-', num2str(P.frameRate), '-', num2str(P.numFramesPerBuffer), '-1-'];

addpath('C:\Users\BOAS-US\Documents\Allen\GitHub\BU-Code\Allen code\Processing\normxcorr3.m')
addpath('C:\Users\BOAS-US\Documents\Allen\GitHub\BU-Code\Allen code\Processing\toolbox_nlmeans_version2')
%% Parameters for processing the data
% Define various processing parameters
% Singular value thresholds
sv_threshold_lower = 10;
sv_threshold_upper = 150;

% % Region of interest
xrange = int16(1:80);
yrange = int16(1:80);
zrange = int16(36:142);

% framerange = 1:200;
% framerange = 1:size(IQf, 3);
% range = {xrange, yrange, zrange, framerange};
range = {xrange, yrange, zrange};

% % Image refinement and localization parameters
irfc = 2;
% imgRefinementFactor = [2, 2, 2]; % z, x pixel refinement factor
imgRefinementFactor = ones(1, 3) .* irfc;

xpix_spacing = P.Trans.spacingMm / 1e3;
ypix_spacing = P.Trans.spacingMm / 1e3;
zpix_spacing = P.wl / 2;
% imgRefinementFactor = [irfc * xpix_spacing/zpix_spacing, irfc * ypix_spacing/zpix_spacing, irfc];

XCThresholdFactor = 0.25;
 
% Load and refine simulated PSF
if ~exist('PSF', 'var')
    load('G:\Allen\Data\01-29-2025 AZ001 ULM\RC15gV\PSF sim\PSF.mat', 'PSF')
    % figure; imagesc(squeeze(abs(PSF(40, :, :)))')
end

% PSFs = PSF(190:210, 118:138, :); % PSF section
PSFs = PSF(30:50, 30:50, 92:110); % PSF section
refPSF = imresize3(PSFs, [size(PSFs, 1) * imgRefinementFactor(1), size(PSFs, 2) * imgRefinementFactor(2), size(PSFs, 3) * imgRefinementFactor(3)]);
% volumeViewer(abs(refPSF))

%% Video writer

vo_ysum = VideoWriter([savepath, 'centers_ysum']);

vo_ysum.Quality = 100;
vo_ysum.FrameRate = 60;
open(vo_ysum);

vo_xsum = VideoWriter([savepath, 'centers_xsum']);

vo_xsum.Quality = 100;
vo_xsum.FrameRate = 60;
open(vo_xsum);

vo_zsum = VideoWriter([savepath, 'centers_xsum']);
vo_zsum.Quality = 100;
vo_zsum.FrameRate = 60;
open(vo_zsum);

vf_ysum = figure; colormap gray
vf_xsum = figure; colormap gray
vf_zsum = figure; colormap gray

findfigs


%% Process the data
IQf_zum_th = 0.35; % normalized threshold on the IQf z sum

for filenum = 47:numFiles
% for filenum = 30
    tic
    load([datapath, IQfolderName, filename_structure, num2str(filenum), '.mat'])  % load each reconstructed buffer/batch/superframe
%     IQr = LA_rollingFrames(IQ);                                                 % rolling method to get more effective frames
    
    IQ = squeeze(IData + 1i .* QData);   % Combine I and Q, which are saved separately. It's easier to save the big reconstructed data with savefast, which doesn't support complex values. The data is already a coherent sum.
    clear IData QData
    
    % SVD proc part 1
%     tic
    [PP, EVs, V_sort] = getSVs2D(IQ);
    disp('SVs decomposed')
%     toc
    % SVD proc part 2
%     tic
    [IQf] = applySVs2D(IQ, PP, EVs, V_sort, sv_threshold_lower, sv_threshold_upper);
    disp('SVD filtered images put together')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find a mask with the filtered IQ on the whole volume
    IQf_zsum = squeeze(sum(abs(IQf), 3)); % sum across z
    max_IQf_zsum = NaN(size(IQf_zsum, 3), 1);
    parfor f = 1:size(IQf_zsum, 3)
        max_IQf_zsum(f) = max(IQf_zsum(:, :, f), [], 'all');
    end
    max_IQf_zsum_M = ones(size(IQf, 1), size(IQf, 2), length(max_IQf_zsum));
    parfor f = 1:length(max_IQf_zsum)
%         max_IQf_zsum_M(:, :, f) = max_IQf_zsum_M(:, :, f) .* max_IQf_zsum;
        max_IQf_zsum_M(:, :, f) = max_IQf_zsum(f);
    end
    IQf_zsum_n = IQf_zsum ./ max_IQf_zsum_M; % normalize
%     figure; imagesc(IQf_zsum_n(:, :, 1))

    IQf_zsum_mask = IQf_zsum_n < IQf_zum_th;
%     figure; imagesc(IQf_zsum_mask(:, :, 1))
    % Could do the bwconncomp to remove little remnants inside the main
    % regions for this non robust thresholding
    
%     [Io, Jo] = find(~IQf_zsum_mask);
    
    IQf_zsum_mask_rm = repmat(IQf_zsum_mask, 1, 1, 1, size(IQ, 3));% repmat version so we have volume masking
    IQf_zsum_mask_rm = permute(IQf_zsum_mask_rm, [1, 2, 4, 3]); % shift the z dimension to the correct dimension
    % Get the IQ only in the hole
    IQh = IQ; IQh(~IQf_zsum_mask_rm) = 0;
%     figure; imagesc(squeeze(abs(IQh(20, :, :, 1)))')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % SVD filter on the hole portion separately
    [PPh, EVsh, V_sorth] = getSVs2D(IQh);

    [IQfh] = applySVs2D(IQh, PPh, EVsh, V_sorth, 1, 199);
    disp('SVs applied in the hole')

%     IQhrf1 = imresize3(IQfh(:, :, :, 1), [size(IQfh(:, :, :, 1), 1) * imgRefinementFactor(1), size(IQfh(:, :, :, 1), 2) * imgRefinementFactor(2), size(IQfh(:, :, :, 1), 3) * imgRefinementFactor(3)]);
%     testXC = normxcorr3(abs(refPSF), abs(IQhrf1));
%     volumeViewer(testXC)
% %     testXCscaled = testXC;
% %     testXCscaled(testXCscaled > 0) = testXCscaled(testXCscaled > 0) .^ 0.3;
% %     volumeViewer(testXCscaled)
% 
%     XCt = testXC; % XC Thresholded
%     XCThresholdAdaptive = XCThresholdFactor * max(XCt, [], 'all');
% %     XCt(XCt < XCThreshold) = 0;
%     XCt(XCt < XCThresholdAdaptive) = 0;
%     centers = imregionalmax(XCt, 6);

    [xp, yp, zp, nf] = size(IQfh);
    range{4} = int16(1:nf); 
    [centers, ~, ~, XCThresholdAdaptive] = localizeBubbles3D(IQfh, refPSF, range, imgRefinementFactor, XCThresholdFactor);
    disp('Localization in the hole done')
%     volumeViewer(abs(IQfh(:, :, :, 1)))
%     figure; imagesc(squeeze(abs(IQfh(20, :, :, 1)))')
%     figure; imagesc(squeeze(abs(IQfh(20, :, :, 1)))' .^ 0.3)

%     figure; imagesc(squeeze(sum(centers(:, :, :, 1), 1)))
%     figure; imagesc(squeeze(sum(centers(:, :, :, 2), 1)))
    
    % make 2D video for testing
    for f = 1:size(centers, 4)
        
        v1f_ysum = sum(centers(:, :, :, f), 1);
        figure(vf_ysum)
        imagesc(squeeze(v1f_ysum)')
        
        cv = getframe(vf_ysum);     % get the current volume
        rgb = frame2im(cv);      % convert the frame to rgb data
    
        writeVideo(vo_ysum, rgb);


        v1f_xsum = sum(centers(:, :, :, f), 2);
        figure(vf_xsum)
        imagesc(squeeze(v1f_xsum)')
        
        cv = getframe(vf_xsum);     % get the current volume
        rgb = frame2im(cv);      % convert the frame to rgb data
    
        writeVideo(vo_xsum, rgb);

        v1f_zsum = sum(centers(:, :, :, f), 3);
        figure(vf_zsum)
        imagesc(squeeze(v1f_zsum)')
        
        cv = getframe(vf_zsum);     % get the current volume
        rgb = frame2im(cv);      % convert the frame to rgb data
    
        writeVideo(vo_zsum, rgb);
    end

    % Framewise difference for extracting the bubble signal
%     IQf = diff(IQ, 1, 4);

    % set frame range
%     if filenum == 1
%         [xp, yp, zp, nf] = size(IQf);
%         range{4} = int16(1:nf); 
%     end
    
    % Denoise
%     IQf_dn = zeros(size(IQf));
%     parfor f = 1:size(IQf, 4) % Go through each frame and denoise
% %         tic
%         IQf_scaled_temp = squeeze(abs(IQf(:, :, :, f)) ./ max(abs(IQf(:, :, :, f)), [], 'all'));
%         IQf_dn(:, :, :, f) = NLMF(IQf_scaled_temp);
% %         toc
%     end

%     IQf_dn_rfn = 
%     parfor f = 1:size(IQf_dn)
%         IQf_dn_rfn = imresize3(IQf_dn, [size(IQf_dn, 1) * imgRefinementFactor(1), size(IQf_dn, 2) * imgRefinementFactor(2), size(IQf_dn, 3) * imgRefinementFactor(3)]);
%     end
%     testXC = normxcorr3(abs(refPSF), IQf_dn_rfn);
%     toc

%     clear PP EVs V_sort

%     save([savepath, 'Filtered-Data-', num2str(filenum)], 'IQr', 'PP', 'EVs', 'V_sort', 'IQf', "-v6")
%     [centers, ~, ~, XCThresholdAdaptive] = localizeBubbles3D(IQf, refPSF, range, imgRefinementFactor, XCThresholdFactor);
%     [coords, img_size, XCThresholdsAdaptive] = localizeBubbles3D_framewise(IQf, refPSF, range, imgRefinementFactor, XCThresholdFactor);
%     save([savepath, 'IQf-', num2str(filenum)], 'IQf', "-v6")

%     save([savepath, 'dataproc-', num2str(filenum)], 'IQf', 'centroidCoordinates', "-v6")
%     savefast([savepath, 'centers-', num2str(filenum)], 'centers', 'XCThresholdAdaptive')
    savefast([savepath, 'centers-', num2str(filenum)], 'centers', 'XCThresholdAdaptive', 'IQf_zsum_mask_rm')
%     savefast([savepath, 'coords-', num2str(filenum)], 'coords', 'img_size', 'XCThresholdsAdaptive')

    disp(strcat("Centroid finding done: file ", num2str(filenum)))
    toc
end
save([savepath, 'proc_params.mat'], 'sv_threshold_lower', 'sv_threshold_upper', 'PSF', 'range', 'imgRefinementFactor', 'XCThresholdFactor', 'xpix_spacing', 'ypix_spacing', 'zpix_spacing')
%%
close(vo_ysum)
close(vo_xsum)
close(vo_zsum)