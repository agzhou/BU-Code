%% Description: look at the motion within an ultrasound acquisition

%% Load an example frame
load('E:\Allen BME-BOAS-27 Data Backup\AZ01 fUS\06-19-2025 L22-14v\run 8 awake wooden stick stim max TGC yes BPF 5 angles from -5 to 5 deg 3 half cycle TW\IQ Data - Verasonics recon\IQ-5-5-5000-1000-1-8.mat')

%%
IQ = squeeze(IData + 1i .* QData);

%%
IQphase = angle(IQ);
%%
figure; imagesc(IQphase(:, :, 1))
figure; imagesc(IQphase(:, :, 2))
figure; imagesc(IQphase(:, :, 3))
figure; imagesc(IQphase(:, :, 4))
figure; imagesc(IQphase(:, :, 5))

% generateTiffStack_acrossframes(IQphase, [8.8, 8.8, 8], 'parula')

%% Use block matching to compute the motion vector of the phase
blockSize = [11, 11];
maxDisp = [5, 5];
blkMatcher = vision.BlockMatcher("ReferenceFrameSource", "Input port", ...
    'OutputValue', 'Horizontal and vertical components in complex form', ...
    'BlockSize', blockSize, 'MaximumDisplacement', maxDisp);
%%
test = IQphase(:, :, [1, 500, 1000]);
figure; imagesc(test(:, :, 1))
figure; imagesc(test(:, :, 2))
V = blkMatcher(test(:, :, 1), test(:, :, 2));
% figure; imagesc(V)
% V = blkMatcher(test);
% V = blkMatcher(IQphase);
img1 = test(:, :, 1);
img2 = test(:, :, 2);
[X,Y] = meshgrid(1:blockSize(1):size(img1, 2), 1:blockSize(2):size(img1, 1));
halphablend = vision.AlphaBlender;
% img12 = halphablend(img2, img1);
% imshow(img12)
% hold on
figure;
quiver(X(:),Y(:),real(V(:)),imag(V(:)),0)
% hold off
%% Block matching for every successive pair of frames
PMV = []; % Phase "motion vector"
for fi = 1:size(IQphase, 3) - 1 % Go through frames
    PMV(:, :, fi) = blkMatcher(IQphase(:, :, fi), IQphase(:, :, fi + 1));
end

%% Block matching for larger gaps of frames
PMV = []; % Phase "motion vector"
PMV_frame_interval = 100;
for fi = 1:size(IQphase, 3)/PMV_frame_interval % Go through frame gaps
    PMV(:, :, fi) = blkMatcher(IQphase(:, :, fi), IQphase(:, :, fi + 1));
end

PMV_diff = diff(IQphase, 1, 3);
lst = find(PMV_diff>3); PMV_diff(lst) = PMV_diff(lst) - 2*3.14159;
lst = find(PMV_diff<-3); PMV_diff(lst) = PMV_diff(lst) + 2*3.14159;
PMV_diff = convn( PMV_diff, ones(1,1,5), 'same');

%% ROI average of phase differences

% phase_change = IQphase - IQphase(:, :, 1);
% ROI = [{30:50}; {150:200}]; % z, x ranges
% PC_ROI = phase_change(ROI{1}, ROI{2}, :);
% PC_ROI_avg = squeeze(mean(mean(PC_ROI, 1), 2));
% 
% figure; plot(PC_ROI_avg)

IQphase_diff = diff(IQphase, 1, 3);
IQphase_diff_uw = PMV_diff; % Unwrapped IQ phase diff
% Ac count for phase wrapping (thanks David)
lst = find(IQphase_diff_uw > 3); IQphase_diff_uw(lst) = IQphase_diff_uw(lst) - 2*pi;
lst = find(IQphase_diff_uw < -3); IQphase_diff_uw(lst) = IQphase_diff_uw(lst) + 2*pi;
IQphase_diff_uw_smoothed = convn( IQphase_diff_uw, ones(1, 1, 5), 'same');

ROI = [{30:50}; {150:200}]; % z, x ranges
IQphase_diff_ROI = squeeze(mean(mean(IQphase_diff(ROI{1}, ROI{2}, :), 1), 2));
IQphase_diff_uw_ROI = squeeze(mean(mean(IQphase_diff_uw(ROI{1}, ROI{2}, :), 1), 2));

IQphase_diff_uw_smoothed_ROI = squeeze(mean(mean(IQphase_diff_uw_smoothed(ROI{1}, ROI{2}, :), 1), 2));
figure; 
frame_times = (1:size(IQphase_diff_uw_smoothed, 3)) ./ P.frameRate;
plot( frame_times, IQphase_diff_uw_ROI ); xlabel('Time [s]'); ylabel('Framewise phase difference [radians]'); title('IQ Phase framewise difference (unwrapped): ROI from z = 30:50, x = 150:200')
hold on
plot( IQphase_diff_uw_smoothed_ROI )
hold off

% avg_IQphase_change = trapz( (1:size(IQphase_diff_uw_smoothed, 3)) ./ P.frameRate, IQphase_diff_uw_smoothed_ROI);
% Since we didn't do a diff and account for the actual time between data
% points, we don't need to do that in the integration
avg_IQphase_change = trapz(IQphase_diff_uw_ROI);
% avg_IQphase_change = trapz(IQphase_diff_uw_smoothed_ROI);

avg_IQphase_cumsum = cumsum(IQphase_diff_uw_ROI); % Progressive integral
figure; plot(frame_times, avg_IQphase_cumsum); xlabel('Time [s]'); ylabel('Cumulative phase change [radians]'); title('Cumulative IQ phase change (unwrapped): ROI from z = 30:50, x = 150:200')




PMV_diff_ROI = PMV_diff(ROI{1}, ROI{2}, :);
PMV_diff_ROI_avg = squeeze(mean(mean(PMV_diff_ROI, 1), 2));

figure; plot(PMV_diff_ROI_avg)
figure; plot(convn( PMV_diff_ROI_avg, ones(1,1,5), 'same'));

figure; plot(convn( squeeze(mean(mean(PMV_diff(40:50, 40:50, :), 1), 2)), ones(1,1,5), 'same') )
%% Phase contrast

dIQ = IQ(:, :, 1:end-1) .* conj(IQ(:, :, 2:end));
aIQ = angle(dIQ);
% generateTiff(aIQ, [-pi, pi]);
generateTiff(aIQ, [-1, 1]);
%% Another phase contrast

dIQ = IQ(:, :, 1) .* conj(IQ(:, :, :));
aIQ = angle(dIQ);
% generateTiff(aIQ, [-pi, pi]);
generateTiff(aIQ, [-1, 1]);

%% Another phase contrast

dIQ = IQ(:, :, end/2) .* conj(IQ(:, :, :));
aIQ = angle(dIQ);
% generateTiff(aIQ, [-pi, pi]);
generateTiff(aIQ, [-1, 1]);
%% Use a moving mean
dIQ_mm = movmean(dIQ, 5, 3);
aIQ_mm = angle(dIQ_mm);
% generateTiff(aIQ_mm, [-pi, pi]);
generateTiff(aIQ_mm, [-1, 1]);

%% Look at the "phase fronts" with a gradient 

% % Testing
% img = IQphase(:, :, 1);
% figure; imagesc(img)
% 
% [Gh, Gv] = imgradient(img); % Horizontal and vertical gradients
% figure; imagesc(Gh)
% figure; imagesc(Gv)
% 
% [Gmag, Gdir] = imgradient(Gh, Gv);

IQphasefronts = zeros(size(IQphase));
for fi = 1:size(IQphase, 3) % Go through frames
    [Gh, Gv] = imgradient(IQphase(:, :, fi)); % Horizontal and vertical gradients
    [Gmag, Gdir] = imgradient(Gh, Gv); % Calculate the gradient magnitude and direction
    IQphasefronts(:, :, fi) = Gdir;

end

%% Make a video for the phase fronts
generateTiff(IQphasefronts);

%% Quiver plot video for the Phase "Motion Vector"
% savepath = uigetdir('D:\Allen\Data\', 'Select the save path');
% savepath = [savepath, '\'];
% 
% PMV_tiff = Tiff([savepath, 'PMV.tif'], 'w');
%     %     xz_tagstruct.ImageLength = nzp;
%     %     xz_tagstruct.ImageWidth = nxp;
% PMV_tagstruct.Photometric = Tiff.Photometric.RGB;
% PMV_tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
% PMV_tagstruct.BitsPerSample = 8;
% PMV_tagstruct.SamplesPerPixel = 3;
% PMV_tagstruct.Software = 'MATLAB';
% setTag(PMV_tiff, PMV_tagstruct) % set the tags
% 
% % adjust size of the image
% %     if exist('actualSize', 'var')
% %         hwRatio_xz = actualSize(3) / actualSize(2); % height to width ratio
% %     else
% %         hwRatio_xz = 1;
% %     end
% 
% PMV_fh = figure;
% [X, Y] = meshgrid(1:blockSize(1):size(IQphase, 2), 1:blockSize(2):size(IQphase, 1));
% 
% for fi = 1:size(PMV, 3) % Go through pairs of frames
%     temp_PMV = PMV(:, :, fi);
%     quiver(X(:), Y(:), real(temp_PMV(:)), imag(temp_PMV(:)), 0)
%     title("Phase motion vector: Frame " + num2str(fi) + " to " + num2str(fi + 1))
% 
%     cf = getframe(PMV_fh);
%     rgb = frame2im(cf);
% 
%     if fi == 1
%         PMV_tagstruct.ImageLength = size(rgb, 1);
%         PMV_tagstruct.ImageWidth = size(rgb, 2);
%         setTag(PMV_tiff, PMV_tagstruct) % set the tags
% %             tf.Position(4) = ceil(tf.Position(3) * hwRatio_xz);
%     else
%         writeDirectory(PMV_tiff)
%         setTag(PMV_tiff, PMV_tagstruct)
%     end
%     write(PMV_tiff, rgb)
% end
% close(PMV_tiff)
% 
% savepath = uigetdir('D:\Allen\Data\', 'Select the save path');
%     savepath = [savepath, '\'];
%     
%     PMV_tiff = Tiff([savepath, 'PMV_diff.tif'], 'w');
%         %     xz_tagstruct.ImageLength = nzp;
%         %     xz_tagstruct.ImageWidth = nxp;
%     PMV_tagstruct.Photometric = Tiff.Photometric.RGB;
%     PMV_tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
%     PMV_tagstruct.BitsPerSample = 8;
%     PMV_tagstruct.SamplesPerPixel = 3;
%     PMV_tagstruct.Software = 'MATLAB';
%     PMV_tagstruct.Compression = Tiff.Compression.Deflate;
%     setTag(PMV_tiff, PMV_tagstruct) % set the tags
%     
%     % adjust size of the image
%     %     if exist('actualSize', 'var')
%     %         hwRatio_xz = actualSize(3) / actualSize(2); % height to width ratio
%     %     else
%     %         hwRatio_xz = 1;
%     %     end
%     
%     PMV_fh = figure;
%     [X, Y] = meshgrid(1:blockSize(1):size(IQphase, 2), 1:blockSize(2):size(IQphase, 1));
%     
%     for fi = 1:size(PMV_diff, 3) % Go through pairs of frames
%         temp_PMV = PMV_diff(:, :, fi);
%         imagesc(temp_PMV)    
%         clim([-0.1, 0.1])
%         colormap jet
%         title("Phase motion vector difference: Frame " + num2str(fi) + " to " + num2str(fi + 1))
%     
%         cf = getframe(PMV_fh);
%         rgb = frame2im(cf);
%     
%         if fi == 1
%             PMV_tagstruct.ImageLength = size(rgb, 1);
%             PMV_tagstruct.ImageWidth = size(rgb, 2);
%             setTag(PMV_tiff, PMV_tagstruct) % set the tags
%     %             tf.Position(4) = ceil(tf.Position(3) * hwRatio_xz);
%         else
%             writeDirectory(PMV_tiff)
%             setTag(PMV_tiff, PMV_tagstruct)
%         end
%         write(PMV_tiff, rgb)
%     end
%     close(PMV_tiff)

%% IQ phase video across all frames
% Choose and load an IQ mat file (from the acquisition)
[IQ_filename, IQ_pathname, ~] = uigetfile('*.mat', 'Select an IQ file', 'E:\Allen BME-BOAS-27 Data Backup\');
%     load([IQ_pathname, IQ_filename], 'P')

% Remove the last (file) number from the IQ filename to define the general
% filename structure
IQ_filename_split = strsplit(IQ_filename, "-");
IQ_filename_join = strjoin(IQ_filename_split(1:end - 1), "-");
IQ_filename_join = IQ_filename_join + "-"; % add the dash at the end
IQfilenameStructure = char(IQ_filename_join);

parameterPrompt = {'Start file number', 'End file number'};
parameterDefaults = {'1', '',};
parameterUserInput = inputdlg(parameterPrompt, 'Input Parameters', 1, parameterDefaults);

startFile = str2double(parameterUserInput{1});
endFile = str2double(parameterUserInput{2});
numFiles = endFile - startFile + 1;

%
savepath = uigetdir('D:\Allen\Data\', 'Select the save path');
savepath = [savepath, '\'];

IQphase_fh = figure;

IQphase_tiff = Tiff([savepath, 'IQphase.tif'], 'w');
    %     xz_tagstruct.ImageLength = nzp;
    %     xz_tagstruct.ImageWidth = nxp;
IQphase_tagstruct.Photometric = Tiff.Photometric.RGB;
IQphase_tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
IQphase_tagstruct.BitsPerSample = 8;
IQphase_tagstruct.SamplesPerPixel = 3;
IQphase_tagstruct.Software = 'MATLAB';
IQphase_tagstruct.Compression = Tiff.Compression.Deflate;
setTag(IQphase_tiff, IQphase_tagstruct) % set the tags

% adjust size of the image
%     if exist('actualSize', 'var')
%         hwRatio_xz = actualSize(3) / actualSize(2); % height to width ratio
%     else
%         hwRatio_xz = 1;
%     end
initTag = false;
for sfi = startFile:endFile % Go through superframes and load each file

% for sfi = 8
    load([IQ_pathname, IQfilenameStructure, num2str(sfi)])
    IQ_temp = squeeze(IData + 1i .* QData);
    IQphase_temp = angle(IQ_temp); % Calculate the IQ phase

    nf = size(IQphase_temp, 3);
    for fi = 1:nf % Go through frames
        IQphase_temp_fi = IQphase_temp(:, :, fi);
        
        % Only do the frame every X frames
        if mod(fi - 1, 10) == 0
            imagesc(IQphase_temp_fi)    
            clim([-pi, pi])
            colorbar
    %         colormap jet
            title("IQ phase: Frame " + num2str((sfi - 1) * nf + fi))
        
            cf = getframe(IQphase_fh);
            rgb = frame2im(cf);
        
%             if sfi == 1 & fi == 1
            if fi == 1 & initTag == false
                IQphase_tagstruct.ImageLength = size(rgb, 1);
                IQphase_tagstruct.ImageWidth = size(rgb, 2);
                setTag(IQphase_tiff, IQphase_tagstruct) % set the tags
                initTag = true;
        %             tf.Position(4) = ceil(tf.Position(3) * hwRatio_xz);
            else
                writeDirectory(IQphase_tiff)
                setTag(IQphase_tiff, IQphase_tagstruct)
            end
            write(IQphase_tiff, rgb)
        end
    end
end
close(IQphase_tiff)

%% Helper functions

% Tiff across frames
function generateTiff(data, varargin)
    if nargin > 1
        clims = varargin{1}; % Optonal input for the colorbar limits
    end
    savepath = uigetdir('D:\Allen\Data\', 'Select the save path');
    savepath = [savepath, '\'];
    
    data_tiff = Tiff([savepath, 'data.tif'], 'w');
    data_tagstruct.Photometric = Tiff.Photometric.RGB;
    data_tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    data_tagstruct.BitsPerSample = 8;
    data_tagstruct.SamplesPerPixel = 3;
    data_tagstruct.Software = 'MATLAB';
    data_tagstruct.Compression = Tiff.Compression.Deflate;
    setTag(data_tiff, data_tagstruct) % set the tags
    
    % adjust size of the image
    %     if exist('actualSize', 'var')
    %         hwRatio_xz = actualSize(3) / actualSize(2); % height to width ratio
    %     else
    %         hwRatio_xz = 1;
    %     end
    
    fh = figure;
%     [X, Y] = meshgrid(1:blockSize(1):size(IQphase, 2), 1:blockSize(2):size(IQphase, 1));
    
    for fi = 1:size(data, 3) % Go through frames
        temp_data = data(:, :, fi);
        imagesc(temp_data)    
%         clim([-0.1, 0.1])
%         colormap jet
        title("Frame " + num2str(fi))
        if exist('clims', 'var')
            clim(clims)
        end
        colorbar
    
        cf = getframe(fh);
        rgb = frame2im(cf);
    
        if fi == 1
            data_tagstruct.ImageLength = size(rgb, 1);
            data_tagstruct.ImageWidth = size(rgb, 2);
            setTag(data_tiff, data_tagstruct) % set the tags
    %             tf.Position(4) = ceil(tf.Position(3) * hwRatio_xz);
        else
            writeDirectory(data_tiff)
            setTag(data_tiff, data_tagstruct)
        end
        write(data_tiff, rgb)
    end
    close(data_tiff)
end