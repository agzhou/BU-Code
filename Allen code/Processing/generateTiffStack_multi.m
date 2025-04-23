%% Make a tiff stack
% Input data is assumed to be a cell array with volumes of (y, x, z)
% Optional inputs: 1. actual image/physical size [y, x, z] in the same physical units
%                  2. colormap (string)
%                  3. MIP window size

function generateTiffStack_multi(volumeData, varargin)

    % volumeData = volumeData ./ max(volumeData, [], 'all'); % Normalize intensities to be between 0 - 1
    showColorbar = false;
    
    mws = 1; % default MIP window size is 1 (no MIP)
    if nargin > 1
        actualSize = varargin{1};
        if nargin > 2
            cmap = varargin{2};
            if nargin > 3
                mws = varargin{3}; % MIP window size
            end
        end
    end
    savepath = uigetdir('D:\Allen\Data\', 'Select the save path');
    savepath = [savepath, '\'];

    numVolumes = length(volumeData); % get the # of volumes from the input
    cr = zeros(numVolumes, 2); % color range
    nyp = zeros(numVolumes, 1); % number of y pixels
    nxp = zeros(numVolumes, 1); % x
    nzp = zeros(numVolumes, 1); % z

    for vi = 1:numVolumes % volume index
        cr(vi, :) = [0, max(volumeData{vi}, [], 'all')];
        nyp(vi) = size(volumeData{vi}, 1);
        nxp(vi) = size(volumeData{vi}, 2);
        nzp(vi) = size(volumeData{vi}, 3);
    end


    %% Generate xz planes
    xz_stack = Tiff([savepath, 'xz_stack.tif'], 'w');
%     xz_tagstruct.ImageLength = nzp;
%     xz_tagstruct.ImageWidth = nxp;
    xz_tagstruct.Photometric = Tiff.Photometric.RGB;
    xz_tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    xz_tagstruct.BitsPerSample = 8;
    xz_tagstruct.SamplesPerPixel = 3;
    xz_tagstruct.Software = 'MATLAB';
    setTag(xz_stack, xz_tagstruct) % set the tags

    tf = figure;
    % set(tf, 'Position', [100, 100, round(nzp * numVolumes * actualSize(3)/actualSize(1)), nxp]); % Define figure size early on
    %%%% this line assumes the same size across each volume in volumeData ****
    set(tf, 'Position', [0, 0, nzp(1) * numVolumes, nxp(1)] .* 5); % Define figure size early on
    tiledlayout(1, numVolumes)
    % subplot(1, numVolumes, 1)
    colormap(cmap)

    % adjust size of the image
%     if exist('actualSize', 'var')
%         hwRatio_xz = actualSize(3) / actualSize(2); % height to width ratio
%     else
%         hwRatio_xz = 1;
%     end
    
        
    for y = 1:size(volumeData{vi}, 1) - mws + 1
        for vi = 1:numVolumes
            nexttile(vi)
            planeTemp = squeeze(max(volumeData{vi}(y:y + mws - 1, :, :), [], 1));
            figure(tf)
            
            % subplot(1, numVolumes, vi, 'Position', [0, 0, (vi - 1)/numVolumes, 0])
            
            imagesc(squeeze(planeTemp)')
            if showColorbar
                colorbar
            end
            clim(cr(vi, :))
            cv = getframe(tf);
            rgb = frame2im(cv);      % convert the frame to rgb data
            if y == 1
                xz_tagstruct.ImageLength = size(rgb, 1);
                xz_tagstruct.ImageWidth = size(rgb, 2);
                setTag(xz_stack, xz_tagstruct) % set the tags
                
    %             tf.Position(4) = ceil(tf.Position(3) * hwRatio_xz);
                
            else
                writeDirectory(xz_stack)
                setTag(xz_stack, xz_tagstruct)
            end
        end

        write(xz_stack, rgb)
    end
    close(xz_stack)

%     %% yz planes
%     yz_stack = Tiff([savepath, 'yz_stack.tif'], 'w');
% %     yz_tagstruct.ImageLength = nzp;
% %     yz_tagstruct.ImageWidth = nxp;
%     yz_tagstruct.Photometric = Tiff.Photometric.RGB;
%     yz_tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
%     yz_tagstruct.BitsPerSample = 8;
%     yz_tagstruct.SamplesPerPixel = 3;
%     yz_tagstruct.Software = 'MATLAB';
%     setTag(yz_stack, yz_tagstruct) % set the tags
% 
%     tf = figure;
%     subplot(1, numVolumes, 1)
%     colormap(cmap)
% 
%     % adjust size of the image
% %     if exist('actualSize', 'var')
% %         hwRatio_yz = actualSize(3) / actualSize(1); % height to width ratio
% %     else
% %         hwRatio_yz = 1;
% %     end
%     for vi = 1:numVolumes
%         for x = 1:size(volumeData{vi}, 2) - mws + 1
%             planeTemp = squeeze(max(volumeData{vi}(:, x:x + mws - 1, :), [], 2));
% 
%             subplot(1, numVolumes, vi)
% 
%             imagesc(squeeze(planeTemp)')
%             if showColorbar
%                 colorbar
%             end
%             clim(cr(vi, :))
%             cv = getframe(tf);
%             rgb = frame2im(cv);      % convert the frame to rgb data
%             if x == 1
%                 yz_tagstruct.ImageLength = size(rgb, 1);
%                 yz_tagstruct.ImageWidth = size(rgb, 2);
%                 setTag(yz_stack, yz_tagstruct) % set the tags
% 
%     %             tf.Position(4) = ceil(tf.Position(3) * hwRatio_yz);
% 
%             else
%                 writeDirectory(yz_stack)
%                 setTag(yz_stack, yz_tagstruct)
%             end
% 
%             write(yz_stack, rgb)
%         end
%     end
%     close(yz_stack)
% 
%     %% xy planes
%     xy_stack = Tiff([savepath, 'xy_stack.tif'], 'w');
% %     xy_tagstruct.ImageLength = nzp;
% %     xy_tagstruct.ImageWidth = nxp;
%     xy_tagstruct.Photometric = Tiff.Photometric.RGB;
%     xy_tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
%     xy_tagstruct.BitsPerSample = 8;
%     xy_tagstruct.SamplesPerPixel = 3;
%     xy_tagstruct.Software = 'MATLAB';
%     setTag(xy_stack, xy_tagstruct) % set the tags
% 
%     tf = figure;
%     subplot(1, numVolumes, 1)
%     colormap(cmap)
% 
%     % adjust size of the image
% %     if exist('actualSize', 'var')
% %         hwRatio_xy = actualSize(3) / actualSize(2); % height to width ratio
% %     else
% %         hwRatio_xy = 1;
% %     end
%     for vi = 1:numVolumes
%         for z = 1:size(volumeData{vi}, 3) - mws + 1
%             planeTemp = squeeze(max(volumeData{vi}(:, :, z:z + mws - 1), [], 3));
% 
%             subplot(1, numVolumes, vi)
% 
%             imagesc(squeeze(planeTemp)')
%             if showColorbar
%                 colorbar
%             end
%             clim(cr(vi, :))
%             cv = getframe(tf);
%             rgb = frame2im(cv);      % convert the frame to rgb data
%             if z == 1
%                 xy_tagstruct.ImageLength = size(rgb, 1);
%                 xy_tagstruct.ImageWidth = size(rgb, 2);
%                 setTag(xy_stack, xy_tagstruct) % set the tags
% 
%     %             tf.Position(4) = ceil(tf.Position(3) * hwRatio_xy);
% 
%             else
%                 writeDirectory(xy_stack)
%                 setTag(xy_stack, xy_tagstruct)
%             end
% 
%             write(xy_stack, rgb)
%         end
%     end
%     close(xy_stack)

    
end