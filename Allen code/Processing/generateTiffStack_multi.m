%% Make tiff stacks of volumetric data
% Input data is assumed to be a cell array with volumes of (y, x, z)
% Optional inputs: 1. actual image/physical size [y, x, z] in the same physical units
%                  2. colormap (string) - choose from Matlab's options
%                  3. MIP window size

% e.g., generateTiffStack_multi(volumeData, [8.8, 8.8, 8], 'hot', 10)

% Output: xz, yz, and xy fly throughs which are saved to a user-specified
%         directory

% Possible issues:
%                  - The stuff I added with aspect ratios might break if an
%                    actualSize is not specified

function generateTiffStack_multi(volumeData, varargin)

    % volumeData = volumeData ./ max(volumeData, [], 'all'); % Normalize intensities to be between 0 - 1
    showColorbar = true;
    scale = 5;
    
    mws = 1; % default MIP window size is 1 (no MIP)
    if nargin > 1
        actualSize = varargin{1}; % Store the actual image region size if that's an input
        if nargin > 2
            cmap = varargin{2};
            if nargin > 3
                mws = varargin{3}; % MIP window size
            end
                if nargin > 4
                    cr = varargin{4}; % colorbar limits/range
                    crUserInput = true;
                else
                    crUserInput = false;
                end
        end
    end
    savepath = uigetdir('D:\Allen\Data\', 'Select the save path');
    savepath = [savepath, '\'];

    numVolumes = length(volumeData); % get the # of volumes from the input
    if ~crUserInput % if a color range is specified in the input
        cr = zeros(numVolumes, 2); % color range
    else
        cr = repmat(cr, numVolumes, 1);
    end
    nyp = zeros(numVolumes, 1); % number of y pixels
    nxp = zeros(numVolumes, 1); % x
    nzp = zeros(numVolumes, 1); % z

    for vi = 1:numVolumes % volume index
        if ~crUserInput % if a color range is specified in the input
            cr(vi, :) = [0, max(volumeData{vi}, [], 'all')]; % Set the color range from 0 to the max of the volumeData for each volume
        end
        nyp(vi) = size(volumeData{vi}, 1);
        nxp(vi) = size(volumeData{vi}, 2);
        nzp(vi) = size(volumeData{vi}, 3);
    end

    %% Generate xz planes
%     xz_stack = Tiff([savepath, 'xz_stack.tif'], 'w');
% %     xz_tagstruct.ImageLength = nzp;
% %     xz_tagstruct.ImageWidth = nxp;
%     xz_tagstruct.Photometric = Tiff.Photometric.RGB;
%     xz_tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
%     xz_tagstruct.BitsPerSample = 8;
%     xz_tagstruct.SamplesPerPixel = 3;
%     xz_tagstruct.Software = 'MATLAB';
%     xz_tagstruct.Compression = Tiff.Compression.Deflate;
%     setTag(xz_stack, xz_tagstruct) % set the tags
% 
%     tf = figure;
%     %%%% this line below assumes the same size across each volume in volumeData ****
%     % set(tf, 'Position', [0, 0, round(nzp(1) * numVolumes * actualSize(3)/actualSize(1)), nxp(1)] .* scale); % Define figure size
%     set(tf, 'Position', round([0, 0, nxp(1) * numVolumes, nzp(1) * actualSize(3)/actualSize(2)] .* scale)); % Define figure size
%     % set(tf, 'Position', [0, 0, nzp(1) * numVolumes, nxp(1)] .* scale); % Define figure size
%     tiledlayout(1, numVolumes)
%     % subplot(1, numVolumes, 1)
%     colormap(cmap)
% 
%     % adjust size of the image
%     if exist('actualSize', 'var')
%         hwRatio_xz = actualSize(3) / actualSize(2); % height to width ratio
%     else
%         hwRatio_xz = 1;
%     end
% 
%     for y = 1:size(volumeData{vi}, 1) - mws + 1
%         for vi = 1:numVolumes
%             nexttile(vi) % cycle between tiles
%             planeTemp = squeeze(max(volumeData{vi}(y:y + mws - 1, :, :), [], 1));
%             figure(tf)
% 
%             % subplot(1, numVolumes, vi, 'Position', [0, 0, (vi - 1)/numVolumes, 0])
% 
%             imagesc(squeeze(planeTemp)')
%             axis tight
%             ax = gca;
%             % set(ax, 'Color', 'k') % set background to black
%             ax.PlotBoxAspectRatio = [1, hwRatio_xz, 1];
% 
%             if showColorbar
%                 colorbar
%             end
%             clim(cr(vi, :)) % set the color range
% 
%             if y == 1 && vi == 1
%                 % Set the tick labels, modifying what the default is
%                 default_xz_xticks = xticks;
%                 if mod(length(default_xz_xticks), 2) == 0 % Modify if the # of ticks is even, so we can get a tick at 0
%                     new_xz_xticks = [1, default_xz_xticks];
%                 else
%                     new_xz_xticks = default_xz_xticks;
%                 end
% 
%                 lateral_range_mm = linspace(-actualSize(2)/2, actualSize(2)/2, length(new_xz_xticks));
% 
%                 new_xz_xticklabels = {};
%                 for i = 1:length(new_xz_xticks)
%                     new_xz_xticklabels{i} = num2str(lateral_range_mm(i));
%                 end
% 
%                 nzp_xz = nzp(1); % # z pixels in the xz image
%                 new_xz_ytick_interval_mm = 1; % Set the interval between y ticks
%                 num_xz_yticks = actualSize(3)/new_xz_ytick_interval_mm; % number of y ticks
%                 new_xz_yticks = round(0:nzp_xz / num_xz_yticks:nzp_xz); new_xz_yticks(1) = 1;
%                 axial_range_mm = 0:new_xz_ytick_interval_mm:actualSize(3);
%                 new_xz_yticklabels = {};
%                 for i = 1:length(new_xz_yticks)
%                     new_xz_yticklabels{i} = num2str(axial_range_mm(i));
%                 end
%             end
% 
%             % Set labels
%             % figure(figHandle)
%             xlabel('x [mm]')
%             ylabel('z [mm]')
%             xticks(new_xz_xticks)
%             xticklabels(new_xz_xticklabels)
%             yticks(new_xz_yticks)
%             yticklabels(new_xz_yticklabels)
%             set(gca, 'TickDir', 'out')
%             set(gca, 'box', 'off')
%         end
% 
%         cv = getframe(tf); % Get the frame once all the tiles are populated for that MIP
%         rgb = frame2im(cv);      % convert the frame to rgb data
% 
%         if y == 1 % Some stuff to get the tif correct
%             xz_tagstruct.ImageLength = size(rgb, 1);
%             xz_tagstruct.ImageWidth = size(rgb, 2);
%             setTag(xz_stack, xz_tagstruct) % set the tags
% 
% %             tf.Position(4) = ceil(tf.Position(3) * hwRatio_xz);
% 
%         else
%             writeDirectory(xz_stack)
%             setTag(xz_stack, xz_tagstruct)
%         end
%         write(xz_stack, rgb)
%     end
%     close(xz_stack)

    %% yz planes
    yz_stack = Tiff([savepath, 'yz_stack.tif'], 'w');
%     yz_tagstruct.ImageLength = nzp;
%     yz_tagstruct.ImageWidth = nxp;
    yz_tagstruct.Photometric = Tiff.Photometric.RGB;
    yz_tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    yz_tagstruct.BitsPerSample = 8;
    yz_tagstruct.SamplesPerPixel = 3;
    yz_tagstruct.Software = 'MATLAB';
    yz_tagstruct.Compression = Tiff.Compression.Deflate;
    setTag(yz_stack, yz_tagstruct) % set the tags

    tf = figure;
    set(tf, 'Position', round([0, 0, nyp(1) * numVolumes, nzp(1) * actualSize(3)/actualSize(1)] .* scale)); % Define figure size
%     set(tf, 'Position', [0, 0, round(nzp(1) * numVolumes * actualSize(3)/actualSize(2)), nyp(1)] .* scale); % Define figure size
%     set(tf, 'Position', [0, 0, nzp(1) * numVolumes, nyp(1)] .* scale); % Define figure size
    tiledlayout(1, numVolumes)
    colormap(cmap)

    % adjust size of the image
    if exist('actualSize', 'var')
        hwRatio_yz = actualSize(3) / actualSize(1); % height to width ratio
    else
        hwRatio_yz = 1;
    end

    for x = 1:size(volumeData{vi}, 2) - mws + 1
        for vi = 1:numVolumes
            nexttile(vi) % cycle between tiles
            planeTemp = squeeze(max(volumeData{vi}(:, x:x + mws - 1, :), [], 2));

            imagesc(squeeze(planeTemp)')
            axis tight
            ax = gca;
            % set(ax, 'Color', 'k') % set background to black
            ax.PlotBoxAspectRatio = [1, hwRatio_yz, 1];
            if showColorbar
                colorbar
            end
            clim(cr(vi, :))

            if x == 1 && vi == 1
                % Set the tick labels, modifying what the default is
                default_yz_xticks = xticks;
                if mod(length(default_yz_xticks), 2) == 0 % Modify if the # of ticks is even, so we can get a tick at 0
                    new_yz_xticks = [1, default_yz_xticks];
                else
                    new_yz_xticks = default_yz_xticks;
                end

                lateral_range_mm = linspace(-actualSize(1)/2, actualSize(1)/2, length(new_yz_xticks));

                new_yz_xticklabels = {};
                for i = 1:length(new_yz_xticks)
                    new_yz_xticklabels{i} = num2str(lateral_range_mm(i));
                end

                nzp_yz = nzp(1); % # z pixels in the xz image
                new_yz_ytick_interval_mm = 1; % Set the interval between y ticks
                num_yz_yticks = actualSize(3)/new_yz_ytick_interval_mm; % number of y ticks
                new_yz_yticks = round(0:nzp_yz / num_yz_yticks:nzp_yz); new_yz_yticks(1) = 1;
                axial_range_mm = 0:new_yz_ytick_interval_mm:actualSize(3);
                new_yz_yticklabels = {};
                for i = 1:length(new_yz_yticks)
                    new_yz_yticklabels{i} = num2str(axial_range_mm(i));
                end
            end

            % Set labels
            % figure(figHandle)
            xlabel('y [mm]')
            ylabel('z [mm]')
            xticks(new_yz_xticks)
            xticklabels(new_yz_xticklabels)
            yticks(new_yz_yticks)
            yticklabels(new_yz_yticklabels)
            set(gca, 'TickDir', 'out')
            set(gca, 'box', 'off')
        end

        cv = getframe(tf);
        rgb = frame2im(cv);      % convert the frame to rgb data
        if x == 1
            yz_tagstruct.ImageLength = size(rgb, 1);
            yz_tagstruct.ImageWidth = size(rgb, 2);
            setTag(yz_stack, yz_tagstruct) % set the tags

%             tf.Position(4) = ceil(tf.Position(3) * hwRatio_yz);

        else
            writeDirectory(yz_stack)
            setTag(yz_stack, yz_tagstruct)
        end

        write(yz_stack, rgb)
    end
    close(yz_stack)

    %% xy planes
    xy_stack = Tiff([savepath, 'xy_stack.tif'], 'w');
%     xy_tagstruct.ImageLength = nzp;
%     xy_tagstruct.ImageWidth = nxp;
    xy_tagstruct.Photometric = Tiff.Photometric.RGB;
    xy_tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    xy_tagstruct.BitsPerSample = 8;
    xy_tagstruct.SamplesPerPixel = 3;
    xy_tagstruct.Software = 'MATLAB';
    xy_tagstruct.Compression = Tiff.Compression.Deflate;
    setTag(xy_stack, xy_tagstruct) % set the tags

    tf = figure;
    set(tf, 'Position', round([0, 0, nxp(1) * numVolumes, nyp(1) * actualSize(1)/actualSize(2)] .* scale)); % Define figure size
%     set(tf, 'Position', [0, 0, nxp(1) * numVolumes, nyp(1)] .* scale); % Define figure size
    tiledlayout(1, numVolumes)
    colormap(cmap)

    % adjust size of the image
    if exist('actualSize', 'var')
        hwRatio_xy = actualSize(1) / actualSize(2); % height to width ratio
    else
        hwRatio_xy = 1;
    end

    for z = 1:size(volumeData{vi}, 3) - mws + 1
        for vi = 1:numVolumes
            nexttile(vi) % cycle between tiles
            planeTemp = squeeze(max(volumeData{vi}(:, :, z:z + mws - 1), [], 3));

            imagesc(squeeze(planeTemp)')
            axis tight
            ax = gca;
            % set(ax, 'Color', 'k') % set background to black
            ax.PlotBoxAspectRatio = [1, hwRatio_xy, 1];
            if showColorbar
                colorbar
            end
            clim(cr(vi, :))

            if z == 1 && vi == 1
                % Set the tick labels, modifying what the default is
                % Set the x ticks
                default_xy_xticks = xticks;
                if mod(length(default_xy_xticks), 2) == 0 % Modify if the # of ticks is even, so we can get a tick at 0
                    new_xy_xticks = [1, default_xy_xticks];
                else
                    new_xy_xticks = default_xy_xticks;
                end

                x_range_mm = linspace(-actualSize(2)/2, actualSize(2)/2, length(new_xy_xticks)); % x on 'x' axis

                new_xy_xticklabels = {};
                for i = 1:length(new_xy_xticks)
                    new_xy_xticklabels{i} = num2str(x_range_mm(i));
                end

                % Set the y ticks
                default_xy_yticks = yticks;
                if mod(length(default_xy_yticks), 2) == 0 % Modify if the # of ticks is even, so we can get a tick at 0
                    new_xy_yticks = [1, default_xy_yticks];
                else
                    new_xy_yticks = default_xy_yticks;
                end

                y_range_mm = linspace(-actualSize(1)/2, actualSize(1)/2, length(new_xy_yticks)); % x on 'x' axis

                new_xy_yticklabels = {};
                for i = 1:length(new_xy_yticks)
                    new_xy_yticklabels{i} = num2str(y_range_mm(i));
                end

                % nyp_xy = nyp(1); % # z pixels in the xz image
                % new_xy_ytick_interval_mm = 1; % Set the interval between y ticks
                % num_xy_yticks = actualSize(1)/new_xy_ytick_interval_mm; % number of y ticks
                % new_xy_yticks = round(0:nyp_xy / num_xy_yticks:nyp_xy); new_xy_yticks(1) = 1;
                % % axial_range_mm = 0:new_xy_ytick_interval_mm:actualSize(1);
                % axial_range_mm = linspace(-actualSize(1)/2, actualSize(1)/2, length(new_xy_xticks)); % x on 'x' axis
                % new_xy_yticklabels = {};
                % for i = 1:length(new_xy_yticks)
                %     new_xy_yticklabels{i} = num2str(axial_range_mm(i));
                % end
            end

            % Set labels
            % figure(figHandle)
            xlabel('x [mm]')
            ylabel('y [mm]')
            xticks(new_xy_xticks)
            xticklabels(new_xy_xticklabels)
            yticks(new_xy_yticks)
            yticklabels(new_xy_yticklabels)
            set(gca, 'TickDir', 'out')
            set(gca, 'box', 'off')

        end
        cv = getframe(tf);
        rgb = frame2im(cv);      % convert the frame to rgb data
        if z == 1
            xy_tagstruct.ImageLength = size(rgb, 1);
            xy_tagstruct.ImageWidth = size(rgb, 2);
            setTag(xy_stack, xy_tagstruct) % set the tags

%             tf.Position(4) = ceil(tf.Position(3) * hwRatio_xy);

        else
            writeDirectory(xy_stack)
            setTag(xy_stack, xy_tagstruct)
        end

        write(xy_stack, rgb)
    end
    close(xy_stack)

    
end