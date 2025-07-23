%% Make a tiff stack
% Input data is assumed to be (y, x, z, frames)
% Optional inputs: 1. actual image/physical size [y, x, z] in the same physical units
%                  2. colormap (string)
%                  3. MIP window range

function generateTiffStack_acrossframes(volumeData, varargin)

    % volumeData = volumeData ./ max(volumeData, [], 'all'); % Normalize intensities to be between 0 - 1
%     cr = [0, max(volumeData, [], 'all')]; % color range
    cr = [min(volumeData, [], 'all'), max(volumeData, [], 'all')]; % change 07/21/25
    showColorbar = true;
%     showColorbar = false;

    if nargin > 1
        actualSize = varargin{1};
        if nargin > 2
            cmap = varargin{2};
            if nargin > 3
                mwr = varargin{3}; % MIP window range
            end
        end
    end
    savepath = uigetdir('D:\Allen\Data\', 'Select the save path');
    savepath = [savepath, '\'];

    if length(size(volumeData)) == 4 % 4D data (3D space + time)
        nyp = size(volumeData, 1);
        nxp = size(volumeData, 2);
        nzp = size(volumeData, 3);
    
    
        %% Generate xz planes
        xz_stack = Tiff([savepath, 'xz_stack.tif'], 'w');
    %     xz_tagstruct.ImageLength = nzp;
    %     xz_tagstruct.ImageWidth = nxp;
        xz_tagstruct.Photometric = Tiff.Photometric.RGB;
        xz_tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        xz_tagstruct.BitsPerSample = 8;
        xz_tagstruct.SamplesPerPixel = 3;
        xz_tagstruct.Software = 'MATLAB';
        xz_tagstruct.Compression = Tiff.Compression.Deflate;
        setTag(xz_stack, xz_tagstruct) % set the tags
    
        tf = figure;
        colormap(cmap)
    
        % adjust size of the image
    %     if exist('actualSize', 'var')
    %         hwRatio_xz = actualSize(3) / actualSize(2); % height to width ratio
    %     else
    %         hwRatio_xz = 1;
    %     end
    
        for f = 1:size(volumeData, 4)
            planeTemp = squeeze(max(volumeData(mwr, :, :, f), [], 1));
    
            imagesc(squeeze(planeTemp)')
            title("Frame " + num2str(f))
            if showColorbar
                colorbar
            end
            % clim(cr)
            cv = getframe(tf);
            rgb = frame2im(cv);      % convert the frame to rgb data
            if f == 1
                xz_tagstruct.ImageLength = size(rgb, 1);
                xz_tagstruct.ImageWidth = size(rgb, 2);
                setTag(xz_stack, xz_tagstruct) % set the tags
                
    %             tf.Position(4) = ceil(tf.Position(3) * hwRatio_xz);
                
            else
                writeDirectory(xz_stack)
                setTag(xz_stack, xz_tagstruct)
            end
    
            write(xz_stack, rgb)
        end
        close(xz_stack)
    
        %% yz planes
    %     yz_stack = Tiff([savepath, 'yz_stack.tif'], 'w');
    % %     yz_tagstruct.ImageLength = nzp;
    % %     yz_tagstruct.ImageWidth = nxp;
    %     yz_tagstruct.Photometric = Tiff.Photometric.RGB;
    %     yz_tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    %     yz_tagstruct.BitsPerSample = 8;
    %     yz_tagstruct.SamplesPerPixel = 3;
    %     yz_tagstruct.Software = 'MATLAB';
    %     yz_tagstruct.Compression = Tiff.Compression.Deflate;
    %     setTag(yz_stack, yz_tagstruct) % set the tags
    % 
    %     tf = figure;
    %     colormap(cmap)
    % 
    %     % adjust size of the image
    % %     if exist('actualSize', 'var')
    % %         hwRatio_yz = actualSize(3) / actualSize(1); % height to width ratio
    % %     else
    % %         hwRatio_yz = 1;
    % %     end
    % 
    %     for x = 1:size(volumeData, 2) - mws + 1
    %         planeTemp = squeeze(max(volumeData(:, x:x + mws - 1, :), [], 2));
    % 
    %         imagesc(squeeze(planeTemp)')
    %         if showColorbar
    %             colorbar
    %         end
    %         clim(cr)
    %         cv = getframe(tf);
    %         rgb = frame2im(cv);      % convert the frame to rgb data
    %         if x == 1
    %             yz_tagstruct.ImageLength = size(rgb, 1);
    %             yz_tagstruct.ImageWidth = size(rgb, 2);
    %             setTag(yz_stack, yz_tagstruct) % set the tags
    % 
    % %             tf.Position(4) = ceil(tf.Position(3) * hwRatio_yz);
    % 
    %         else
    %             writeDirectory(yz_stack)
    %             setTag(yz_stack, yz_tagstruct)
    %         end
    % 
    %         write(yz_stack, rgb)
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
    %     xy_tagstruct.Compression = Tiff.Compression.Deflate;
    %     setTag(xy_stack, xy_tagstruct) % set the tags
    % 
    %     tf = figure;
    %     colormap(cmap)
    % 
    %     % adjust size of the image
    % %     if exist('actualSize', 'var')
    % %         hwRatio_xy = actualSize(3) / actualSize(2); % height to width ratio
    % %     else
    % %         hwRatio_xy = 1;
    % %     end
    % 
    %     for z = 1:size(volumeData, 3) - mws + 1
    %         planeTemp = squeeze(max(volumeData(:, :, z:z + mws - 1), [], 3));
    % 
    %         imagesc(squeeze(planeTemp)')
    %         if showColorbar
    %             colorbar
    %         end
    %         clim(cr)
    %         cv = getframe(tf);
    %         rgb = frame2im(cv);      % convert the frame to rgb data
    %         if z == 1
    %             xy_tagstruct.ImageLength = size(rgb, 1);
    %             xy_tagstruct.ImageWidth = size(rgb, 2);
    %             setTag(xy_stack, xy_tagstruct) % set the tags
    % 
    % %             tf.Position(4) = ceil(tf.Position(3) * hwRatio_xy);
    % 
    %         else
    %             writeDirectory(xy_stack)
    %             setTag(xy_stack, xy_tagstruct)
    %         end
    % 
    %         write(xy_stack, rgb)
    %     end
    %     close(xy_stack)
    else % If the input data is 3D (2D space + time)
        nzp = size(volumeData, 1);
        nxp = size(volumeData, 2); 
    
        % Generate planes
        xz_stack = Tiff([savepath, 'xz_stack.tif'], 'w');
    %     xz_tagstruct.ImageLength = nzp;
    %     xz_tagstruct.ImageWidth = nxp;
        xz_tagstruct.Photometric = Tiff.Photometric.RGB;
        xz_tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        xz_tagstruct.BitsPerSample = 8;
        xz_tagstruct.SamplesPerPixel = 3;
        xz_tagstruct.Software = 'MATLAB';
        xz_tagstruct.Compression = Tiff.Compression.Deflate;
        setTag(xz_stack, xz_tagstruct) % set the tags
    
        tf = figure;
        colormap(cmap)
    
        % adjust size of the image
    %     if exist('actualSize', 'var')
    %         hwRatio_xz = actualSize(3) / actualSize(2); % height to width ratio
    %     else
    %         hwRatio_xz = 1;
    %     end
    
        for f = 1:size(volumeData, 3)
            planeTemp = squeeze(volumeData(:, :, f))';
    
            imagesc(squeeze(planeTemp)')
            title("Frame " + num2str(f))
            if showColorbar
                colorbar
            end
            clim(cr)
            cv = getframe(tf);
            rgb = frame2im(cv);      % convert the frame to rgb data
            if f == 1
                xz_tagstruct.ImageLength = size(rgb, 1);
                xz_tagstruct.ImageWidth = size(rgb, 2);
                setTag(xz_stack, xz_tagstruct) % set the tags
                
    %             tf.Position(4) = ceil(tf.Position(3) * hwRatio_xz);
                
            else
                writeDirectory(xz_stack)
                setTag(xz_stack, xz_tagstruct)
            end
    
            write(xz_stack, rgb)
        end
        close(xz_stack)
    
    end

end