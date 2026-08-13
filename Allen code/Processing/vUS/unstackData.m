% Unstack the spatial dimensions for 2D or 3D data
% ex: 2D data could be [nPixels, nt] and it outputs [nz, nx, nt] data
% The code works for any size in the 2nd dimension, e.g., # of Fourier
%     Transform points, # of time lags, etc.

% Required inputs within the Processing Params struct:
%   PP.dimensionality: 2 for 2D data, 3 for 3D

function data_unstacked = unstackData(data_stacked, PP)
    if iscell(data_stacked) % If input data is a cell array, apply the stacking to all cells individually and output cells
        data_unstacked = cell(size(data_stacked)); % Initialize cell output
        for j = 1:numel(data_stacked)
            ds = size(data_stacked{j}); % Data size [nVox, nFrames or nSamples etc.]
            switch PP.dimensionality
                case 2
                    data_unstacked{j} = reshape(data_stacked{j}, [PP.zp, PP.xp, ds(2)]);
                case 3
                    data_unstacked{j} = reshape(data_stacked{j}, [PP.xp, PP.yp, PP.zp, ds(2)]);
            end
        end
    else % If input data is not a cell (normal matrix)
        ds = size(data_stacked); % Data size [nVox, nFrames or nSamples etc.]
        switch PP.dimensionality
            case 2
                data_unstacked = reshape(data_stacked, [PP.zp, PP.xp, ds(2)]);
            case 3
                data_unstacked = reshape(data_stacked, [PP.xp, PP.yp, PP.zp, ds(2)]);
        end
    end
    % warning('Verify if this code works, because I have never run it')
    % ds = size(data_stacked); % Data size
    % switch PP.dimensionality
    %     case 2
    %         data_unstacked = reshape(data_stacked, [ds(PP.zDim), ds(PP.xDim), ds(PP.fDim)]);
    %     case 3
    %         data_unstacked = reshape(data_stacked, [ds(PP.xDim), ds(PP.yDim), ds(PP.zDim), ds(PP.fDim)]);
    % end
end