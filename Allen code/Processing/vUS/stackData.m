% Stack the spatial dimensions for 2D or 3D data
% ex: 2D data could be [nz, nx, nt] and it outputs [nPixels, nt] data

% Required inputs within the Processing Params struct:
%   PP.dimensionality: 2 for 2D data, 3 for 3D

function data_stacked = stackData(data_unstacked, PP)
    % warning('Verify if this code works, because I have never run it')

    if iscell(data_unstacked) % If input data is a cell array, apply the stacking to all cells individually and output cells
        data_stacked = cell(size(data_unstacked)); % Initialize cell output
        for j = 1:numel(data_unstacked)
            ds = size(data_unstacked{j}); % Data size
            switch PP.dimensionality
                case 2
                    data_stacked{j} = reshape(data_unstacked{j}, [ds(PP.zDim) * ds(PP.xDim), ds(PP.fDim)]);
                case 3
                    data_stacked{j} = reshape(data_unstacked{j}, [ds(PP.xDim) * ds(PP.yDim) * ds(PP.zDim), ds(PP.fDim)]);
            end
        end
    else % If input data is not a cell (normal matrix)
        ds = size(data_unstacked); % Data size
        if length(ds) < PP.fDim % Add info in case the expected time dimension doesn't exist
            ds(PP.fDim) = 1;
        end

        switch PP.dimensionality
            case 2
                data_stacked = reshape(data_unstacked, [ds(PP.zDim) * ds(PP.xDim), ds(PP.fDim)]);
            case 3
                data_stacked = reshape(data_unstacked, [ds(PP.xDim) * ds(PP.yDim) * ds(PP.zDim), ds(PP.fDim)]);
        end
    end
end