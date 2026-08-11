% Unstack the spatial dimensions for 2D or 3D data
% ex: 2D data could be [nPixels, nt] and it outputs [nz, nx, nt] data

% Required inputs within the Processing Params struct:
%   PP.dimensionality: 2 for 2D data, 3 for 3D

function data_unstacked = unstackData(data_stacked, PP)
    warning('Verify if this code works, because I have never run it')
    ds = size(data_stacked);
    switch PP.dimensionality
        case 2
            data_unstacked = reshape(data_stacked, [PP.zDim, PP.xDim, ds(2:end)]);
        case 3
            data_unstacked = reshape(data_stacked, [PP.xDim, PP.yDim, PP.zDim, ds(2:end)]);
    end