% Description: plot a FC correlation matrix

% Inputs:
%   - CM: correlation matrix
%   - roi: struct with ROI info

function plotCM(CM, roi, varargin)
    separate_hemis = false;
    if nargin > 2
        separate_hemis = varargin{1};
    end
    if separate_hemis   % Separate ROIs by hemisphere
        figure; imagesc(CM); colormap jet; axis square; colorbar; clim([-1, 1])
        xticks(1:roi.num_regions * 2)
        xticklabels(roi.acronyms_hemis_interleaved) % Set the tick labels to be the ROI acronyms
        yticks(1:roi.num_regions * 2)
        yticklabels(roi.acronyms_hemis_interleaved) % Set the tick labels to be the ROI acronyms
    else                % No hemisphere separation of ROIs
        figure; imagesc(CM); colormap jet; axis square; colorbar; clim([-1, 1])
        xticks(1:roi.num_regions)
        xticklabels(roi.acronyms)
        yticks(1:roi.num_regions)
        yticklabels(roi.acronyms) % Set the tick labels to be the ROI acronyms
    end
end