% Input: IQf_separated (3x1 cell array with the negative, positive, and all
% frequency components of the filtered IQ (IQf)

% Output: PDI (same except Power Doppler index)
function [PDI] = calcPowerDoppler(IQf_separated)
    if iscell(IQf_separated)
        PDI = cell(size(IQf_separated));
    
        frameDim = length(size(IQf_separated{1})); % Usually the frame dimension is the last dimension. 3 for 2D data and 4 for 3D data.
    
        for i = 1:length(PDI)
            PDI{i} = mean(abs(IQf_separated{i}) .^ 2, frameDim);
        end
    else
        frameDim = length(size(IQf_separated));
        PDI = mean(abs(IQf_separated) .^ 2, frameDim);
    end

end