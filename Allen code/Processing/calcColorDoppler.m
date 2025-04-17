% Input: IQf_separated (3x1 cell array with the negative, positive, and all
% frequency components of the filtered IQ (IQf)

% Output: CDI (same except Color Doppler index)
function [CDI] = calcColorDoppler(IQf_FT_separated, P)
    CDI = cell(size(IQf_FT_separated));

    freqDim = length(size(IQf_FT_separated{1})); % Usually the frequency dimension is the last dimension. 3 for 2D data and 4 for 3D data.
    for i = 1:length(CDI)
%     for i = 3
        fi = linspace(-P.frameRate/2, P.frameRate/2, size(IQf_FT_separated{i}, freqDim)); % frequencies
        
        % Note: this line only works for 3D right now
        fi = repmat(permute(fi, [1, 4, 3, 2]), size(IQf_FT_separated{i}, 1), size(IQf_FT_separated{i}, 2), size(IQf_FT_separated{i}, 3), 1);
        CDI{i} = sum( abs(IQf_FT_separated{i} .^ 2) .* fi, freqDim) ./ sum( abs(IQf_FT_separated{i} .^ 2), freqDim);
    end

end