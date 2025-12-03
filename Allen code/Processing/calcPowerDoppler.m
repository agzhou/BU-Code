% Input: IQf_separated (3x1 cell array with the negative, positive, and all
% frequency components of the filtered IQ (IQf)
% Optional input #2: noise from SVD for noise equalization

% Output: PDI (same except Power Doppler index)
function [PDI] = calcPowerDoppler(IQf_separated, varargin)
    useNoiseEq = false;
    if nargin > 1 % Use the noise equalization 
        useNoiseEq = true;
        noise = varargin{1};
    end

    if iscell(IQf_separated)
        PDI = cell(size(IQf_separated));
    
        frameDim = length(size(IQf_separated{1})); % Usually the frame dimension is the last dimension. 3 for 2D data and 4 for 3D data.
    
        for i = 1:length(PDI)
            if useNoiseEq
                PDI{i} = getPDI(IQf_separated{i}, frameDim) ./ noise;
            else
                PDI{i} = getPDI(IQf_separated{i}, frameDim);
            end
        end
    else
        frameDim = length(size(IQf_separated));

        if useNoiseEq
            PDI = getPDI(IQf_separated, frameDim) ./ noise;
        else
            PDI = getPDI(IQf_separated, frameDim);
        end
    end

end

% Helper function
function PDI = getPDI(IQf, frameDim)
    PDI = mean(abs(IQf) .^ 2, frameDim) ./ size(IQf, frameDim);
end