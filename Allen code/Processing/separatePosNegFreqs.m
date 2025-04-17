% Separate the negative and positive frequency components of the IQf data
% (filtered to have the blood signal)

% Output: 
%   IQf_separated (negative, positive, all)
function [IQf_separated, IQf_FT_separated]  = separatePosNegFreqs(IQf)
    % [negativeComponent, positiveComponent, allComponents]
    frameDim = length(size(IQf)); % Usually the frame dimension is the last dimension. 3 for 2D data and 4 for 3D data.
    nf = size(IQf, frameDim); % # of frames in the IQf data

    np = 2^nextpow2(2*size(IQf, frameDim));           % # of Fourier Transform points

    IQf_FT = fft(IQf, np, frameDim);
    IQf_FT_shifted = fftshift(IQf_FT, frameDim);

    if frameDim == 4 % 3D data
        negativeFTComponent = IQf_FT_shifted(:, :, :, 1:np/2);
        positiveFTComponent = IQf_FT_shifted(:, :, :, np/2 + 1:np);
        allFTComponents = IQf_FT_shifted;
    elseif frameDim == 3 % 2D data
        negativeFTComponent = IQf_FT_shifted(:, :, 1:np/2);
        positiveFTComponent = IQf_FT_shifted(:, :, np/2 + 1:np);
        allFTComponents = IQf_FT_shifted;
    end
    IQf_FT_separated = [{negativeFTComponent}; {positiveFTComponent}; {allFTComponents}];
    
    negativeComponent = ifft(negativeFTComponent, np, frameDim);
    positiveComponent = ifft(positiveFTComponent, np, frameDim);
    allComponents = ifft(allFTComponents, np, frameDim);

    negativeComponent = negativeComponent(:, :, :, 1:nf);
    positiveComponent = positiveComponent(:, :, :, 1:nf);
    allComponents = allComponents(:, :, :, 1:nf);

    IQf_separated = [{negativeComponent}; {positiveComponent}; {allComponents}];

end