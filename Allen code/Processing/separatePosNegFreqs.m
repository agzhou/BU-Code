% Separate the negative and positive frequency components of the IQf data
% (filtered to have the blood signal)

% Output: 
%   IQf_separated (negative, positive, all)
function [IQf_separated, IQf_FT_separated, varargout] = separatePosNegFreqs(IQf)
    % [negativeComponent, positiveComponent, allComponents]
    frameDim = length(size(IQf)); % Usually the frame dimension is the last dimension. 3 for 2D data and 4 for 3D data.
    nf = size(IQf, frameDim); % # of frames in the IQf data

    % Pad with trailing zeros to improve computation speed
    % np = 2^nextpow2(2*size(IQf, frameDim));           % # of Fourier Transform points
    % Change 6/29/26:
    np = 2^nextpow2(size(IQf, frameDim));           % # of Fourier Transform points

    IQf_FT = fft(IQf, np, frameDim);
    IQf_FT_shifted = fftshift(IQf_FT, frameDim);

    posMask = ones(np, 1); posMask(1:np/2) = 0;
    negMask = ones(np, 1); negMask(np/2 + 1:np) = 0;

    if frameDim == 4 % 3D data
        % negativeFTComponent = IQf_FT_shifted(:, :, :, 1:np/2);
        % positiveFTComponent = IQf_FT_shifted(:, :, :, np/2 + 1:np);
        negativeFTComponent = IQf_FT_shifted; negativeFTComponent(:, :, :, ~negMask) = 0;
        positiveFTComponent = IQf_FT_shifted; positiveFTComponent(:, :, :, ~posMask) = 0;
        allFTComponents = IQf_FT_shifted;
    elseif frameDim == 3 % 2D data
        % negativeFTComponent = IQf_FT_shifted(:, :, 1:np/2);
        % positiveFTComponent = IQf_FT_shifted(:, :, np/2 + 1:np);
        negativeFTComponent = IQf_FT_shifted; negativeFTComponent(:, :, ~negMask) = 0;
        positiveFTComponent = IQf_FT_shifted; positiveFTComponent(:, :, ~posMask) = 0;
        allFTComponents = IQf_FT_shifted;
    end
    IQf_FT_separated = [{negativeFTComponent}; {positiveFTComponent}; {allFTComponents}];
    
    negativeComponent = ifft(negativeFTComponent, np, frameDim);
    positiveComponent = ifft(positiveFTComponent, np, frameDim);
    allComponents = ifft(allFTComponents, np, frameDim);

    if frameDim == 4 % 3D data
        negativeComponent = negativeComponent(:, :, :, 1:nf);
        positiveComponent = positiveComponent(:, :, :, 1:nf);
        allComponents = allComponents(:, :, :, 1:nf);
    elseif frameDim == 3 % 2D data
        negativeComponent = negativeComponent(:, :, 1:nf);
        positiveComponent = positiveComponent(:, :, 1:nf);
        allComponents = allComponents(:, :, 1:nf);
    end

    IQf_separated = [{negativeComponent}; {positiveComponent}; {allComponents}];

    % Optional outputs
    if nargout > 2
        varargout{1} = np; % # of Fourier transform points
    end
end