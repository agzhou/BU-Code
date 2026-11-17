% Separate the negative and positive frequency components of the IQf data
% (filtered to have the blood signal)

% Input:
%   IQf: filtered IQ, frames in the last dimension (3 = 2D data, 4 = 3D data)
%   Optional name-value pairs (defaults reproduce the original hard 0/1 split):
%     'DCTaper'   : width of a raised-cosine ramp at the DC edge, as a fraction
%                   of the Nyquist band (0 to 1). The ramp lies entirely on each
%                   component's own side of 0 Hz, so there is no leakage between
%                   the negative and positive components. Default 0 (hard edge).
%     'NyqTaper'  : width of a raised-cosine ramp down to 0 at +/-Nyquist, as a
%                   fraction of the Nyquist band (0 to 1). Softens the edge that
%                   otherwise gives the white noise floor a half-band
%                   autocorrelation (zero at even lags, ~2/(pi*n) at odd lags),
%                   i.e. a zig-zag in g1. Default 0 (hard edge).
%     'TaperAll'  : if true, also apply the Nyquist taper to the "all
%                   frequencies" component (default false: left untouched).
%
% Output:
%   IQf_separated: 3x1 cell array (negative, positive, all) of directionally-separated IQf values
%   IQf_FT_separated: 3x1 cell array (negative, positive, all) of directionally-separated, Fourier transformed IQf
%                     (NOTE: THIS IS ZERO-SHIFTED!!!)
%   varargout{1}: # of Fourier transform points

function [IQf_separated, IQf_FT_separated, varargout] = separatePosNegFreqs_taper(IQf, varargin)
    ip = inputParser;
    ip.addParameter('DCTaper', 0, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
    ip.addParameter('NyqTaper', 0, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
    ip.addParameter('TaperAll', false, @(x) islogical(x) || isnumeric(x));
    ip.parse(varargin{:});
    dcTaper = ip.Results.DCTaper;
    nyqTaper = ip.Results.NyqTaper;
    taperAll = logical(ip.Results.TaperAll);

    % [negativeComponent, positiveComponent, allComponents]
    frameDim = length(size(IQf)); % Usually the frame dimension is the last dimension. 3 for 2D data and 4 for 3D data.
    nf = size(IQf, frameDim); % # of frames in the IQf data

    % Pad with trailing zeros to improve computation speed
    % np = 2^nextpow2(2*size(IQf, frameDim));           % # of Fourier Transform points
    % Change 6/29/26:
    np = 2^nextpow2(size(IQf, frameDim));           % # of Fourier Transform points

    IQf_FT = fft(IQf, np, frameDim);
    IQf_FT_shifted = fftshift(IQf_FT, frameDim);

    % Frequency weights on the fftshifted axis. Bin index m = -np/2 : np/2 - 1
    % and u = |f|/(fs/2) in [0, 1]. With both tapers at 0 these are exactly the
    % original hard masks: negative = bins m < 0 (includes the -Nyquist bin),
    % positive = bins m >= 0 (includes DC).
    m = (-np/2:np/2 - 1)';
    u = abs(m) ./ (np/2);
    rc = @(x) 0.5 * (1 - cos(pi * min(max(x, 0), 1))); % raised cosine ramp 0 -> 1 for x in [0, 1]
    if dcTaper > 0
        rise = rc(u ./ dcTaper);
    else
        rise = ones(np, 1);
    end
    if nyqTaper > 0
        fall = rc((1 - u) ./ nyqTaper);
    else
        fall = ones(np, 1);
    end
    negWeight = (m < 0)  .* rise .* fall;
    posWeight = (m >= 0) .* rise .* fall;
    if taperAll
        allWeight = fall;
    else
        allWeight = ones(np, 1);
    end

    % Put the weights along the frame (frequency) dimension
    wShape = ones(1, frameDim); wShape(frameDim) = np;
    negWeight = reshape(negWeight, wShape);
    posWeight = reshape(posWeight, wShape);
    allWeight = reshape(allWeight, wShape);

    negativeFTComponent = IQf_FT_shifted .* negWeight;
    positiveFTComponent = IQf_FT_shifted .* posWeight;
    allFTComponents = IQf_FT_shifted .* allWeight;
    IQf_FT_separated = [{negativeFTComponent}; {positiveFTComponent}; {allFTComponents}];

    % Undo the fftshift and perform inverse FFTs to recover the signals
    negativeComponent = ifft(ifftshift(negativeFTComponent, frameDim), np, frameDim);
    positiveComponent = ifft(ifftshift(positiveFTComponent, frameDim), np, frameDim);
    allComponents = ifft(ifftshift(allFTComponents, frameDim), np, frameDim);

    % Crop to the # of frames, since there was padding above to get the #
    % of FT points to a power of 2 for computational efficiency
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
