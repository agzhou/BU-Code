% g1T calculation in 1D, 2D, or 3D
% Input: 1. some data, likely IQ (coherently summed across angles).
%           The data should have the spatial dimensions first and then a
%           time/frame dimension last.
%        2. (Optional) number of tau values to calculate g1 at, starting
%           from 0
% Output: temporal g1

% Notes: 

function [g1] = g1T(data, varargin)
    tstart = clock;
    frameDim = length(size(data)); % Usually the frame dimension is the last dimension. 3 for 2D data and 4 for 3D data.
    nf = size(data, frameDim); % # of frames
    
    if nargin > 1 % If the # of points to calculate g1 at is specified
        np = varargin{1};
    else
        np = nf;
    end

    dataSize = size(data);
    dataSize(frameDim) = np; % Set the # of points to calculate g1 at

    g1 = zeros(dataSize); % g1 is calculated for each spatial dimension and tau step
                                % e.g., in 3D: (x, y, z, tau step)
    numer = zeros(dataSize);
    denom = mean((conj(data) .* data), frameDim); % temporal (frame) average
    switch frameDim
        case 2 % 1D
            for f = 1:np % go through each frame to get the g1 at each tau step
                numer = mean(conj(data(:, 1:(nf - f + 1))) .* data(:, f:end), frameDim);
                g1(:, f) = numer ./ denom;
            end
        case 3 % 2D
            for f = 1:np % go through each frame to get the g1 at each tau step
                numer = mean(conj(data(:, :, 1:(nf - f + 1))) .* data(:, :, f:end), frameDim);
                g1(:, :, f) = numer ./ denom;
            end
        case 4 % 3D
            for f = 1:np % go through each frame to get the g1 at each tau step
                numer = mean(conj(data(:, :, :, 1:(nf - f + 1))) .* data(:, :, :, f:end), frameDim);
                g1(:, :, :, f) = numer ./ denom;
            end
    end

    tend = clock;
%     disp(strcat("Temporal g_{1} processing done, elapsed time is ", num2str(etime(tend, tstart)), "s"))
end

% function [g1] = g1T(data)
%     tstart = clock;
%     frameDim = length(size(data)); % Usually the frame dimension is the last dimension. 3 for 2D data and 4 for 3D data.
%     nf = size(data, frameDim); % # of frames
%     dataSize = size(data);
% 
%     g1 = zeros(dataSize); % g1 is calculated for each spatial dimension and tau step
%                                 % e.g., in 3D: (x, y, z, tau step)
%     numer = zeros(dataSize);
%     denom = mean((conj(data) .* data), frameDim); % temporal (frame) average
%     switch frameDim
%         case 2 % 1D
%             for f = 1:nf % go through each subframe to get the g1 at each tau step
%                 numer = mean(conj(data(:, 1:(nf - f + 1))) .* data(:, f:end), frameDim);
%                 g1(:, f) = numer ./ denom;
%             end
%         case 3 % 2D
%             for f = 1:nf % go through each subframe to get the g1 at each tau step
%                 numer = mean(conj(data(:, :, 1:(nf - f + 1))) .* data(:, :, f:end), frameDim);
%                 g1(:, :, f) = numer ./ denom;
%             end
%         case 4 % 3D
%             for f = 1:nf % go through each subframe to get the g1 at each tau step
%                 numer = mean(conj(data(:, :, :, 1:(nf - f + 1))) .* data(:, :, :, f:end), frameDim);
%                 g1(:, :, :, f) = numer ./ denom;
%             end
%     end
% 
%     tend = clock;
%     disp(strcat("Temporal g_{1} processing done, elapsed time is ", num2str(etime(tend, tstart)), "s"))
% end