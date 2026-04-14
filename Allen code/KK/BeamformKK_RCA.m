% Description: KK beamforming for row column array

% Inputs:
%   - RawDataKK is the KK-transformed RF data, with dimensions [# samples, # TX angles, # RX angles]
%   - anglesRX is a list of the RX angles. Column 1 = theta_x, column 2 = theta_y [rad]
%   - anglesTX is a list of the TX angles. Column 1 = theta_x, column 2 = theta_y [rad]
%   - BFgrid: a struct with X, Y, and Z fields corresponding to the beamforming grid
%   - param: a struct with at least the parameters 'fs' (sampling frequency
%            in Hz), 'c' (speed of sound [m/s]), 't0' (initial time offset [s])
%   - (Optional input) whatToReturn: 'compounded' volume or 'allAngles'
%                                    (return all individual volumes for each pair of TX-RX plane waves)

% test = BeamformKK_MatrixArray(RawDataKK, anglesRX, BFgrid, param);

function [BFData, varargout] = BeamformKK_RCA(RawDataKK, anglesRX, anglesTX, BFgrid, param, interp_method, varargin)
    
    % Determine if all angles will be returned, or just the compounded result
    whatToReturn = 'compounded'; % Default is to return only the compounded result
    if nargin > 5
        whatToReturn = varargin{1};
    end
    
    ns = size(RawDataKK, 1); % # of samples
    naTX = size(RawDataKK, 2); % # of TX angles
    naRX = size(RawDataKK, 3); % # of RX angles
    [nx, ny, nz] = size(BFgrid.X); % Get the size of the grid

    if naRX ~= size(anglesRX, 1)
        error('Number of RX angles is inconsistent between RawDataKK and anglesRX')
    end

    % **** BREAK UP ANGLES INTO RC AND CR **** %
    indsRX_CR = 1:naRX/2;
    indsRX_RC = naRX/2 + 1:naRX;
    anglesRX_CR = anglesRX(indsRX_CR, 2);
    anglesRX_RC = anglesRX(indsRX_RC, 1);

    indsTX_CR = 1:naTX/2;
    indsTX_RC = naTX/2 + 1:naTX;
    anglesTX_CR = anglesTX(indsTX_CR, 1);
    anglesTX_RC = anglesTX(indsTX_RC, 2);

    anglesTX_vec = [anglesTX_CR(:); anglesTX_RC(:)]; % Make anglesTX a vector of only the relevant angles (since one of theta_x or theta_y is zero with the RCA OPW)
    anglesRX_vec = [anglesRX_CR(:); anglesRX_RC(:)]; % Make anglesTX a vector of only the relevant angles (since one of theta_x or theta_y is zero with the RCA OPW)

    % Make look-up tables for each receive angle and transmit angle
    % LUTRX = cell(naRX, 1);
    % LUTTX = cell(naTX, 1);
    LUTRX = zeros(nx, ny, nz, naRX);
    LUTTX = zeros(nx, ny, nz, naTX);

    for tai = 1:naTX     % Transmit angle index
    % for tai = 1
        % LUTTX{tai} = genLUT(anglesTX(tai, :), BFgrid, param.c);
        % LUTTX(:, :, :, tai) = genLUT(anglesTX(tai, :), BFgrid, param.c);
        LUTTX(:, :, :, tai) = genLUT(anglesTX(tai, :), BFgrid, param.c);
    end
    disp('TX LUTs generated')

    for rai = 1:naRX % Receive angle index
    % for rai = 1
        % LUTRX{rai} = genLUT(anglesRX(rai, :), BFgrid, param.c);
        % LUTRX(:, :, :, rai) = genLUT(anglesRX(rai, :), BFgrid, param.c);
        LUTRX(:, :, :, rai) = genLUT(anglesRX(rai, :), BFgrid, param.c);
    end
    disp('RX LUTs generated')

    % Go through each transmit angle and beamform with its constituent
    % receive angles

    switch whatToReturn
        %% Compounded output
        case 'compounded'
            BFData = zeros(nx, ny, nz); % Initialize final beamformed volume
        
            temp = zeros(nx, ny, nz); % Initialize a volume to keep adding to
            for tai = 1:naTX/2     % Transmit angle index
                disp(tai)
                % tempLUTTX = LUTTX{tai}; % Temporarily store the TX time delays for angle index tai
                tempLUTTX = squeeze(LUTTX(:, :, :, tai)); % Temporarily store the TX time delays for angle index tai
                % tempLUTTX = genLUT(anglesTX(tai, :), BFgrid, param.c, param.t0);
                for rai = 1:naRX/2
                    tempLUTRX = squeeze(LUTRX(:, :, :, rai)); % Temporarily store the RX time delays for angle index rai
                    % tempLUTRX = genLUT(anglesRX(rai, :), BFgrid, param.c, param.t0);
                    [verytemp] = applyTimeDelays(nx, ny, nz, ns, tempLUTTX, tempLUTRX, param, squeeze(RawDataKK(:, tai, rai)), interp_method);
                    temp = temp + verytemp;
                end
                
            end
            BFData = BFData + temp;
            
            temp = zeros(nx, ny, nz); % Initialize a volume to keep adding to
            for tai = naTX/2 + 1:naTX     % Transmit angle index
                disp(tai)
                % temp = zeros(nx, ny, nz); % Initialize a volume to keep adding to
                % tempLUTTX = LUTTX{tai}; % Temporarily store the TX time delays for angle index tai
                tempLUTTX = squeeze(LUTTX(:, :, :, tai)); % Temporarily store the TX time delays for angle index tai
                % tempLUTTX = genLUT(anglesTX(tai, :), BFgrid, param.c, param.t0);
                for rai = naRX/2 + 1:naRX
                    tempLUTRX = squeeze(LUTRX(:, :, :, rai)); % Temporarily store the RX time delays for angle index rai
                    % tempLUTRX = genLUT(anglesRX(rai, :), BFgrid, param.c, param.t0);
                    [verytemp] = applyTimeDelays(nx, ny, nz, ns, tempLUTTX, tempLUTRX, param, squeeze(RawDataKK(:, tai, rai)), interp_method);
                    temp = temp + verytemp;
                end
                
            end
            BFData = BFData + temp;
        %% Non-compounded output
        case 'allAngles'
            BFData = zeros(nx, ny, nz, naTX, naRX); % Initialize final beamformed volume
            for tai = 1:naTX/2     % Transmit angle index
                disp(tai)
                
                % tempLUTTX = LUTTX{tai}; % Temporarily store the TX time delays for angle index tai
                tempLUTTX = squeeze(LUTTX(:, :, :, tai)); % Temporarily store the TX time delays for angle index tai
                % tempLUTTX = genLUT(anglesTX(tai, :), BFgrid, param.c, param.t0);
                for rai = 1:naRX/2
                    tempLUTRX = squeeze(LUTRX(:, :, :, rai)); % Temporarily store the RX time delays for angle index rai
                    % tempLUTRX = genLUT(anglesRX(rai, :), BFgrid, param.c, param.t0);
                    
                    [verytemp] = applyTimeDelays(nx, ny, nz, ns, tempLUTTX, tempLUTRX, param, squeeze(RawDataKK(:, tai, rai)), interp_method);
                    % temp = temp + verytemp;
                    BFData(:, :, :, tai, rai) = verytemp; % Save each TX and RX angle's BF data separately
                end
                
            end
            
            % temp = zeros(nx, ny, nz); % Initialize a volume to keep adding to
            for tai = naTX/2 + 1:naTX     % Transmit angle index
                disp(tai)
                % temp = zeros(nx, ny, nz); % Initialize a volume to keep adding to
                % tempLUTTX = LUTTX{tai}; % Temporarily store the TX time delays for angle index tai
                tempLUTTX = squeeze(LUTTX(:, :, :, tai)); % Temporarily store the TX time delays for angle index tai
                % tempLUTTX = genLUT(anglesTX(tai, :), BFgrid, param.c, param.t0);
                for rai = naRX/2 + 1:naRX
                    tempLUTRX = squeeze(LUTRX(:, :, :, rai)); % Temporarily store the RX time delays for angle index rai
                    % tempLUTRX = genLUT(anglesRX(rai, :), BFgrid, param.c, param.t0);
                    
                    [verytemp] = applyTimeDelays(nx, ny, nz, ns, tempLUTTX, tempLUTRX, param, squeeze(RawDataKK(:, tai, rai)), interp_method);
                    % temp = temp + verytemp;
                    BFData(:, :, :, tai, rai) = verytemp; % Save each TX and RX angle's BF data separately
                end
                
            end
    end

    % Return the TX/RX LUTs as optional outputs
    if nargout > 1
        varargout{1} = LUTTX;
        if nargout > 2
            varargout{2} = LUTRX;
        end
    end
    
    % % Return the TX angles, RX angles, and TX/RX LUTs as optional outputs
    % if nargout > 1
    %     varargout{1} = anglesTX_vec;
    %     if nargout > 2
    %         varargout{2} = anglesRX_vec;
    %         if nargout > 3
    %             varargout{3} = LUTTX;
    %             if nargout > 4
    %                 varargout{4} = LUTRX;
    %             end
    %         end
    %     end
    % end

end

%% Helper functions

% Description: Generate LUTs for a given angle and grid, for the plane wave
%              time delay. Assumes the probe is aligned with z = 0.
% Inputs:
%   - theta: angle [rad]
%   - BFgrid: a struct with X, Y, and Z fields corresponding to the beamforming grid
%   - c: speed of sound [m/s]
% Outputs:
%   - LUT: a matrix of time delays. Dimensions are the same as the grid.
% function [LUT] = genLUT(theta, BFgrid, c)
function [LUT] = genLUT(theta, BFgrid, c)

    [nx, ny, nz] = size(BFgrid.X); % Get the size of the grid
    LUT = zeros(nx, ny, nz); % Initialize the LUT
    u = [sin(theta(2)), -sin(theta(1))*cos(theta(2)), cos(theta(1))*cos(theta(2))]; % y rotation and then x rotation

    % Get the distance version of the time delays
    for xi = 1:nx
        for yi = 1:ny
            for zi = 1:nz
                        
                LUT(xi, yi, zi) = abs(dot([BFgrid.X(xi, yi, zi), BFgrid.Y(xi, yi, zi), BFgrid.Z(xi, yi, zi)], u));
            end
        end
    end
    LUT = LUT ./ c; % Convert from distances to time delays
    % LUT = LUT ./ c + t0; % Convert from distances to time delays


end

% Get the volume given some beaforming grid, KK-transformed RF Data, and
% time delays for each voxel in the volume
%   - interp_method: 'round' or 'linear' as how to get the time-delayed
%                    value for each voxel
function [vol] = applyTimeDelays(nx, ny, nz, ns, tempLUTTX, tempLUTRX, param, RawDataKK_vec, interp_method)    
    vol = zeros(nx, ny, nz); % Initialize a volume to keep adding to
    for xi = 1:nx
        for yi = 1:ny
            for zi = 1:nz
                % sampleDelay = ( tempLUTTX(xi, yi, zi) + tempLUTRX(xi, yi, zi) ).*param.fs;
                % sampleDelay = round( (tempLUTTX(xi, yi, zi) + tempLUTRX(xi, yi, zi)).*param.fs + 1);
                % sampleDelay = round( (tempLUTTX(xi, yi, zi) + tempLUTRX(xi, yi, zi)).*param.fs) + 1;

                sampleDelay = (tempLUTTX(xi, yi, zi) + tempLUTRX(xi, yi, zi) + param.t0).*param.fs + 1;
                switch interp_method
                    case 'round'
                        sampleDelay = round(sampleDelay);
                        if (sampleDelay < ns - 1) && (sampleDelay >= 1) % if statement for out of bounds delays
                            % disp('flag')
                            vol(xi, yi, zi) = RawDataKK_vec(sampleDelay);
                        end
                    case 'linear'
                        if (sampleDelay < ns - 1) && (sampleDelay >= 1) % if statement for out of bounds delays
                            % disp('flag')
                            vol(xi, yi, zi) = interpLinear(RawDataKK_vec, sampleDelay);
                        end
                end

                % sampleDelay = (tempLUTTX(xi, yi, zi) + tempLUTRX(xi, yi, zi)).*param.fs + 1; % delay in [samples]
                % disp(sampleDelay)
                % verytemp(xi, yi, zi) = interpLinear(squeeze(RawDataKK(:, tai, rai)), sampleDelay);
            end
        end
    end
end

% Linear interpolation
% 'data' should be all the samples for RawDataKK for one TX angle and RX angle pair
% 'delay' is in units of [samples]
function [value] = interpLinear(data, delay)
    ns = length(data); % # samples in the data
    % lowerSample = mod(floor(delay), ns);
    % upperSample = mod(ceil(delay), ns);

    lowerSample = floor(delay);
    % disp(lowerSample)
    upperSample = ceil(delay);

    % I guess I added some "circshifting" for if the delay is out of range
    if lowerSample > ns
        % disp(lowerSample)
        lowerSample = lowerSample - ns;
        % disp(lowerSample)
    end
    if upperSample > ns
        upperSample = upperSample - ns;
    end

    lower = data(lowerSample);
    upper = data(upperSample);
    
    % lower = data(floor(delay));
    % upper = data(ceil(delay));
    value = lower + (delay - floor(delay))*(upper - lower);

end