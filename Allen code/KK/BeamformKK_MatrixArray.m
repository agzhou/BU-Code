% Description: KK beamforming for matrix array

% Inputs:
%   - RawDataKK is the KK-transformed RF data, with dimensions [# samples, # TX angles, # RX angles]
%   - anglesRX is a list of the RX angles. Column 1 = theta_x, column 2 = theta_y [rad]
%   - anglesTX is a list of the TX angles. Column 1 = theta_x, column 2 = theta_y [rad]
%   - BFgrid: a struct with X, Y, and Z fields corresponding to the beamforming grid


% test = BeamformKK_MatrixArray(RawDataKK, anglesRX, BFgrid, param);

function [BFData, varargout] = BeamformKK_MatrixArray(RawDataKK, anglesRX, anglesTX, BFgrid, param)

    % **** ADD IQ DEMODULATION STEP **** %
    % more so for speedup

    ns = size(RawDataKK, 1); % # of samples
    naTX = size(RawDataKK, 2); % # of TX angles
    naRX = size(RawDataKK, 3); % # of RX angles
    [nx, ny, nz] = size(BFgrid.X); % Get the size of the grid

    if naRX ~= size(anglesRX, 1)
        error('Number of RX angles is inconsistent between RawDataKK and anglesRX')
    end

    % Make look-up tables for each receive angle and transmit angle
    % LUTRX = cell(naRX, 1);
    % LUTTX = cell(naTX, 1);
    LUTRX = zeros(nx, ny, nz, naRX);
    LUTTX = zeros(nx, ny, nz, naTX);

    for tai = 1:naTX     % Transmit angle index
    % for tai = 1
        % LUTTX{tai} = genLUT(anglesTX(tai, :), BFgrid, param.c);
        % LUTTX(:, :, :, tai) = genLUT(anglesTX(tai, :), BFgrid, param.c);
        LUTTX(:, :, :, tai) = genLUT(anglesTX(tai, :), BFgrid, param.c, param.t0);
    end
    disp('TX LUTs generated')
    for rai = 1:naRX % Receive angle index
    % for rai = 1
        % LUTRX{rai} = genLUT(anglesRX(rai, :), BFgrid, param.c);
        % LUTRX(:, :, :, rai) = genLUT(anglesRX(rai, :), BFgrid, param.c);
        LUTRX(:, :, :, rai) = genLUT(anglesRX(rai, :), BFgrid, param.c, param.t0);
    end
    disp('RX LUTs generated')

    % Go through each transmit angle and beamform with its constituent
    % receive angles
    BFData = zeros(nx, ny, nz); % Initialize final beamformed volume
    % BFData = zeros(nx, ny, nz, naTX, naRX); % Initialize final beamformed volume

    for tai = 1:naTX     % Transmit angle index
    % for tai = 1
        disp(tai)
        temp = zeros(nx, ny, nz); % Initialize a volume to keep adding to
        % tempLUTTX = LUTTX{tai}; % Temporarily store the TX time delays for angle index tai
        tempLUTTX = squeeze(LUTTX(:, :, :, tai)); % Temporarily store the TX time delays for angle index tai
        % tempLUTTX = genLUT(anglesTX(tai, :), BFgrid, param.c, param.t0);
        for rai = 1:naRX
            tempLUTRX = squeeze(LUTRX(:, :, :, rai)); % Temporarily store the RX time delays for angle index rai
            % tempLUTRX = genLUT(anglesRX(rai, :), BFgrid, param.c, param.t0);

            verytemp = zeros(nx, ny, nz); % Initialize a volume to keep adding to
            for xi = 1:nx
                for yi = 1:ny
                    for zi = 1:nz
                        % sampleDelay = ( tempLUTTX(xi, yi, zi) + tempLUTRX(xi, yi, zi) ).*param.fs;
                        % sampleDelay = round( (tempLUTTX(xi, yi, zi) + tempLUTRX(xi, yi, zi)).*param.fs + 1);
                        sampleDelay = round( (tempLUTTX(xi, yi, zi) + tempLUTRX(xi, yi, zi)).*param.fs) + 1;
                        % INTERPOLATE???????????????
                        
                        if (sampleDelay < ns - 1) && (sampleDelay >= 1) % if statement for out of bounds delays
                            % disp('flag')
                            verytemp(xi, yi, zi) = RawDataKK(sampleDelay, tai, rai);
                        end

                        % sampleDelay = (tempLUTTX(xi, yi, zi) + tempLUTRX(xi, yi, zi)).*param.fs + 1; % delay in [samples]
                        % disp(sampleDelay)
                        % verytemp(xi, yi, zi) = interpLinear(squeeze(RawDataKK(:, tai, rai)), sampleDelay);
                    end
                end
            end
            temp = temp + verytemp;
            % BFData(:, :, :, tai, rai) = verytemp; % Save each TX and RX angle's BF data separately
        end
        BFData = BFData + temp;
    end

    % Return the LUTs as optional outputs
    if nargout > 1
        varargout{1} = LUTTX;
        if nargout > 2
            varargout{2} = LUTRX;
        end
    end
end

%% Helper functions

% Description: Generate LUTs for a given angle and grid, for the plane wave time delay
% Inputs:
%   - theta: angle [rad]
%   - BFgrid: a struct with X, Y, and Z fields corresponding to the beamforming grid
%   - c: speed of sound [m/s]
% Outputs:
%   - LUT: a matrix of time delays. Dimensions are the same as the grid.
% function [LUT] = genLUT(theta, BFgrid, c)
function [LUT] = genLUT(theta, BFgrid, c, t0)

    [nx, ny, nz] = size(BFgrid.X); % Get the size of the grid
    LUT = zeros(nx, ny, nz); % Initialize the LUT
    % u = [sin(theta(2)), -sin(theta(1)), cos(theta(1))*cos(theta(2))]; % Unit direction vector for the plane wave = [sin(theta_y), -sin(theta_x), cos(theta_x) * cos(theta_y)]
    u = [sin(theta(2)), -sin(theta(1))*cos(theta(2)), cos(theta(1))*cos(theta(2))]; % y rotation and then x rotation

    % Get the distance version of the time delays
    for xi = 1:nx
        for yi = 1:ny
            for zi = 1:nz
                % LUT(xi, yi, zi) = dot([BFgrid.X(xi, yi, zi), BFgrid.Y(xi, yi, zi), BFgrid.Z(xi, yi, zi)], u);
                % DO I NEED TO ACCOUNT FOR THE EXTRA SHIFT THAT IS DONE TO
                % ENSURE THE ORIGINAL TIME DELAYS >= 0??? 
                % time_delays0 = time_delays0(:) - min(time_delays0(:)); % Shift so the lowest time delay is zero
                
                % if LUT(xi, yi, zi) < 0
                %     % disp([xi, yi, zi, LUT(xi, yi, zi)])
                %     disp([BFgrid.X(xi, yi, zi), BFgrid.Y(xi, yi, zi), BFgrid.Z(xi, yi, zi), LUT(xi, yi, zi)])
                % end
                LUT(xi, yi, zi) = abs(dot([BFgrid.X(xi, yi, zi), BFgrid.Y(xi, yi, zi), BFgrid.Z(xi, yi, zi)], u));
            end
        end
    end
    % LUT = LUT ./ c; % Convert from distances to time delays
    LUT = LUT ./ c + t0; % Convert from distances to time delays


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