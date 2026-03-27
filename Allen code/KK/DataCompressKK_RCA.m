% Description: Compresses data(samples, elements, TXangles) --> RawDataKK(samples, TXangles, RXangles)  
% Inputs:
%   - data: RFdata, with dimensions [# samples, # elements, # TX angles]
%   - RXangles: receive angles in radians, with dimensions [# receive angles, 2 (x and y)]
%   - ratio = fs/c0 (sampling frequency / speed of sound in medium)
%   - Elem: element coords with dimensions [2, total # elements] (x, y)
%   - time_delays_TX: [# elements, # TX angles] matrix of time delays [s]
function RawDataKK = DataCompressKK_RCA(data, RXangles, ratio, Elem, time_delays_TX, RF_fs)
    
    % Assign parameters
    numSamples = size(data, 1);
    numElements = size(data, 2); % Total # of elements in the matrix array
    % numEl1dim = sqrt(numElements);
    disp('Note: the code accounts for only a full square matrix array')
    numTXAngles = size(data, 3);
    numRXAngles = size(RXangles, 1);
    data_nonshifted = data;

    % Shift the input RF data to correct for the unequal time delays bias terms
    % caused by not being able to do negative time delays
    % TX_shift_compensation = numSamples - round(max(time_delays_TX, [], 1) * RF_fs / 2); % [samples] Note: subtracting the bias terms from the total numSamples allows us to shift backwards with circshift, which only takes in positive numbers for the shifts
    TX_shift_compensation = numSamples - round( (max(time_delays_TX, [], 1) - min(time_delays_TX, [], 1)) * RF_fs / 2); % [samples] Note: subtracting the bias terms from the total numSamples allows us to shift backwards with circshift, which only takes in positive numbers for the shifts. This version accounts for an additional delay past zero, e.g., a positive startDepth
    for tai = 1:numTXAngles
        data(:, :, tai) = circshift(data_nonshifted(:, :, tai), TX_shift_compensation(tai), 1);
    end
    

    % Initialize output
    RawDataKK = zeros(numSamples, numTXAngles, numRXAngles);

    % nShiftAll = zeros(numElements, numRXAngles);

    % Go through and perform the shifting/basis transformation
    % dataTemp = zeros(numSamples, numElements); % Temp variable for shifting the RF Data for each TX/RX combo
    for rai = 1:numRXAngles % Receive angle index
        % u = [sin(RXangles(rai, 2)), -sin(RXangles(rai, 1))]; % Unit direction vector for the plane wave = [sin(theta_RX_y), -sin(theta_RX_x)]

        % slope = s*sin(RXangles(rai))/2; % [slope_x, slope_y]
        % slope = s .* u ./ 2;
        % Create the time delays for each element in the matrix array
        nShiftNoCompensation = ( Elem(1, :).*sin(RXangles(rai, 2)) - Elem(2, :).*sin(RXangles(rai, 1)) .*cos(RXangles(rai, 2)) )* ratio; % Plane wave
        % disp(min(nShift(:)))
        % nShift = round(nShiftNoCompensation(:) - min(nShiftNoCompensation(:))); % Shift so the lowest time delay is zero (% Compensate for not being able to have negative shifts since the origin is at the center of the probe)
        % nShiftRX_compensation = round(nShiftNoCompensation(:) - nShift);
        nShift = round(nShiftNoCompensation); % TESTING
        nShiftAll(:, rai) = nShift;

        for tai = 1:numTXAngles % Transmit angle index
                dataTemp = zeros(numSamples, numElements); % Temp variable for shifting the RF Data for each TX/RX combo

               for ei = 1:numElements % Element index
                   % nShift = zeros(size(slope)); % How many samples to shift by for element ei
                   % 
                   % 
                   % % fix below (replace with the time_delay generation)
                   % for dim = 1:2
                   %     if slope(dim) > 0
                   %         nShift(dim) = ((ei - 1) * abs(slope(dim)));
                   %     else
                   %         % nShift(dim) = ((numElements - ei)*abs(slope)); % old code 
                   %         nShift(dim) = ((numEl1dim - ei) * abs(slope(dim)));
                   %     end
                   % end
                    % nShift = time_delays(ei);
                   
                   
                   % **** BELOW IS STILL UNCHANGED **** %
                   % dataTemp(:, ei) = circshift( data(:, ei, tai), nShift(ei) );

                   if nShift(ei) < 0
                       numShift = numSamples + nShift(ei);
                   else
                       numShift = nShift(ei);
                   end
                   % disp(numShift)
                   dataTemp(:, ei) = circshift( data(:, ei, tai), numShift ); 

               end 
               
               RawDataKK(:, tai, rai) = sum(dataTemp, 2);
                
        end

    end


    % nShiftAll = reshape(nShiftAll,[16,16,numRXAngles]);
    % genSliderV2(nShiftAll)

end

