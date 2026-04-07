% Description: Compresses data(samples, elements, TXangles) --> RawDataKK(samples, TXangles, RXangles)  
% Inputs:
%   - data: RFdata, with dimensions [# samples, # elements, total # TX angles]
%   - RXangles: receive angles in radians, with dimensions [total # receive angles, 1]
%   - ratio = fs/c0 (sampling frequency / speed of sound in medium)
%   - ElemPos: element coords with dimensions [2, total # elements] (x, y)
%   - time_delays_TX: [total # TX angles, total # elements] matrix of time delays. Units are [samples]

% Use this with BeamformKK_RCA.m

function RawDataKK = DataCompressKK_RCA(data, anglesRX, ratio, ElemPos, time_delays_TX)
    
    % Assign parameters
    numSamples = size(data, 1);
    numElements = size(data, 2); % # of elements in one dimension of the RCA (# rows = # columns)
    % numEl1dim = sqrt(numElements);
    % disp('Note: the code accounts for only a full square matrix array')
    numTXAngles = size(data, 3);
    numRXAngles = size(anglesRX, 1);

    % Initialize output
    RawDataKK = zeros(numSamples, numTXAngles, numRXAngles);

    % nShiftAll = zeros(numElements, numRXAngles);

    % Go through and perform the shifting/basis transformation
    % dataTemp = zeros(numSamples, numElements); % Temp variable for shifting the RF Data for each TX/RX combo
    
    % **** BREAK UP ANGLES INTO RC AND CR **** %
    inds_CR = 1:numRXAngles/2;
    inds_RC = numRXAngles/2 + 1:numRXAngles;
    anglesRX_CR = anglesRX(inds_CR, 2);
    anglesRX_RC = anglesRX(inds_RC, 1);

    elem_inds_CR = 1:numElements;
    elem_inds_RC = numElements + 1:numElements*2;
    ElemPos_CR = ElemPos(1, elem_inds_CR);
    ElemPos_RC = ElemPos(2, elem_inds_RC);

    % Shift the input RF data to correct for the unequal time delays bias terms
    % caused by not being able to do negative time delays
    % TX_shift_compensation = numSamples - round(max(time_delays_TX, [], 1) * RF_fs / 2); % [samples] Note: subtracting the bias terms from the total numSamples allows us to shift backwards with circshift, which only takes in positive numbers for the shifts
    data_nonshifted = data;
    TX_shift_compensation = numSamples - round( (max(time_delays_TX, [], 2) - min(time_delays_TX, [], 2)) / 2); % [samples] Note: subtracting the bias terms from the total numSamples allows us to shift backwards with circshift, which only takes in positive numbers for the shifts. This version accounts for an additional delay past zero, e.g., a positive startDepth
    for tai = 1:numTXAngles
        data(:, :, tai) = circshift(data_nonshifted(:, :, tai), TX_shift_compensation(tai), 1);
    end

    for rai = inds_CR % Receive angle index
    
        % Create the time delays for each element in the RCA
        % nShift = ( ElemPos(1, :).*sin(RXangles(rai, 2)) - Elem(2, :).*sin(RXangles(rai, 1)) .*cos(RXangles(rai, 2)) )* ratio; % Plane wave
        nShift = round(ElemPos_CR .* sin(anglesRX_CR(rai)) * ratio);
        % nShift = nShift - min(nShift);
        nShiftAll(:, rai) = nShift;

        for tai = 1:numTXAngles % Transmit angle index
               dataTemp = zeros(numSamples, numElements); % Temp variable for shifting the RF Data for each TX/RX combo

               for ei = 1:numElements % Element index
                   % nShift = zeros(size(slope)); % How many samples to shift by for element ei
                   % 
                   %     if slope(dim) > 0
                   %         nShift(dim) = ((ei - 1) * abs(slope(dim)));
                   %     else
                   %         % nShift(dim) = ((numElements - ei)*abs(slope)); % old code 
                   %         nShift(dim) = ((numEl1dim - ei) * abs(slope(dim)));
                   %     end
                    % nShift = time_delays(ei);
                   
                   
                   % **** BELOW IS STILL UNCHANGED **** %
                   dataTemp(:, ei) = circshift( data(:, ei, tai), nShift(ei) );


               end 
               RawDataKK(:, tai, rai) = sum(dataTemp, 2);
                
        end

    end



    for rai = inds_RC % Receive angle index
    
        % Create the time delays for each element in the RCA
        % nShift = round(( ElemPos(1, :).*sin(anglesRX(rai, 2)) - ElemPos(2, :).*sin(anglesRX(rai, 1)) .*cos(anglesRX(rai, 2)) )* ratio); % Plane wave
        nShift = round(ElemPos_RC .* sin(anglesRX_RC(rai - numRXAngles/2)) * ratio);
        % nShift = nShift - min(nShift);
        nShiftAll(:, rai) = nShift;

        for tai = 1:numTXAngles % Transmit angle index
               dataTemp = zeros(numSamples, numElements); % Temp variable for shifting the RF Data for each TX/RX combo

               for ei = 1:numElements % Element index
                   % nShift = zeros(size(slope)); % How many samples to shift by for element ei
                   % 
                   %     if slope(dim) > 0
                   %         nShift(dim) = ((ei - 1) * abs(slope(dim)));
                   %     else
                   %         % nShift(dim) = ((numElements - ei)*abs(slope)); % old code 
                   %         nShift(dim) = ((numEl1dim - ei) * abs(slope(dim)));
                   %     end
                    % nShift = time_delays(ei);
                   
                   
                   % **** BELOW IS STILL UNCHANGED **** %
                   dataTemp(:, ei) = circshift( data(:, ei, tai), nShift(ei) );


               end 
               RawDataKK(:, tai, rai) = sum(dataTemp, 2);
                
        end

    end

end

