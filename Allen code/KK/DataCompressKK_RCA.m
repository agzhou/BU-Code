% Description: Compresses data(samples, elements, TXangles) --> RawDataKK(samples, TXangles, RXangles)  
% Inputs:
%   - data: RFdata, with dimensions [# samples, # elements, total # TX angles]
%   - RXangles: receive angles in radians, with dimensions [total # receive angles, 1]
%   - ratio = fs/c0 (sampling frequency / speed of sound in medium)
%   - ElemPos: element coords with dimensions [2, total # elements] (x, y)
%   - time_delays_TX: [# elements, # TX angles] matrix of time delays [s]

function RawDataKK = DataCompressKK_RCA(data, anglesRX, ratio, ElemPos, time_delays_TX, RF_fs)
    
    % Assign parameters
    numSamples = size(data, 1);
    numElements = size(data, 2); % Total # of elements in the RCA (# rows + # columns)
    % numEl1dim = sqrt(numElements);
    % disp('Note: the code accounts for only a full square matrix array')
    numTXAngles = size(data, 3);
    numRXAngles = size(anglesRX, 1);

    % Initialize output
    RawDataKK = zeros(numSamples, numTXAngles, numRXAngles);

    % nShiftAll = zeros(numElements, numRXAngles);

    % Go through and perform the shifting/basis transformation
    % dataTemp = zeros(numSamples, numElements); % Temp variable for shifting the RF Data for each TX/RX combo
    
    % **** NEED TO BREAK THIS UP INTO RC AND CR **** %

    for rai = 1:numRXAngles % Receive angle index

        % Create the time delays for each element in the matrix array
        nShift = round(( ElemPos(1, :).*sin(anglesRX(rai, 2)) - ElemPos(2, :).*sin(anglesRX(rai, 1)) .*cos(anglesRX(rai, 2)) )* ratio); % Plane wave

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

