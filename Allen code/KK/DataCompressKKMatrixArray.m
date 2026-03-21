% Description: Compresses data(samples, elements, TXangles) --> RawDataKK(samples, TXangles, RXangles)  
% Inputs:
%   - data: RFdata, with dimensions [# samples, # elements, # TX angles]
%   - RXangles: receive angles in radians, with dimensions [# receive angles, 2 (x and y)]
%   - s: "aspect ratio" (From TestKKadaptive.m: s=2*Pitch*SamplingFrequency/c;  %aspect ratio dx/dz)
%       - My interpretation is that s is a conversion factor between actual
%         time and samples... though maybe with another factor?
% ratio = fs/c0
% Elem: element coords with dimensions [2, total # elements] (x, y)
function RawDataKK = DataCompressKKMatrixArray(data, RXangles, ratio, Elem)
    
    % Assign parameters
    numSamples = size(data, 1);
    numElements = size(data, 2); % Total # of elements in the matrix array
    % numEl1dim = sqrt(numElements);
    disp('Note: the code accounts for only a full square matrix array')
    numTXAngles = size(data, 3);
    numRXAngles = size(RXangles, 1);
    
    % Initialize output
    RawDataKK = zeros(numSamples, numTXAngles, numRXAngles);

    nShiftAll = zeros(numElements, numRXAngles);

    % Go through and perform the shifting/basis transformation
    % dataTemp = zeros(numSamples, numElements); % Temp variable for shifting the RF Data for each TX/RX combo
    for rai = 1:numRXAngles % Receive angle index
        % u = [sin(RXangles(rai, 2)), -sin(RXangles(rai, 1))]; % Unit direction vector for the plane wave = [sin(theta_RX_y), -sin(theta_RX_x)]

        % slope = s*sin(RXangles(rai))/2; % [slope_x, slope_y]
        % slope = s .* u ./ 2;
        % Create the time delays for each element in the matrix array
        nShift = ( Elem(1, :).*sin(RXangles(rai, 2)) - Elem(2, :).*sin(RXangles(rai, 1)) .*cos(RXangles(rai, 2)) )* ratio; % Plane wave
        % disp(min(nShift(:)))
        % nShift = round(nShift(:) - min(nShift(:))); % Shift so the lowest time delay is zero
        nShift = round(nShift); % TESTING
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


    nShiftAll = reshape(nShiftAll,[16,16,numRXAngles]);
    genSliderV2(nShiftAll)

end

