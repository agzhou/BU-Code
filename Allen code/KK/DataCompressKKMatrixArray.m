% Description: Compresses data(samples, elements, TXangles) --> RawDataKK(samples, TXangles, RXangles)  
% Inputs:
%   - data: RFdata, with dimensions [# samples, # elements, # TX angles]
%   - RXangles: receive angles in radians, with dimensions [# receive angles, 2 (x and y)]
%   - s: "aspect ratio" (From TestKKadaptive.m: s=2*Pitch*SamplingFrequency/c;  %aspect ratio)
%       - My interpretation is that s is a conversion factor between actual
%         time and samples... though maybe with another factor?

function RawDataKK = DataCompressKKMatrixArray(data, RXangles, s)
    
    % Assign parameters
    numSamples = size(data, 1);
    numElements = size(data, 2); % Total # of elements in the matrix array
    numEl1dim = sqrt(numElements);
    disp('Note: the code accounts for only a full square matrix array')
    numTXAngles = size(data, 3);
    numRXAngles = size(RXangles, 1);
    
    % Initialize output
    RawDataKK = zeros(numSamples, numTXAngles, numRXAngles);
    
    % Go through and perform the shifting/basis transformation
    dataTemp = zeros(numSamples, numElements); % Temp variable for shifting the RF Data for each TX/RX combo
    for rai = 1:numRXAngles % Receive angle index
        u = [sin(RXangles(rai, 2)), -sin(RXangles(rai, 1))]; % Unit direction vector for the plane wave = [sin(theta_RX_y), -sin(theta_RX_x)]

        % slope = s*sin(RXangles(rai))/2; % [slope_x, slope_y]
        slope = s .* u ./ 2;
    
        for tai = 1:numTXAngles % Transmit angle index

               for ei = 1:numElements % Element index
                   nShift = zeros(size(slope)); % How many samples to shift by for element ei
                   for dim = 1:2
                       if slope(dim) > 0
                           nShift(dim) = ((ei - 1) * abs(slope(dim)));
                       else
                           % nShift(dim) = ((numElements - ei)*abs(slope)); % old code 
                           nShift(dim) = ((numEl1dim - ei) * abs(slope(dim)));
                       end
                   end
                   
                   % **** BELOW IS STILL UNCHANGED **** %
                   dataTemp(:, ei) = circshift( data(:, ei, tai), round(nShift) );
               end 
               
               RawDataKK(:, tai, rai) = sum(dataTemp,2);
                
        end

    end
end

