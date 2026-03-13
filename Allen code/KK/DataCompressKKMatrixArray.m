% Description: Compresses data(samples, elements, TXangles) --> RawDataKK(samples, TXangles, RXangles)  
% Inputs:
%   - s: "aspect ratio" (From TestKKadaptive.m: s=2*Pitch*SamplingFrequency/c;  %aspect ratio)
%       - My interpretation is that s is a conversion factor between actual
%         time and samples... though maybe with another factor?
%   - data: RFdata, with dimensions [# samples, # elements, # TX angles]

function RawDataKK = DataCompressKKMatrixArray(data, RXangles, s)
    
    % Assign parameters
    numSamples = size(data, 1);
    numElements = size(data, 2); % Total # of elements in the matrix array
    numEl1dim = sqrt(numElements);
    disp('Note: the code accounts for only a full square matrix array')
    numTXAngles = size(data, 3);
    numRXAngles = size(RXangles, 1);
    
    % Initialize outputs and temp variables
    RawDataKK = zeros(numSamples, numTXAngles, numRXAngles);
    dataTemp = zeros(numSamples, numElements);
    
    for rai = 1:numRXAngles % Receive angle index
        
        slope = s*sin(RXangles(rai))/2; % [slope_x, slope_y]
    
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
                   
                   dataTemp(:,ei)=circshift(data(:,ei,tai),round(nShift));
               end 
               
               RawDataKK(:,tai,rai)=sum(dataTemp,2);
                
        end

    end
end

