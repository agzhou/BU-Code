
% Description:
%   Create a struct with all the function inputs assigned as themselves in the
%   struct

function [s] = createStruct(varargin)
    ni = length(varargin); % # of inputs
    s = struct();
    for j = 1:ni
        s.(inputname(j)) = varargin{j}; % Get the names of the variables that were passed in, and create a field in the struct with it, with its associated value
    end

end