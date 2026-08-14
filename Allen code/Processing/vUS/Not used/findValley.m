% Description: Find the first local minimum in some data
% Data should be of shape [nVoxels, nTau]
function [value, ind] = findValley(data)
    diff_data = diff(data, 1, 2);
    sign_diff = sign(diff_data.').' == 1;
    [value, ind] = max(sign_diff, [], 2);

end