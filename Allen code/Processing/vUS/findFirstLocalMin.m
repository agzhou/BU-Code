% Find the index of the first local min of the |vectorized data g1v|
%   g1v: [# voxels, # time lags (tau)]
%   tau: 1 x # time lags (tau) vector

%   tau_first_min_index: index at the first local min

% [tau_first_min_index, tau_first_min]
function tau_first_min_index = findFirstLocalMin(g1v, numg1pts)
    ag1v = abs(g1v); 
%     ag1v_diff = diff(ag1v, 1, 2); % Find the local min by seeing when the diff changes sign
%     [test, tau_index] = find(ag1v_diff > 0, size(ag1v_diff, 1), 'first');
%     tau_matrix = repmat(tau, size(ag1v, 1), 1);
    tau_index_matrix = repmat(1:numg1pts, size(ag1v, 1), 1);
    TF = islocalmin(ag1v, 2, 'MaxNumExtrema', 1); % Get the first local minimum of |g1| per vectorized voxel
%     tau_first_min = tau_matrix(TF);

    % Account for voxels where there is no local minimum (|g1| keeps decreasing)
    bad_voxels = find(sum(TF, 2) == 0);

    TF_corrected = TF;
%     TF_corrected(bad_voxels, :) = NaN(1, numg1pts);
    TF_corrected(bad_voxels, :) = repmat([zeros(1, numg1pts - 1), 1], length(bad_voxels), 1); % Set the initial guess for tau_V (the first min) equal to the last point for these bad voxels




    tau_first_min_index = tau_index_matrix(TF_corrected);
end