% Description: find the (first) time lag at which the field autocorrelation
% function's magnitude has decayed to its steady-state value

% Inputs:
%   g1: [nVoxels, nTau] complex field autocorrelation function
%   tau: [nTau, 1] time lag vector
%   absolute_tau_ss_cutoff_s: [scalar, in seconds] safe guess for when the steady state
%                           has already occurred
%   tau1_ind: [scalar] index at which tau=tau1 occurs
%   too_fast_decay_s: [scalar] if the |g1| decays before this time lag,
%                     assume the voxel is noise (not physiologically
%                     possible to have flow that fast). This should be in
%                     the same units as the tau vector.

% Outputs:
%   tau_decayed_ind: [nVoxels, 1] (first) time lag index vector where the g1 magnitude has
%                    likely finished its decay
%   voxel_quality: [nVoxels, 1] boolean value per voxel if the decay is too
%                  fast to be reasonable: true if voxel is good, false is
%                  voxel is bad.

function [tau_decayed_ind, voxel_quality] = findTauDecayed(g1, tau, tau1_ind, absolute_tau_ss_cutoff_s, too_fast_decay_s)
    % Crop data to consider only tau1:nTau
    tau = tau(tau1_ind:end);
    g1 = g1(:, tau1_ind:end);

    if size(g1, 2) == 1
        warning('Check the g1 input: should be of size [nVoxels, nTau]')
    end
    numVoxels = size(g1, 1);

    % Convert cutoffs/thresholds to indices
    absolute_tau_ss_cutoff_ind = find(tau >= absolute_tau_ss_cutoff_s, 1, 'first');
    too_fast_decay_s_ind = find(tau >= too_fast_decay_s, 1, 'first');

    nTau = length(tau); % # of time lags provided
    ag1 = abs(g1); % Get magnitude
    ag1_smoothed = movmean(ag1, 3, 2); % Smooth the |g1| over time a little

    % First metric: when the |g1| reaches the "steady-state" value for the
    %               first time
    tau_ss = max(round(nTau*2/3), absolute_tau_ss_cutoff_ind); % A safe guess for when the g1 has approached steady state
    ag1_ss = median(ag1(:, tau_ss:end), 2); % Define the steady-state value for each voxel as the median over some tau range

    tau_decayed_ind_1 = ones(numVoxels, 1);
    for vi = 1:numVoxels
        tau_decayed_ind_1(vi) = find(ag1_smoothed(vi, :) < ag1_ss(vi), 1, 'first');
    end
    
    % Second metric: when the |g1| decays to 10% of its initial value at
    %                tau = tau1
    ag1_decay_threshold = 0.1 .* ag1(:, 1); % tau1 is at index 1, since it was all cropped earlier in the function
    tau_decayed_ind_2 = ones(numVoxels, 1);
    for vi = 1:numVoxels
        tau_decayed_ind_2(vi) = find(ag1_smoothed(vi, :) < ag1_decay_threshold(vi), 1, 'first');
    end

    tau_decayed_ind = min(tau_decayed_ind_1, tau_decayed_ind_2); % Take the minimum of guesses #1 and #2

    % Check for bad voxels --> the decay is too fast
    voxel_quality = tau_decayed_ind >= too_fast_decay_s_ind;

    % Compensate for the tau cropping
    tau_decayed_ind = tau_decayed_ind + (tau1_ind - 1);
    
end