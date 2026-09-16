% Description: find the (first) time lag at which the field autocorrelation
% function's magnitude has decayed to its steady-state value

% Inputs:
%   g1: [nVoxels, nTau] complex field autocorrelation function
%   tau: [nTau, 1] time lag vector
%   absolute_tau_ss_cutoff: [scalar] safe guess for when the steady state
%                           has already occurred
% Outputs:
%   tau_decayed_ind: [scalar] (first) time lag index where the g1 magnitude has
%   likely decayed

function tau_decayed_ind = findTauDecayed(g1, tau, absolute_tau_ss_cutoff)
    nTau = length(tau); % # of time lags provided
    ag1 = abs(g1); % Get magnitude
    ag1_smoothed = movmean(ag1, 3, 2); % Smooth the |g1| over time a little
    tau_ss = max(round(nTau*2/3), absolute_tau_ss_cutoff); % A safe guess for when the g1 has approached steady state
    ag1_ss = median(ag1(:, tau_ss:end)); % Define the steady-state value as the median over some tau range

    tau_decayed_ind = 
    
end