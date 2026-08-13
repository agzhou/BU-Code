% Description: calculate an initial guess for the axial component of group
%              velocity (Eq. 16 in the vUS paper)

% Inputs:
%   g1_stacked: [nVoxels, nTau] complex g1 matrix (spatial dimensions stacked)
%   PP: Processing Parameters struct
%   tauInterpFactor: factor by which to upsample time lags

function [Vz0, tau_V] = guessVz0(g1_stacked, PP, tauInterpFactor)
    % Set smoothing parameters
    smoothing_window_ms = 10; % Set a time-based [milliseconds] smoothing window (see upsampling step below)
    smoothing_window_samples = round(PP.frameRate*smoothing_window_ms/1e3); % Convert the above to samples, according to the acquisition frame rate [Hz]
    
    % Define time lag vectors
    tau = (0:PP.nTau - 1)' ./ PP.frameRate; % [nTau, 1] time lag vector (seconds)
    tauUS = linspace(0, (PP.nTau - 1)./ PP.frameRate, PP.nTau*tauInterpFactor).'; % UpSampled time lag vector (seconds)

    % Upsample g1 in time (and then smooth in time)
    g1US_R = movmean(interp1(tau, real(g1_stacked).', tauUS, 'linear').', smoothing_window_samples, 2); % Real component of the linearly UpSampled g1 [nVoxels, nTau]
    g1US_I = movmean(interp1(tau, imag(g1_stacked).', tauUS, 'linear').', smoothing_window_samples, 2); % Real component of the linearly UpSampled g1 [nVoxels, nTau]
    % ind = sub2ind([PP.zp, PP.xp], 60, 136);
    % figure; plot(tau, abs(g1_stacked(ind, :)), tauUS, abs(complex(g1US_R(ind, :), g1US_I(ind, :))))
    % figure; plot(tau, real(g1_stacked(ind, :)), tauUS, g1US_R(ind, :))
    % figure; plot(tau, imag(g1_stacked(ind, :)), tauUS, g1US_I(ind, :))

    % First tau_V estimate: the first minimum of the upsampled g1's real component
    % tau_V_halfcycle_v1 = islocalmin(g1US_R, 2); % Get indices of the first valley in the real component
    %     % Account for voxels where there is no local minimum (|g1| keeps decreasing)
    % bad_voxels = find(sum(tau_V_halfcycle_v1, 2) == 0);
    % tau_V_halfcycle_v1(bad_voxels, :) = repmat([zeros(1, PP.nTau - 1), 1], length(bad_voxels), 1); % Set the initial guess equal to the last point for these bad voxels
    diff_g1US_R = diff(g1US_R, 1, 2); % First-order diff along the time dimension
    test2 = sign(diff_g1US_R) == 1;
    % tau_V_halfcycle_v1

    % First tau_V estimate: autocorrelation of the smoothed, original g1's real component
    g1R_AC = g1T(movmean(real(g1_stacked), 5, 2), PP.nTau);
    [value, ind] = findValley(g1R_AC);
    % tau_V_halfcycle_v2

    % g1_R_halfcycle_v1



end
