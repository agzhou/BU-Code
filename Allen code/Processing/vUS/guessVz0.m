% Description: calculate an initial guess for the axial component of group
%              velocity (Eq. 16 in the vUS paper)

% Inputs:
%   g1_stacked: [nVoxels, nTau] complex g1 matrix (spatial dimensions stacked)
%   PP: Processing Parameters struct
%   tauInterpFactor: factor by which to upsample time lags

% Outputs:
%   Vz0: initial guesses for v_zgp [nVox, 1]
%   tau_V_halfcycle: tau_V (half-cycle version) [nVox, 1]
function [Vz0, tau_V_halfcycle] = guessVz0(g1_stacked, PP, tauInterpFactor)
    % Set smoothing parameters
    smoothing_window_ms = 10; % Set a time-based [milliseconds] smoothing window (see upsampling step below)
    smoothing_window_samples = round(PP.frameRate*smoothing_window_ms/1e3); % Convert the above to samples, according to the acquisition frame rate [Hz]
    
    % Set guardrail parameter
    % guard_index = 5/5000*PP.frameRate; % 1 ms, converted to index, based on frame rate
    guard_time = 5/5000; % 1 ms, converted to index, based on frame rate

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
    [value1, index1] = findValley(g1US_R);
    tau_V_halfcycle_v1 = index1 ./ PP.frameRate; % Convert the indices for each local min to time lags [s]
    tau_V_halfcycle_v1(tau_V_halfcycle_v1 < guard_time) = tauUS(end); % If the minimum was detected too early (e.g., within the first 1 ms, or whatever guard_index is set to), set tau_V_halfcycle to be the last possible time lag

    % First tau_V estimate: autocorrelation of the smoothed, original g1's real component
    g1R_AC = g1T(movmean(real(g1_stacked), 5, 2), PP.nTau);
    [value2, index2] = findValley(g1R_AC);
    tau_V_halfcycle_v2 = index2 ./ PP.frameRate; % Convert the indices for each local min to time lags [s]
    tau_V_halfcycle_v2(value2 == 0) = tauUS(end); % If the minimum was detected at a time lag where the autocorrelation = 0, set tau_V_halfcycle to be the last possible time lag

    % Choose between v1 and v2
    diff2 = diff( (sign(diff(g1US_R, 1, 2)).' == 1).', 1, 2);
    diff2(diff2 < 1) = 0; % only valleys
    numWiggles = sum(diff2, 2); % Total # of valleys for each voxel
    % How many valleys (wiggles) are present --> a measure of noisiness (more = more noisy)

    wiggleThreshold = 5; % Threshold for # of wiggles (less than this: use the first valley method. More: use the autocorrelation method)
    tau_V_halfcycle = (numWiggles < wiggleThreshold) .* tau_V_halfcycle_v1 + (numWiggles >= wiggleThreshold) .* tau_V_halfcycle_v2;
    
    Vz0 = PP.wl ./ (4.*tau_V_halfcycle); % Apply Eq. 16 from the vUS paper
    Vz0 = Vz0 .* sign( mean(g1US_I(:, 1:3), 2) ); % Multiply v_zgp by its sign (average the first few values of imag(g1) to get the direction of the spiral)
    

end
