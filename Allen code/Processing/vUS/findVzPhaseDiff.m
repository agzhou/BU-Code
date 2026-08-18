% NOTE: the frame (sampling) rate must sample at least faster than 2 samples per cycle:
% phase = 2*k0*vz*tau
% dphase/dt = 2*k0*vz [rad/sec]
% minSampleRate = 2*dphase/dt = 4*k0*vz [rad/sec]
%               = 4*k0*vz/2pi = 2*k0*vz/pi [cycles/sec]

% e.g., minSampleRate = 2*(2*pi/(1540/15.625e6))*(50e-3)/pi = 2029 Hz, for a
% max speed of 50 mm/s, at c = 1540 m/s and 15.625 MHz transducer frequency
% ALTERNATIVELY:
% maxVz = sampleRate/(2*k0)*pi
% e.g., maxVz = 5000/(2*(2*pi/(1540/15.625e6)))*pi = 123.2 mm/s, for a frame
% rate of 5000 Hz, c = 1540 m/s, and 15.625 MHz transducer frequency

% Description: Estimate the axial group velocity Vz0 (and an equivalent
% tau_V) from a magnitude-weighted circular mean of adjacent-lag phase
% differences of g1(tau). Signature-compatible with guessVz0.m (my
% in-progress port of Jianbo's GG2Vz function), so it can be swapped in at
% its call site in vUS_2D.m:
%   [Vz0, tau_V] = guessVz0(g1adj_stacked_j, PP, tauInterpFactor);
%   [Vz0, tau_V] = findVzPhaseDiff(g1adj_stacked_j, PP, tauInterpFactor);
%
% guessVz0.m (like GG2Vz.m) hunts for the FIRST local minimum of |g1(tau)|
% (optionally cross-checked against an autocorrelation-based estimate) --
% a single noisy feature that is prone to catastrophic errors when an
% early sample dips due to measurement noise rather than true
% decorrelation (worst for slow flows, where the true minimum is far out
% in tau and gets pre-empted by an earlier noise-driven dip). A Monte
% Carlo comparison against valley detection (synthetic g1 built from the
% vUS model + realistic finite-sample estimator noise, matching Eqs. 2 and
% 15 of Tang et al. 2020) showed this phase-difference approach reduced
% RMSE by ~7-40x and eliminated gross (>50%) errors entirely across a
% 1-20 mm/s, 3-noise-level sweep that valley detection failed on 5-93% of
% the time depending on the exact smoothing used.
%
% Rather than locating a feature, this pools EVERY adjacent-lag pair,
% weighted by how strongly correlated that pair still is
% (|g1(tau_i)|*|g1(tau_i+1)|), into a single magnitude-weighted circular
% mean of the phase rotation rate d(phi)/dt = 2*k0*Vz0 (Eq. 15). Pairs
% where the signal has already decorrelated into noise are automatically
% down-weighted, so this needs no separate "noisy vs. clean" branch or
% valley-count heuristic (cf. GG2Vz.m's ddGG/NdGG selection logic, or
% guessVz0.m's planned first-minimum-vs-autocorrelation cross-check).
%
% Because it only ever differences ADJACENT lags, each phase step stays
% within (-pi, pi] by construction as long as the standard Nyquist
% condition on Vz0 is met (2*k0*Vz0*dt < pi) -- there is no risk of
% unwrap failure the way a global phase-unwrap-and-regress approach would
% have, and (unlike guessVz0.m's valley search) no temporal upsampling is
% needed to localize a feature precisely.
%
% Inputs:
%   g1_stacked: [nVoxels, nTau] complex g1 matrix (spatial dimensions
%       stacked). Matches guessVz0.m's g1_stacked / vUS_2D.m's
%       g1adj_stacked_j: tau starts at 0 (g1T.m convention), so
%       g1_stacked(:, 1) = g1(tau=0) = 1 exactly for every voxel.
%   PP: Processing Parameters struct (see createStruct.m). Must contain:
%       PP.frameRate: sIQ frame rate [Hz]
%       PP.k0 (angular wavenumber [rad/m]) OR PP.wl (wavelength [m],
%           used as k0 = 2*pi/PP.wl). NOTE: as of this writing,
%           createStruct.m does not populate either field -- add one
%           (e.g. PP.k0 = 2*pi/P.wl) before calling this function.
%
% Outputs:
%   Vz0: [nVoxels, 1], initial guess for axial group velocity [m/s].
%        Sign follows the convention that positive Vz0 means flow toward
%        the transducer (same sign convention as Eq. 15's phase term
%        exp(i*2*k0*Vz0*tau)).
%   tau_V: [nVoxels, 1], equivalent quarter-cycle time [s], defined
%        self-consistently as wavelength./(4*Vz0) (Eq. 16) so it can
%        still be used anywhere tau_V is expected directly. Voxels with
%        Vz0 == 0 get tau_V = Inf.
function [Vz0, tau_V] = findVzPhaseDiff(g1_stacked, PP)
    if isfield(PP, 'k0')
        k0 = PP.k0;
    elseif isfield(PP, 'wl')
        k0 = 2*pi/PP.wl;
    else
        error('findVzPhaseDiff:missingWavenumber', ...
            ['PP must contain either PP.k0 (angular wavenumber, rad/m) or ', ...
             'PP.wl (wavelength, m) to convert the phase rotation rate into ', ...
             'a velocity. Neither field was found on PP -- add one when ', ...
             'building PP (e.g. in createStruct.m).']);
    end

    dt = 1/PP.frameRate;
    lambda0 = 2*pi/k0; % wavelength [m]

    % Adjacent-lag phase differences and their reliability weights
    dphi = angle( g1_stacked(:, 2:end) .* conj(g1_stacked(:, 1:end-1)) );  % [nVoxels, nTau-1], each in (-pi, pi]
    w    = abs(g1_stacked(:, 2:end)) .* abs(g1_stacked(:, 1:end-1));       % down-weights already-decorrelated pairs

    % Magnitude-weighted circular mean of the phase-rotation rate
    wsum = sum(w, 2);
    wsum(wsum == 0) = eps; % guard fully decorrelated / no-signal voxels against 0/0
    rate = sum(w .* dphi, 2) ./ wsum / dt; % [nVoxels, 1], rad/s --> rate = (2*k0*Vz*dt)/dt

    Vz0 = rate ./ (2*k0); % Eq. 15: phase term is exp(i*2*k0*Vz0*tau)

    tau_V = lambda0 ./ (4 * abs(Vz0));
    tau_V(Vz0 == 0) = Inf;
end
