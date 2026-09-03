function [g1_split] = vUS_2D_erf_windowed_vec_split(x, tau, k0, sigma, f_min, f_max)
%% Description:
%   Partial-volume generalization of vUS_2D_erf_vec_split.m: instead of
%   assuming the PSF samples the FULL vessel cross-section (uniform
%   velocity distribution over f in [0,2], f=2*v(r)/Vmax), this allows
%   the sampled velocity range to be a sub-window [f_min, f_max] of that
%   full range -- representing a PSF that only sees part of a vessel
%   (near-wall, near-center, or offset), rather than the whole thing.
%
%   f_min=0, f_max=2 recovers vUS_2D_erf_vec_split.m exactly (verified:
%   same prefactor 1/2, same erf arguments A/B).
%
%   f_min, f_max are FIXED inputs, not fitted parameters -- per the
%   design discussion, these should come from independent vessel
%   geometry (diameter/position estimates), not be freely optimized,
%   since this model is already weakly identified in v_xgp with just 4
%   free parameters. Adding 2 more free parameters would make that
%   worse, not better.
%
% Inputs:
%   x: [v_xgp, v_zgp, F, DC] (DC real-valued) -- same 4 parameters as
%      the full-window model
%   tau: vector of time lags [s]
%   k0: wavenumber [rad/m]
%   sigma: [sigma_x, sigma_z] [m]
%   f_min, f_max: sampled window bounds within the full [0,2] velocity
%      range (f=0: vessel wall/zero speed; f=2: vessel center/max speed)
%
% Outputs:
%   g1_split: [numel(tau), 2] = [real(g1), imag(g1)]

    v_xgp = x(1); v_zgp = x(2); F = x(3); DC = x(4);
    tau = tau(:);

    M = v_xgp.^2./sigma(1)^2 + v_zgp.^2./sigma(2)^2;
    sqrtM = sqrt(M);

    A_max = (tau.*sqrtM/2).*f_max - 2i*k0*v_zgp./sqrtM;
    A_min = (tau.*sqrtM/2).*f_min - 2i*k0*v_zgp./sqrtM;

    g1 = DC + F ./ (f_max-f_min) .* sqrt(pi)./(tau.*sqrtM) .* exp(-4*k0^2*v_zgp.^2./M) .* ...
         (erfz(A_max) - erfz(A_min));

    g1 = g1(:);
    g1_split = [real(g1), imag(g1)];
end
