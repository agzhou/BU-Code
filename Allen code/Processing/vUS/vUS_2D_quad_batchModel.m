function [Y, J] = vUS_2D_quad_batchModel(X, tau, k0, sigma, s, w)
%% Description:
%   2D (v_xgp, v_zgp) version of vUS_3D_quad_batchModel.m: the corrected
%   Tangelder-profile g1 model, evaluated for Bc voxels at once with fixed
%   Gauss-Legendre quadrature over s = (r/R)^2, in the plug-in form that
%   fitBatchedLM.m expects:
%       model = @(X, cols) vUS_2D_quad_batchModel(X, tau(cols), k0, sigma, s, w);
%
%       g1(tau) = DC + F * int_0^1 exp( -M tau^2 f(s)^2/4 + i 2 k0 v_zgp tau f(s) ) ds
%       f(s)    = fmax * (1 - a^k s^(k/2)),   fmax = (k+2) / (k+2-2 a^k)
%       M       = v_xgp^2/sigma_x^2 + v_zgp^2/sigma_z^2
%
%   This is the model of vUS_2D_num_vec.m / vUS_2D_num_wrapper.m (they agree to
%   ~4e-8, the accuracy of integral()), with g1(0) = DC + F. Fixing k = 2 and a = 1
%   (lb = ub) gives the 4-parameter erf model of vUS_2D_erf_vec.m that
%   vUS_2D_newmodel.m fits with lsqnonlin.
%
%   Same model as the 3D version with sigma_y = sigma_x and v_tgp = v_xgp (there is
%   no y here), so it is evaluated by vUS_3D_quad_batchModel.m: one copy of the
%   Jacobian formulas. The only 2D-specific parts are the parameter naming and the
%   2-element sigma of vUS_2D_newmodel.m, which is checked -- a 3-element sigma from
%   the 3D scripts would otherwise silently put sigma_y where sigma_z belongs.
%
%   Example -- replaces the "for vi = 1:num_voxels ... lsqnonlin ..." loop of
%   vUS_2D_newmodel.m section 5, keeping its start point and bounds and freeing k and a:
%       [s, w] = gaussLegendre01(48);
%       model  = @(X, cols) vUS_2D_quad_batchModel(X, tau(cols), PP.k0, sigma, s, w);
%       vi     = find(maskToUse);  nv = numel(vi);
%       x0     = [Vx0(vi), Vz0(vi), FR0_j(vi), DCR0_j(vi), 2.5*ones(nv,1), ones(nv,1)];   % the script's x0 = [Vx0 Vz0 FR0 DCR0] + [k a]
%       lb     = [zeros(nv,1), Vz0(vi) - 0.01, zeros(nv,2), 2*ones(nv,1), zeros(nv,1)];   % the script's lb/ub + bounds on k and a
%       ub     = [30e-3*ones(nv,1), Vz0(vi) + 0.01, ones(nv,2), 3*ones(nv,1), ones(nv,1)];
%       opts   = struct('window', tau_decayed_ind(vi), 't1i', t1i, 'jacobian', 'analytic', 'xscale', [1e-2 1e-2 1 1 1 1]);
%       [x, cost] = fitBatchedLM(model, x0, g1_exp{j}(vi, :), lb, ub, opts);
%       v_xgp_stacked(vi) = x(:,1);  v_zgp_stacked(vi) = x(:,2);  F_stacked(vi) = x(:,3);
%       DC_stacked(vi)    = x(:,4);  k_stacked(vi)     = x(:,5);  a_stacked(vi) = x(:,6);
%   To reproduce the script's existing 4-parameter (erf) fit instead, fix k = 2 and a = 1:
%       x0(:,5:6) = repmat([2 1], nv, 1);  lb(:,5:6) = repmat([2 1], nv, 1);  ub(:,5:6) = repmat([2 1], nv, 1);
%   (with a Jacobian-free model, parameters fixed for every voxel are not perturbed, so this costs no extra
%   model calls for k and a).
%
%   Caveat when fitting with the batched LM solver and the script's uniform Vx0 = 5 mm/s start with lb(v_xgp) = 0:
%   the model is even in v_xgp (it enters only through v_xgp^2), so v_xgp = 0 is an exact stationary point on the
%   bound, and the LM solver can end there. On 300 synthetic voxels (sigma = [60.55 73.95] um, true v_xgp 2-20 mm/s, six free
%   parameters) 16% ended >1% worse than a per-voxel lsqnonlin, 37 of the 47 exactly at v_xgp = 0; lsqnonlin's reflective
%   trust region does not get stuck there. This is the same behaviour as vUS_3D_quad_fitBatched (the fits are identical
%   through the sigma_y = sigma_x mapping) and is not affected by tolerance (1e-14: 15%), a symmetric bound (worse) or
%   the start (multi-start over 3-4 v_xgp starts: ~9-10%; stepCap 0.1: 9%). With k = 2, a = 1 fixed it did not occur
%   (0 of 40 voxels worse than lsqnonlin on the erf model). Treat voxels with v_xgp exactly on lb as suspect. Not tested
%   on real 2D data.
%
%   Works for any numeric class of X (double, single, gpuArray).
%
% Inputs:
%   X: [6, Bc] parameters in PHYSICAL units, order [v_xgp, v_zgp, F, DC, k, a] (as in vUS_2D_num_vec.m)
%       v_xgp: x group velocity [m/s];  v_zgp: z group velocity [m/s]
%       F: dynamic fraction;  DC: static fraction (real-valued)
%       k: profile bluntness (typically 2 to 3);  a in [0, 1]: wall-speed coefficient
%   tau: [L, 1] time lags [s];  k0: wavenumber [rad/m]
%   sigma: [sigma_x, sigma_z] [m], 2 elements (the 2D convention of vUS_2D_newmodel.m)
%   s, w: Gauss-Legendre nodes [1, N] and weights [N, 1] on [0, 1], from gaussLegendre01.m
%
% Outputs:
%   Y: [L, Bc] complex g1 model
%   J: [L, 6, Bc] complex Jacobian dY/dX in physical units (only computed if requested)

    if numel(sigma) ~= 2
        error('vUS_2D_quad_batchModel:sigmaSize', ...
            ['The 2D model takes sigma = [sigma_x, sigma_z] (2 elements, as in vUS_2D_newmodel.m), got %d elements. ', ...
             'A 3D sigma = [sigma_x, sigma_y, sigma_z] belongs to vUS_3D_quad_batchModel.'], numel(sigma))
    end
    sigma3 = [sigma(1), sigma(1), sigma(2)];              % the 3D transverse-combined model with sigma_y = sigma_x, v_tgp = v_xgp
    if nargout < 2
        Y = vUS_3D_quad_batchModel(X, tau, k0, sigma3, s, w);
    else
        [Y, J] = vUS_3D_quad_batchModel(X, tau, k0, sigma3, s, w);
    end
end
