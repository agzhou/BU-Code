function [Y, J] = vUS_3D_combined_batchModel(X, tau, k0)
%% Description:
%   Combined-parameter model (see vUS_3D_combined_vec.m / vUS_3D_combined_split.m) for Bc voxels at once,
%   in the plug-in form fitBatchedLM.m expects:
%       model = @(X, cols) vUS_3D_combined_batchModel(X, tau(cols), k0);
%       vi   = find(maskToUse);  nv = numel(vi);
%       x0   = [C0(vi), Vz0(vi), ones(nv,1), zeros(nv,1)];
%       lb   = [0, -50e-3, 0, 0];   ub = [Inf, 50e-3, 1, 1];
%       opts = struct('window', tau_decayed_ind(vi), 't1i', 2, 'jacobian', 'analytic', 'xscale', [1, 1e-2, 1, 1]);
%       [x, cost] = fitBatchedLM(model, x0, g1_exp{j}(vi, :), lb, ub, opts);
%   Direct analog of vUS_3D_quad_batchModel.m / vUS_3D_quad_RCA_xy_batchModel.m, but far simpler: the
%   model has no velocity-distribution integral (no k, a, no quadrature nodes) and no PSF-shape
%   parameters at all, so no moment-weighting matrix is needed -- every Jacobian column is an elementary,
%   closed-form multiple of the same [L, Bc] complex array E. There is no combined-model-specific batched
%   solver: fitBatchedLM.m / batchLM.m are model-agnostic and work with this model directly.
%
%       E = exp(-C tau^2 + i 2 k0 v_zgp tau),     Y = DC + F.*E
%       dY/dC     = F.*E.*(-tau^2)
%       dY/dv_zgp = F.*E.*(i 2 k0 tau)
%       dY/dF     = E,   dY/dDC = 1
%
%   Works for any numeric class of X (double, single, gpuArray).
%
% Inputs:
%   X: [4, Bc] parameters in PHYSICAL units, order [C, v_zgp, F, DC] (see vUS_3D_combined_vec.m)
%   tau: [L, 1] time lags [s];  k0: wavenumber [rad/m]
%
% Outputs:
%   Y: [L, Bc] complex model
%   J: [L, 4, Bc] complex Jacobian dY/dX in physical units (only computed if requested)

    Bc = size(X, 2);  L = numel(tau);
    tau = tau(:);
    C  = reshape(X(1,:),1,Bc);  vz = reshape(X(2,:),1,Bc);
    F  = reshape(X(3,:),1,Bc);  DC = reshape(X(4,:),1,Bc);

    T2 = tau.^2;                                                      % Lx1
    E  = exp(complex(-T2.*C, (2*k0*tau).*vz));                        % LxBc (the only transcendental)

    if nargout < 2
        Y = DC + F.*E;
        return
    end

    Y  = DC + F.*E;
    JC  = F.*E.*(-T2);
    Jvz = F.*E.*(1i*2*k0.*tau);
    J   = cat(2, reshape(JC,L,1,Bc), reshape(Jvz,L,1,Bc), reshape(E,L,1,Bc), ones(L,1,Bc,'like',real(E)));   % Lx4xBc
end
