% Description: initialize vUS (2D) parameters with the mesh --> R^2 method
% as in the vUS paper (the 2D analog of the v_xgp/p mesh search inside
% Jianbo's GG2vUS.m). Only v_xgp is searched here -- v_zgp (from
% findVzPhaseDiff.m), DC, and F already have trusted independent initial
% estimates and are held FIXED during this search, matching GG2vUS.m's
% own approach of holding its Ms/Mf guesses fixed while meshing vx/p.
%
% Inputs:
%   g1: [nPix, nTauFit] complex g1(tau) data (spatial dims stacked). Pass
%       already sliced to whatever tau range should be fit -- e.g.
%       tau(2:end) and g1(:,2:end), excluding tau=0 -- to match
%       g1vUS2D_residJac.m's convention: g1(tau=0)=1 is trivially
%       satisfied by every candidate model and doesn't discriminate
%       between them, so including it would just dilute R^2.
%   v_zgp0: [nPix,1] axial group velocity [m/s], e.g. from findVzPhaseDiff.m
%   DC0, F0: [nPix,1] initial DC-offset and dynamic-fraction estimates
%   PP: Processing Parameters struct (see createStruct.m). Must contain
%       either PP.k0 (angular wavenumber, rad/m) or PP.wl (wavelength, m)
%       -- same fallback convention as findVzPhaseDiff.m. PP.zp/PP.xp are
%       used, if present and consistent with size(g1,1), to spatially
%       median-filter the v_zgp0/F0 guesses before use.
%   sigma: [sigma_x, sigma_z], 1/e PSF widths [m]
%   tau: [nTauFit,1] (or [1,nTauFit]) time lag vector [s], matching g1's columns
%
% Outputs:
%   v_xgp [nVox,1] -- the v_xgp maximizing R^2 per pixel
%   R2: [nVox,1] -- the achieved R^2 (Eq. 17) at the chosen v_xgp0 per pixel
function [v_xgp, R2] = InitVx0WithMesh2D(g1, v_zgp0, DC0, F0, PP, sigma, tau)
    nPix = size(g1, 1);
    tau_mask = 2:PP.nTau;
    % Avoid NaNs on tau = 0
    tau = tau(tau_mask);
    g1 = g1(:, tau_mask);

    tau = reshape(tau, 1, []); % row, [1, nTauFit], for broadcasting against [nPix, nTauFit] data

    if isfield(PP, 'k0')
        k0 = PP.k0;
    elseif isfield(PP, 'wl')
        k0 = 2*pi/PP.wl;
    else
        error('InitvUS2DParamsWithMesh:missingWavenumber', ...
            ['PP must contain either PP.k0 (angular wavenumber, rad/m) or ', ...
             'PP.wl (wavelength, m) to build the g1 model.']);
    end

    v_zgp0 = v_zgp0(:); DC0 = DC0(:); F0 = F0(:);

    % Smooth the trusted initial guesses spatially before use (only valid
    % if PP's spatial dimensions actually match the data passed in --
    % guards against a hard reshape error if this is ever called on a
    % masked/subsetted pixel list rather than the full image)
    if nPix > 1 && isfield(PP,'zp') && isfield(PP,'xp') && PP.zp*PP.xp == nPix
        v_zgp0 = reshape( medfilt2(reshape(v_zgp0, [PP.zp, PP.xp]), [5, 5]), [nPix, 1] );
        F0     = reshape( medfilt2(reshape(F0,     [PP.zp, PP.xp]), [5, 5]), [nPix, 1] );
    end

    v_zgp = v_zgp0; % We trust the v_zgp guess a lot from the phase difference method, so don't change it
    DC = DC0;
    F = F0;

    % Mesh for v_xgp, p
    v_xgp_bounds = [0, 30].*1e-3; % Bounds for v_xgp values [m/s]
    v_xgp_step = 1*1e-3; % Step for v_xgp grid [m/s]
    v_xgp_vec = v_xgp_bounds(1):v_xgp_step:v_xgp_bounds(2);
    nVx = numel(v_xgp_vec); % # of Vx guesses, per voxel

    % R^2 denominator (Eq. 17) depends only on the observed data -- compute once
    % ydataDev = sum(abs(g1 - mean(g1, 2)).^2, 2); % [nPix,1]
    ydataDev = mean(abs(g1 - mean(g1, 2)), 2).^2; % [nPix,1]
    ydataDev(ydataDev == 0) = eps;

    bestR2 = -Inf(nPix, 1);
    v_xgp = zeros(nPix, 1);

    for vxi = 1:nVx % Loop through Vx values
        Vx_vxi = v_xgp_vec(vxi);

        g1_model = vUS_2D_erf(tau, k0, sigma, Vx_vxi, v_zgp, F, DC); % [nPix, nTauFit]

        R2_vxi = 1 - mean(abs(g1 - g1_model).^2, 2) ./ ydataDev; % [nPix,1]

        improved_voxels = R2_vxi > bestR2;
        bestR2(improved_voxels) = R2_vxi(improved_voxels);
        v_xgp(improved_voxels) = Vx_vxi;
    end

    R2 = bestR2;
end
