% Description: initialize vUS (2D) parameters with the mesh --> R^2 method
% as in the vUS paper (the 2D analog of the v_xgp/p mesh search inside
% Jianbo's GG2vUS.m). Only v_xgp and p are searched here -- v_zgp (from
% findVzPhaseDiff.m), DC, and F already have trusted independent initial
% estimates and are held FIXED during this search, matching GG2vUS.m's
% own approach of holding its Ms/Mf guesses fixed while meshing vx/p.
%
% Efficiency note: naively broadcasting a [nPix, nVxCandidates,
% nPCandidates, nTau] model tensor for a full image is prohibitively
% large (e.g. ~613M elements / ~10GB as complex doubles for a 20,000-pixel
% image with the default grid below) -- this instead loops over the
% (v_xgp, p) mesh points (a few hundred, not per-pixel) and reduces each
% candidate's [nPix, nTau] model curve down to a per-pixel R^2
% immediately, keeping a running best-R^2 (and the (v_xgp, p) achieving
% it) per pixel rather than ever storing the full tensor. Peak memory is
% O(nPix x nTau), the same order as the input g1 data itself, regardless
% of how fine the v_xgp/p grid is.
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
%   v_zgp: [nPix,1] -- v_zgp0, median-filtered (if applicable); unchanged
%          otherwise. This function does not search over v_zgp.
%   v_xgp, p: [nPix,1] -- the (v_xgp, p) mesh point maximizing R^2 per pixel
%   DC, F: [nPix,1] -- DC0/F0 passed through (F is median-filtered like v_zgp0)
%   R2: [nPix,1] -- the achieved R^2 (Eq. 17) at the chosen (v_xgp, p) per pixel
function [v_zgp, v_xgp, p, DC, F, R2] = InitvUS2DParamsWithMesh(g1, v_zgp0, DC0, F0, PP, sigma, tau)
    nPix = size(g1, 1);
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
    p_vec = linspace(0, 1, 10); % p vector of values from 0-1
    nVx = numel(v_xgp_vec); nP = numel(p_vec);

    % R^2 denominator (Eq. 17) depends only on the observed data -- compute once
    ydataDev = sum(abs(g1 - mean(g1,2)).^2, 2); % [nPix,1]
    ydataDev(ydataDev == 0) = eps;

    bestR2 = -inf(nPix,1);
    v_xgp = zeros(nPix,1);
    p = zeros(nPix,1);

    zVec = v_zgp.*tau; % [nPix, nTauFit], reused across every (v_xgp,p) candidate
    envelopeZphase = exp(1i*2*k0*zVec); % [nPix, nTauFit], v_xgp/p-independent, reused across every candidate

    for iVx = 1:nVx
        vxCand = v_xgp_vec(iVx);
        envelopeX = exp(-(vxCand.*tau).^2 ./ (4*sigma(1)^2)); % [1, nTauFit] -- p-independent, hoisted out of the inner loop
        for iP = 1:nP
            pCand = p_vec(iP);
            envelope = envelopeX .* exp(-(v_zgp.*tau).^2 ./ (4*sigma(2)^2) - (pCand.*zVec.*k0).^2); % [nPix, nTauFit]
            g1_model = DC + F.*envelope.*envelopeZphase; % [nPix, nTauFit]

            R2cand = 1 - sum(abs(g1 - g1_model).^2, 2) ./ ydataDev; % [nPix,1]

            improved = R2cand > bestR2;
            bestR2(improved) = R2cand(improved);
            v_xgp(improved) = vxCand;
            p(improved) = pCand;
        end
    end

    R2 = bestR2;
end
