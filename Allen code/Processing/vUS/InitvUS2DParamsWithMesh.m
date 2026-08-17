
% Description: initialize vUS (2D) parameters with the mesh --> R^2 method as in the vUS paper

% Inputs:
% vectorized

function [v_zgp, v_xgp, p, DC, F, R2] = InitvUS2DParamsWithMesh(g1, v_zgp0, DC0, F0, PP, sigma, tau)
    nPix = PP.zp*PP.xp;

    % Smooth the initial guesses
    if PP.xp*PP.zp > 1
        v_zgp0 = reshape( medfilt2(reshape(v_zgp0, [PP.zp, PP.xp]), [5, 5]), [PP.zp * PP.xp, 1] ); % 2D median filter with kernel size (5, 5)
        F0 = reshape( medfilt2(reshape(F0, [PP.zp, PP.xp]), [5, 5]), [PP.zp * PP.xp, 1] ); % 2D median filter with kernel size (5, 5)
    end
    
    v_zgp = v_zgp0; % We trust the v_zgp guess a lot from the phase difference method, so don't change it
    
    % Create mesh for Vx, p
    v_xgp_bounds = [0, 30].*1e-3; % Bounds for v_xgp values [m/s]
    v_xgp_step = 1*1e-3; % Step for v_xgp grid [m/s]
    v_xgp_vec = v_xgp_bounds(1):v_xgp_step:v_xgp_bounds(2);
    p_vec = linspace(0, 1, 10); % p vector of values from 0-1
    
    [v_xgp_grid, p_grid] = ndgrid(v_xgp_vec, p_vec);
    % v_xgp_mesh = repmat(permute(v_xgp_grid, [3, 4, 1, 2]), [nPix, PP.nTau, 1, 1]); % Dimensions [nPix, nTau, # v_xgp grid points, # p grid points]
    % p_mesh = repmat(permute(p_grid, [3, 4, 1, 2]), [nPix, PP.nTau, 1, 1]); % Dimensions [nPix, nTau, # v_xgp grid points, # p grid points]
    v_xgp_mesh = repmat(permute(v_xgp_grid, [3, 1, 2]), [nPix, 1, 1]); % Dimensions [nPix, # v_xgp grid points, # p grid points]
    p_mesh = repmat(permute(p_grid, [3, 1, 2]), [nPix, 1, 1]); % Dimensions [nPix, # v_xgp grid points, # p grid points]

    % Reshape the other values to get the right shape
    s = size(p_mesh);
    v_zgp_mesh = repmat(v_zgp, [1, s(2), s(3)]);
    DC_mesh = repmat(DC0, [1, s(2), s(3)]);
    F_mesh = repmat(F0, [1, s(2), s(3)]);

    tau_mesh = permute(repmat(tau, [1, nPix]), [2, 1]);

    % Calculate g1 model
    g1_model = @(tau, sigma, k0, v_zgp, v_xgp, p, DC, F) DC + F.*exp(-(v_xgp.*tau).^2 ./ (4*sigma(1)^2) - (v_zgp.*tau).^2 ./ (4*sigma(2)^2) -(p.*v_zgp.*k0.*tau).^2) .* exp(1i .* 2 .* k0 .* v_zgp .* tau);
    k0 = 2*pi/PP.wl;
    g1_model_values = g1_model(tau, sigma, k0, v_zgp_mesh, v_xgp_mesh, p_mesh, DC_mesh, F_mesh);
end