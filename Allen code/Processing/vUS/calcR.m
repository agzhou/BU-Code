% Description: calculate the coefficient of determination R for all voxels at once. (Equation 17 of the 2020 Jianbo Tang vUS paper) 

% Inputs:
%   g1exp: experimental g1(tau) matrix with shape: [# voxels, # time lags]
%   tau: vector of time delays tau [1, # time lags]
%   F: vector with shape [# voxels, 1]
%   v_xgp: vector of x velocities per voxel, with shape [# voxels, 1]
%   v_ygp: vector of y velocities per voxel, with shape [# voxels, 1]
%   v_zgp: vector of z velocities per voxel, with shape [# voxels, 1]
%   sigma: 1x3 vector of the system PSF width [sigma_x, sigma_y, sigma_z]
%   p: vector with shape [# voxels, 1]
%   k0: angular wavenumber [rad/m] (a constant)

% Outputs:
%   R: [# voxels, 1] vector of R values

function R = calcR(g1exp, tau, F, v_xgp, v_ygp, v_zgp, sigma, p, k0)

    % Could probably vectorize these expressions
    % so that it isn't voxel by voxel...........
    % (make R a matrix and change the dimension of the mean operator)
    % But think about how v_xgp, p0, etc. may be mesh grids
    % numer = mean( abs( g1exp - ( F.*exp(-(v_xgp .* tau).^2 ./ (4 * sigma(1)^2) - (v_ygp .* tau).^2 ./ (4 * sigma(2)^2) - (v_zgp .* tau).^2 ./ (4 * sigma(3)^2)) .* exp(-(p .* v_zgp .* k0 .* tau).^2) .* exp(2.*1i.*k0.*tau.*v_zgp) ) ).^2 );
    % denom = mean( abs(g1exp - mean(g1exp)) ) .^ 2; % SStotal
    % R = 1 - numer./denom;

    num_voxels = size(g1exp, 1); % # of voxels to consider
    % Could add some input checking (e.g., if the # of voxels for v_xgp = num_voxels) here...

    % Expand the time lag vector so we can do multiplications across voxels
    if length(size(tau)) ~= 2 & size(tau, 1) ~= 1
        error('Time lag (tau) vector should be of shape [1, # time lags]')
    end
    tau_mat = repmat(tau, num_voxels, 1);

    numer = mean( abs( g1exp - ( F.*exp(-(v_xgp .* tau_mat).^2 ./ (4 * sigma(1)^2) - (v_ygp .* tau_mat).^2 ./ (4 * sigma(2)^2) - (v_zgp .* tau_mat).^2 ./ (4 * sigma(3)^2)) .* exp(-(p .* v_zgp .* k0 .* tau_mat).^2) .* exp(2.*1i.*k0.*tau_mat.*v_zgp) ) ).^2 );
    denom = mean( abs(g1exp - mean(g1exp)) ) .^ 2; % SStotal
    R = 1 - numer./denom;


end
