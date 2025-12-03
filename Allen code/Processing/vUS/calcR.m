% Equation 17 of the 2020 Jianbo Tang paper: coefficient of determination R
% Inputs:
%   g1exp: experimental g1(tau) (vector)
%   tau: vector of time delays tau
%   ..
%   sigma: 1x3 vector of the system PSF width [sigma_x, sigma_y, sigma_z]

function R = calcR(g1exp, tau, F, v_xgp, v_ygp, v_zgp, sigma, p0, k0)

    % Could probably vectorize these expressions
    % so that it isn't voxel by voxel...........
    % (make R a matrix and change the dimension of the mean operator)
    % But think about how v_xgp, p0, etc. may be mesh grids
    numer = mean( abs( g1exp - ( F.*exp(-(v_xgp .* tau).^2 ./ (4 * sigma(1)^2) - (v_ygp .* tau).^2 ./ (4 * sigma(2)^2) - (v_zgp .* tau).^2 ./ (4 * sigma(3)^2)) .* exp(-(p0 .* v_zgp .* k0 .* tau).^2) .* exp(2.*i.*k0.*tau.*v_zgp) ) ).^2 );
    denom = mean( abs(g1exp - mean(g1exp)) ) .^ 2; % SStotal
    R = 1 - numer./denom;


end
