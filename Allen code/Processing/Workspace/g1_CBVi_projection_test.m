% Use g1(tau2) and g1(tau1) to extrapolate back to a "g1(tau = 0)"
tau1_ind = 2;
tau2_ind = 3;
% g1_tau1_tau2_slope = (g1(:, :, :, tau2_ind) - g1(:, :, :, tau1_ind)) ./ taustep;
% g1_tau0 = g1(:, :, :, tau1_ind) + g1_tau1_tau2_slope .* (-taustep);
g1_tau1_tau2_slope = (abs(g1(:, :, :, tau2_ind)) - abs(g1(:, :, :, tau1_ind))) ./ taustep;

g1_tau1_tau2_slope_mask = g1_tau1_tau2_slope < 0; % Keep only voxels where g1(tau2) < g1(tau1)
% g1_tau0 = abs(g1(:, :, :, tau1_ind)) + g1_tau1_tau2_slope .* (-taustep);
g1_tau0 = abs(g1(:, :, :, tau1_ind));
g1_tau0(g1_tau1_tau2_slope_mask) = g1_tau0(g1_tau1_tau2_slope_mask) + g1_tau1_tau2_slope(g1_tau1_tau2_slope_mask) .* (-taustep);

CBVi_proj_test = abs(squeeze(g1_tau0)) ./ (1 - abs(squeeze(g1_tau0)));


figure; imagesc(squeeze(max(CBVi_proj_test, [], 1))')