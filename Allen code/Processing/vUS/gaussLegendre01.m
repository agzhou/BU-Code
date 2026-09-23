function [s, w] = gaussLegendre01(N)
%% Description:
%   Gauss-Legendre nodes and weights on [0, 1] (Golub-Welsch algorithm), for
%   the fixed-quadrature form of the corrected Tangelder-profile g1 model
%   (vUS_3D_quad_*.m).
%
% Inputs:
%   N: number of nodes (48 matched direct integration to <= 1.6e-8 for lags
%      up to 20 ms in the vUS_3D_quad_* tests)
%
% Outputs:
%   s: [1, N] nodes in (0, 1) (row, so tau(:) .* s expands to [nTau, N])
%   w: [N, 1] weights (column, so E*w integrates over s)

    n = 1:N-1;
    beta = n ./ sqrt(4*n.^2 - 1);
    T = diag(beta, 1) + diag(beta, -1);
    [V, D] = eig(T);
    [xk, idx] = sort(diag(D));
    w = 2 * V(1, idx).^2;      % weights on [-1, 1]
    s = ((xk + 1) / 2).';      % nodes mapped to [0, 1]
    w = w(:) / 2;              % weights mapped to [0, 1]
end
