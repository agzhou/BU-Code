function P = vUS_3D_quad_batchPack(c, g1, tdi, t1i, tau, sigma, k0, s, w, xscale, cls)
%% Description:
%   Pack the voxels c into one padded batch for vUS_3D_quad_batchLM.m.
%   Data are laid out [L, 1, B] with L = the longest fit window in the batch;
%   shorter windows are zero-padded and marked in P.MASK. Sort voxels by window
%   length before chunking (as vUS_3D_quad_fitBatched.m does) so padding stays small.
%
% Inputs:
%   c: [B, 1] voxel indices (rows of g1) to pack
%   g1: [nVox, nTau] complex g1, tau starting at 0 (g1T convention)
%   tdi: [nVox, 1] index of the last lag to fit, per voxel (fit window = t1i:tdi)
%   t1i: index of the first lag to fit (2 = tau1)
%   tau: [nTau, 1] time lags [s]
%   sigma: [sigma_x, sigma_y, sigma_z] [m];  k0: wavenumber [rad/m]
%   s, w: Gauss-Legendre nodes [1,N] and weights [N,1] from gaussLegendre01.m
%   xscale: [1, 6] parameter scale (parameters are fitted as x./xscale)
%   cls: 'double' or 'single'
%
% Outputs:
%   P: struct with fields sigma k0 s logs w tauw xscale YR YI MASK (see vUS_3D_quad_batchEval.m)

    B = numel(c); win = tdi(c) - t1i + 1; L = max(win);
    Gw = g1(c, t1i:t1i+L-1);                                  % [B, L]
    msk = (1:L) <= win(:);                                    % [B, L]
    Gw(~msk) = 0;
    P.YR = cast(reshape(real(Gw).', L, 1, B), cls);  P.YI = cast(reshape(imag(Gw).', L, 1, B), cls);
    P.MASK = cast(reshape(msk.', L, 1, B), cls);
    P.tauw = cast(tau(t1i:t1i+L-1), cls);
    P.sigma = cast(sigma, cls); P.k0 = cast(k0, cls);
    P.s = cast(s(:).', cls); P.logs = log(P.s); P.w = cast(w(:), cls); P.xscale = cast(xscale(:), cls);
end
