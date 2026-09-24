function [r, J, cost] = vUS_3D_quad_batchEval(X, sel, P)
%% Description:
%   Batched residual + analytic Jacobian of the corrected Tangelder-profile
%   g1 model (fixed-quadrature form, see vUS_3D_quad_vec.m) for Bc voxels at
%   once. This is the "all voxels together" counterpart of
%   vUS_3D_quad_residJac.m / vUS_3D_quad_complex_Jac.m, and agrees with them
%   to ~1e-15.
%
%   Only ONE [L, N, Bc] complex exp array is formed (L lags, N quadrature
%   nodes, Bc voxels). Every Jacobian column is a weighted moment of it,
%   computed for all voxels with a single batched matrix product (pagemtimes):
%       G_m  = sum_j w_j * fw_m(j) * E(:, j)   for 7 weight sets fw_m:
%              1, f, f^2, df/dk, f*df/dk, df/da, f*df/da
%       g1   = DC + F*G_1
%       dg1/dv_t = F * (-tau^2 v_t/(2 sigma_t^2))                * G_{f^2}
%       dg1/dv_z = F * (-tau^2 v_z/(2 sigma_z^2) G_{f^2} + i 2 k0 tau G_f)
%       dg1/dF   = G_1,   dg1/dDC = 1
%       dg1/dk   = F * (-(M tau^2/2) G_{f df/dk} + i 2 k0 v_z tau G_{df/dk})
%       dg1/da   = F * (-(M tau^2/2) G_{f df/da} + i 2 k0 v_z tau G_{df/da})
%   with df/da, df/dk from the symbolic derivation in generate_vUS_3D_quad_Jac.m.
%
%   Works for any numeric class of X (double, single, gpuArray).
%
% Inputs:
%   X: [6, Bc] SCALED parameters (x./xscale), order [v_tgp, v_zgp, F, DC, k, a]
%   sel: [1, Bc] indices of the active voxels into the chunk data P.YR/P.YI/P.MASK
%   P: struct from vUS_3D_quad_batchPack.m with fields
%       sigma k0 s [1,N] logs [1,N] w [N,1] tauw [L,1] xscale [6,1] YR YI MASK [L,1,B]
%
% Outputs (already masked to each voxel's fit window):
%   r: [2L, 1, Bc] [real; imag] residual (model - data)
%   J: [2L, 6, Bc] Jacobian w.r.t. the SCALED parameters
%   cost: [1, Bc] sum of squared residuals

    Bc = size(X, 2);
    x  = X .* P.xscale;
    vt = reshape(x(1,:),1,1,Bc);  vz = reshape(x(2,:),1,1,Bc);
    F  = reshape(x(3,:),1,1,Bc);  DC = reshape(x(4,:),1,1,Bc);
    k  = reshape(x(5,:),1,1,Bc);  a  = max(reshape(x(6,:),1,1,Bc), realmin('like', X));   % a = 0 (plug flow) is finite here

    M    = vt.^2/P.sigma(1)^2 + vz.^2/P.sigma(3)^2;        % 1x1xBc
    ak   = a.^k;  D = k + 2 - 2*ak;  fmax = (k + 2)./D;
    sk   = P.s.^(k/2);                                       % 1xNxBc
    u    = ak.*sk;
    f    = fmax.*(1 - u);
    la   = log(a);
    fa   = k.*(k+2).*a.^(k-1).*(2 - (k+2).*sk)./D.^2;       % df/da
    fk   = 2*ak.*((k+2).*la - 1).*(1 - u)./D.^2 - fmax.*u.*(la + 0.5*P.logs);   % df/dk

    T  = P.tauw;  T2 = T.^2;                                 % Lx1
    E  = exp(complex(-(M/4).*T2.*f.^2, (2*P.k0*vz).*T.*f));  % LxNxBc  (the only big transcendental)

    fN = permute(f,[2 1 3]); fkN = permute(fk,[2 1 3]); faN = permute(fa,[2 1 3]);   % Nx1xBc
    W  = cat(2, P.w.*ones(1,1,Bc,'like',fN), P.w.*fN, P.w.*fN.^2, P.w.*fkN, P.w.*fkN.*fN, P.w.*faN, P.w.*faN.*fN);  % Nx7xBc
    G  = pagemtimes(E, W);                                   % Lx7xBc : all moments in one batched GEMM
    G0 = G(:,1,:); G1 = G(:,2,:); G2 = G(:,3,:); Gk0 = G(:,4,:); Gk1 = G(:,5,:); Ga0 = G(:,6,:); Ga1 = G(:,7,:);

    g   = DC + F.*G0;
    Jvt = F.*(-(T2.*vt)./(2*P.sigma(1)^2)).*G2;
    Jvz = F.*(-(T2.*vz)./(2*P.sigma(3)^2).*G2 + 1i*2*P.k0.*T.*G1);
    Jk  = F.*(-(M/2).*T2.*Gk1 + 1i*2*P.k0.*vz.*T.*Gk0);
    Ja  = F.*(-(M/2).*T2.*Ga1 + 1i*2*P.k0.*vz.*T.*Ga0);
    Jc  = cat(2, Jvt, Jvz, G0, ones(size(G0),'like',real(G0)), Jk, Ja);          % Lx6xBc complex

    msk = P.MASK(:,:,sel);  m2 = [msk; msk];
    r   = [real(g) - P.YR(:,:,sel); imag(g) - P.YI(:,:,sel)] .* m2;
    J   = [real(Jc); imag(Jc)] .* reshape(P.xscale,1,6) .* m2;
    cost = reshape(sum(r.^2, 1), 1, Bc);
end
