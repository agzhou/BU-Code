function [Y, J] = vUS_3D_quad_batchModel(X, tau, k0, sigma, s, w)
%% Description:
%   Corrected Tangelder-profile g1 model (fixed-quadrature form, see
%   vUS_3D_quad_vec.m) for Bc voxels at once, in the plug-in form that
%   fitBatchedLM.m expects:
%       model = @(X, cols) vUS_3D_quad_batchModel(X, tau(cols), k0, sigma, s, w);
%   The model part of vUS_3D_quad_batchEval.m (same formulas), with the data,
%   mask and parameter scaling left to fitBatchedLM.m.
%
%   Only ONE [L, N, Bc] complex exp array is formed (L lags, N quadrature
%   nodes, Bc voxels). With one output only the zeroth moment is formed. With
%   two outputs every Jacobian column is a weighted moment of the same array,
%   computed for all voxels with a single batched matrix product (pagemtimes),
%   see vUS_3D_quad_batchEval.m for the formulas.
%
%   Works for any numeric class of X (double, single, gpuArray).
%
% Inputs:
%   X: [6, Bc] parameters in PHYSICAL units, order [v_tgp, v_zgp, F, DC, k, a] (see vUS_3D_quad_vec.m)
%   tau: [L, 1] time lags [s];  k0: wavenumber [rad/m]
%   sigma: [sigma_x, sigma_y, sigma_z] [m]; sigma_x must equal sigma_y
%   s, w: Gauss-Legendre nodes [1, N] and weights [N, 1] on [0, 1], from gaussLegendre01.m
%
% Outputs:
%   Y: [L, Bc] complex g1 model
%   J: [L, 6, Bc] complex Jacobian dY/dX in physical units (only computed if requested)

    if sigma(1) ~= sigma(2)
        error('vUS_3D_quad_batchModel:sigmaMismatch', ...
            'For the transverse-combined model, sigma_x and sigma_y must be equal.')
    end
    Bc = size(X, 2);  L = numel(tau);
    tau = tau(:);  s = s(:).';  w = w(:);
    vt = reshape(X(1,:),1,1,Bc);  vz = reshape(X(2,:),1,1,Bc);
    F  = reshape(X(3,:),1,1,Bc);  DC = reshape(X(4,:),1,1,Bc);
    k  = reshape(X(5,:),1,1,Bc);  a  = max(reshape(X(6,:),1,1,Bc), realmin('like', X));   % a = 0 (plug flow) is finite here

    M    = vt.^2/sigma(1)^2 + vz.^2/sigma(3)^2;             % 1x1xBc
    ak   = a.^k;  D = k + 2 - 2*ak;  fmax = (k + 2)./D;
    sk   = s.^(k/2);                                        % 1xNxBc
    u    = ak.*sk;
    f    = fmax.*(1 - u);

    T  = tau;  T2 = T.^2;                                   % Lx1
    E  = exp(complex(-(M/4).*T2.*f.^2, (2*k0*vz).*T.*f));   % LxNxBc  (the only big transcendental)

    if nargout < 2
        Y = reshape(DC + F.*pagemtimes(E, w), L, Bc);
        return
    end

    la = log(a);
    fa = k.*(k+2).*a.^(k-1).*(2 - (k+2).*sk)./D.^2;         % df/da
    fk = 2*ak.*((k+2).*la - 1).*(1 - u)./D.^2 - fmax.*u.*(la + 0.5*log(s));   % df/dk

    fN = permute(f,[2 1 3]); fkN = permute(fk,[2 1 3]); faN = permute(fa,[2 1 3]);   % Nx1xBc
    W  = cat(2, w.*ones(1,1,Bc,'like',fN), w.*fN, w.*fN.^2, w.*fkN, w.*fkN.*fN, w.*faN, w.*faN.*fN);   % Nx7xBc
    G  = pagemtimes(E, W);                                  % Lx7xBc : all moments in one batched GEMM
    G0 = G(:,1,:); G1 = G(:,2,:); G2 = G(:,3,:); Gk0 = G(:,4,:); Gk1 = G(:,5,:); Ga0 = G(:,6,:); Ga1 = G(:,7,:);

    Yp  = DC + F.*G0;
    Jvt = F.*(-(T2.*vt)./(2*sigma(1)^2)).*G2;
    Jvz = F.*(-(T2.*vz)./(2*sigma(3)^2).*G2 + 1i*2*k0.*T.*G1);
    Jk  = F.*(-(M/2).*T2.*Gk1 + 1i*2*k0.*vz.*T.*Gk0);
    Ja  = F.*(-(M/2).*T2.*Ga1 + 1i*2*k0.*vz.*T.*Ga0);
    J   = cat(2, Jvt, Jvz, G0, ones(size(G0),'like',real(G0)), Jk, Ja);   % Lx6xBc complex
    Y   = reshape(Yp, L, Bc);
end
