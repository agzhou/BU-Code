function [Y, J] = vUS_3D_quad_RCA_xy_batchModel(X, tau, k0, sigma, s, w)
%% Description:
%   Two-step RCA-PSF g1 model (fixed-quadrature form, see vUS_3D_quad_RCA_xy_vec.m) for Bc voxels at
%   once, in the plug-in form fitBatchedLM.m expects:
%       model = @(X, cols) vUS_3D_quad_RCA_xy_batchModel(X, tau(cols), k0, sigma, s, w);
%       vi   = find(maskToUse);  nv = numel(vi);
%       x0   = [5e-3*ones(nv,1), 5e-3*ones(nv,1), Vz0(vi), ones(nv,1), zeros(nv,1), 2.5*ones(nv,1), ones(nv,1)];
%       lb   = [-50e-3, -50e-3, -50e-3, 0, 0, 2, 0];   ub = [50e-3, 50e-3, 50e-3, 1, 1, 3, 1];
%       opts = struct('window', tau_decayed_ind(vi), 't1i', 2, 'jacobian', 'analytic', ...
%                      'xscale', [1e-2 1e-2 1e-2 1 1 1 1]);
%       [x, cost] = fitBatchedLM(model, x0, g1_exp{j}(vi, :), lb, ub, opts);
%   Direct RCA-PSF analog of vUS_3D_quad_batchModel.m: same batching strategy (one big [L, N, Bc] complex
%   exp array per lateral term, one shared moment-weighting matrix, Jacobian columns as weighted moments
%   via pagemtimes), used here for THREE lateral terms (see below) instead of one. There is no
%   RCA-specific batched solver: fitBatchedLM.m / batchLM.m are model-agnostic and work with this model
%   directly (see their own headers). Agrees with vUS_3D_quad_RCA_xy_residJac.m to round-off
%   (test_vUS_3D_quad_RCA_xy_batch.m, check 1) and reduces exactly to vUS_3D_quad_batchModel.m when
%   sigma_narrow == sigma_wide (check 4).
%
%   With S = sigma_narrow^2 + sigma_wide^2, c = 2 sigma_narrow sigma_wide/S, w1 = 1/(2(1+c)), w23 = c/(1+c):
%       Integrand = EZ .* (w1*D1 + w1*D4 + w23*D23)
%       EZ = exp(i 2 k0 v_zgp tau f - (Mz/4) tau^2 f^2),           Mz  = v_zgp^2/sigma_axial^2
%       D1 = exp(-(Mt1/4) tau^2 f^2),  Mt1  = v_xgp^2/sigma_narrow^2 + v_ygp^2/sigma_wide^2
%       D4 = exp(-(Mt4/4) tau^2 f^2),  Mt4  = v_xgp^2/sigma_wide^2  + v_ygp^2/sigma_narrow^2
%       D23 = exp(-(Mt23/4) tau^2 f^2), Mt23 = 2(v_xgp^2+v_ygp^2)/S
%       g1  = DC + F * int_0^1 Integrand ds
%   Only THREE [L, N, Bc] complex arrays are formed, P_j = EZ.*D_j for j in {1, 4, 23} (L lags, N
%   quadrature nodes, Bc voxels) -- three times the single-Gaussian model's one, since the RCA model is a
%   weighted sum of three lateral Gaussian terms instead of one; everything past that point (the 7-column
%   moment-weighting matrix W, the pagemtimes moments) is identical in structure to
%   vUS_3D_quad_batchModel.m, done three times (once per P_j) and then linearly recombined with the
%   FIXED (parameter-independent) weights w1, w1, w23. With one output only the zeroth moment is formed.
%
%   Works for any numeric class of X (double, single, gpuArray).
%
% Inputs:
%   X: [7, Bc] parameters in PHYSICAL units, order [v_xgp, v_ygp, v_zgp, F, DC, k, a] (see vUS_3D_quad_RCA_xy_vec.m)
%   tau: [L, 1] time lags [s];  k0: wavenumber [rad/m]
%   sigma: [sigma_narrow, sigma_wide, sigma_axial] [m]
%   s, w: Gauss-Legendre nodes [1, N] and weights [N, 1] on [0, 1], from gaussLegendre01.m
%
% Outputs:
%   Y: [L, Bc] complex g1 model
%   J: [L, 7, Bc] complex Jacobian dY/dX in physical units, order matching X (only computed if requested)

    Bc = size(X, 2);  L = numel(tau);
    tau = tau(:);  s = s(:).';  w = w(:);
    vx = reshape(X(1,:),1,1,Bc);  vy = reshape(X(2,:),1,1,Bc);  vz = reshape(X(3,:),1,1,Bc);
    F  = reshape(X(4,:),1,1,Bc);  DC = reshape(X(5,:),1,1,Bc);
    k  = reshape(X(6,:),1,1,Bc);  a  = max(reshape(X(7,:),1,1,Bc), realmin('like', X));   % a = 0 (plug flow) is finite here

    sn = sigma(1); sw = sigma(2); sa = sigma(3);
    S  = sn^2 + sw^2;
    c  = 2*sn*sw/S;
    w1 = 1/(2*(1 + c));  w23 = c/(1 + c);            % scalars: w1 = w4, and 2*w1 + w23 = 1

    Mt1  = vx.^2/sn^2 + vy.^2/sw^2;                  % 1x1xBc
    Mt4  = vx.^2/sw^2 + vy.^2/sn^2;
    Mt23 = 2*(vx.^2 + vy.^2)/S;
    Mz   = vz.^2/sa^2;

    ak   = a.^k;  D_ = k + 2 - 2*ak;  fmax = (k + 2)./D_;
    sk   = s.^(k/2);                                 % 1xNxBc
    u    = ak.*sk;
    f    = fmax.*(1 - u);

    T  = tau;  T2 = T.^2;                            % Lx1
    EZ = exp(complex(-(Mz/4).*T2.*f.^2, (2*k0*vz).*T.*f));          % LxNxBc (phase + axial: shared by all three terms)
    D1 = exp(-(Mt1/4).*T2.*f.^2); D4 = exp(-(Mt4/4).*T2.*f.^2); D23 = exp(-(Mt23/4).*T2.*f.^2);   % LxNxBc, real
    P1 = EZ.*D1; P4 = EZ.*D4; P23 = EZ.*D23;         % LxNxBc, complex -- the "only big transcendental" arrays, x3

    if nargout < 2
        G0 = w1.*pagemtimes(P1, w) + w1.*pagemtimes(P4, w) + w23.*pagemtimes(P23, w);   % Lx1xBc
        Y = reshape(DC + F.*G0, L, Bc);
        return
    end

    la = log(a);
    fa = k.*(k+2).*a.^(k-1).*(2 - (k+2).*sk)./D_.^2;                                 % df/da
    fk = 2*ak.*((k+2).*la - 1).*(1 - u)./D_.^2 - fmax.*u.*(la + 0.5*log(s));         % df/dk

    fN = permute(f,[2 1 3]); fkN = permute(fk,[2 1 3]); faN = permute(fa,[2 1 3]);   % Nx1xBc
    W  = cat(2, w.*ones(1,1,Bc,'like',fN), w.*fN, w.*fN.^2, w.*fkN, w.*fkN.*fN, w.*faN, w.*faN.*fN);   % Nx7xBc, SHARED by all three P_j
    G1  = pagemtimes(P1,  W);   % Lx7xBc: (:,1)=int P1, (:,2)=int P1 f, (:,3)=int P1 f^2, (:,4)=int P1 fk,
    G4  = pagemtimes(P4,  W);   %         (:,5)=int P1 fk f, (:,6)=int P1 fa, (:,7)=int P1 fa f -- likewise for P4, P23
    G23 = pagemtimes(P23, W);

    % Combined (F- and vz-facing) moments: weighted sum over the three lateral terms
    G0  = w1.*G1(:,1,:) + w1.*G4(:,1,:) + w23.*G23(:,1,:);     % int Integrand       (= dY/dF)
    G1m = w1.*G1(:,2,:) + w1.*G4(:,2,:) + w23.*G23(:,2,:);     % int Integrand f
    G2m = w1.*G1(:,3,:) + w1.*G4(:,3,:) + w23.*G23(:,3,:);     % int Integrand f^2

    Yp  = DC + F.*G0;

    % v_xgp, v_ygp: each lateral term has its own sigma normalization (D1, D4 swap sigma_narrow/wide roles;
    % D23 is symmetric in v_xgp, v_ygp), so unlike v_zgp these need the PER-TERM second moment, not G2m.
    term_vx = w1.*G1(:,3,:)./sn^2 + w1.*G4(:,3,:)./sw^2 + w23.*G23(:,3,:).*(2/S);
    term_vy = w1.*G1(:,3,:)./sw^2 + w1.*G4(:,3,:)./sn^2 + w23.*G23(:,3,:).*(2/S);
    Jvx = F.*(-(T2.*vx)./2).*term_vx;
    Jvy = F.*(-(T2.*vy)./2).*term_vy;

    Jvz = F.*( -(T2.*vz)./(2*sa^2).*G2m + 1i*2*k0.*T.*G1m );

    % k, a: chain rule through f(s) needs, per lateral term, the fk/fa moments (columns 4-7 of G1/G4/G23)
    % both as a plain weighted sum (Gk0sum, Gk1sum) and weighted by that term's own Mt_j (Gk1weighted) --
    % see generate_vUS_3D_quad_RCA_xy_Jac.m's header / RCA_g1s_derivation.md for the derivation. This
    % reduces exactly to vUS_3D_quad_batchModel.m's "Jk = F.*(-(M/2).*T2.*Gk1 + 1i*2*k0.*vz.*T.*Gk0)" when
    % sigma_narrow == sigma_wide (Mt1 = Mt4 = Mt23, so Gk1weighted = Mt1.*Gk1sum and Mz+Mt1 = M).
    Gk0sum = w1.*G1(:,4,:) + w1.*G4(:,4,:) + w23.*G23(:,4,:);
    Gk1sum = w1.*G1(:,5,:) + w1.*G4(:,5,:) + w23.*G23(:,5,:);
    Gk1w   = w1.*Mt1.*G1(:,5,:) + w1.*Mt4.*G4(:,5,:) + w23.*Mt23.*G23(:,5,:);
    Ga0sum = w1.*G1(:,6,:) + w1.*G4(:,6,:) + w23.*G23(:,6,:);
    Ga1sum = w1.*G1(:,7,:) + w1.*G4(:,7,:) + w23.*G23(:,7,:);
    Ga1w   = w1.*Mt1.*G1(:,7,:) + w1.*Mt4.*G4(:,7,:) + w23.*Mt23.*G23(:,7,:);

    Jk = F.*( -(Mz/2).*T2.*Gk1sum + 1i*2*k0.*vz.*T.*Gk0sum - (T2/2).*Gk1w );
    Ja = F.*( -(Mz/2).*T2.*Ga1sum + 1i*2*k0.*vz.*T.*Ga0sum - (T2/2).*Ga1w );

    J = cat(2, Jvx, Jvy, Jvz, G0, ones(size(G0),'like',real(G0)), Jk, Ja);   % Lx7xBc complex
    Y = reshape(Yp, L, Bc);
end
