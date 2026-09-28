function [g1, J] = vUS_3D_combined_complex_Jac(x, tau, k0)
%% Description:
%   [g1, J] for the combined-parameter model (vUS_3D_combined_vec.m / vUS_3D_combined_split.m), with the
%   analytic Jacobian w.r.t. all four parameters [C, v_zgp, F, DC]. Hand-derived directly (no symbolic
%   generator, unlike the RCA-PSF models' generate_vUS_3D_quad_RCA_xy_Jac.m): the model is a single
%   exponential of a quadratic-in-tau exponent, so all four derivatives are elementary:
%
%       g1 = DC + F*E,           E = exp(-C tau^2 + i 2 k0 v_zgp tau)
%       dg1/dC    = F*E*(-tau^2)
%       dg1/dv_zgp = F*E*(i 2 k0 tau)
%       dg1/dF    = E
%       dg1/dDC   = 1
%
% Inputs:
%   x: [C, v_zgp, F, DC];  tau: [nTau, 1] [s];  k0: [rad/m]  (see vUS_3D_combined_vec.m)
%
% Outputs:
%   g1: [nTau, 1] complex
%   J:  [nTau, 4] complex Jacobian of g1 w.r.t. [C, v_zgp, F, DC]

    C = x(1); v_zgp = x(2); F = x(3); DC = x(4);
    tau = tau(:);

    E = exp(-C.*tau.^2 + 1i*2*k0.*v_zgp.*tau);
    g1 = DC + F.*E;

    if nargout > 1
        J = [F.*E.*(-tau.^2), F.*E.*(1i*2*k0.*tau), E, ones(size(tau))];
    end
end
