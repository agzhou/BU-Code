% Use calculated g1 to get CBF and CBV index with the time-lagged fUS
% analysis

%% to do
% add directional filtering? I don't think it's applicable
%%

function [CBFi, CBVi] = g1_to_CBi_NEW(g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV, n_CBV)
%     CBF = @(tau_2, tau_1, g1_tau_2, g1_tau_1) 1./(sqrt(tau_2 ^2 - tau_1 ^2) * 1) .* sqrt(abs(log(abs(g1_tau_1) ./ abs(g1_tau_2))));
%     CBV = @(tau_1, g1_tau_1) abs(squeeze(g1_tau_1)) ./ (1 - abs(squeeze(g1_tau_1)));
    CBF = @(tau, tau_1, g1_tau_2, g1_tau_1) 1./(sqrt(tau_2 ^2 - tau_1 ^2) * 1) .* sqrt(abs(log(abs(g1_tau_1) ./ abs(g1_tau_2))));
    CBV = @(tau_1, g1_tau_1) abs(squeeze(g1_tau_1)) ./ (1 - abs(squeeze(g1_tau_1)));


    if numel(size(g1)) == 3     % 2D data (linear array)
%         CBFi = CBF(tau(tau2_index_CBF), tau(tau1_index_CBF), g1(:, :, tau2_index_CBF), g1(:, :, tau1_index_CBF));
%         % CBV calculation relies on noise...
%         CBVi = CBV(tau(tau1_index_CBV), g1(:, :, tau1_index_CBV));
    else                        % 3D data
%         CBFi = sqrt(abs(log(abs(g1(:, :, :, tau1_index_CBF)) ./ abs(g1(:, :, :, tau2_index_CBF))))) ./ sqrt(tau(tau2_index_CBF)^2 - tau(tau1_index_CBF)^2);
%         CBFi = CBF(tau(tau2_index_CBF), tau(tau1_index_CBF), g1(:, :, :, tau2_index_CBF), g1(:, :, :, tau1_index_CBF)); % row column array
%         % CBV calculation relies on noise...
%         CBVi = CBV(tau(tau1_index_CBV), g1(:, :, :, tau1_index_CBV));
        CBFi = 1./(sqrt(tau(tau2_index_CBF) ^2 - tau(tau1_index_CBF) ^2) * 1) .* sqrt( abs( log( abs(g1(:, :, :, tau1_index_CBF)) ./ abs(g1(:, :, :, tau2_index_CBF)) ) ) );
        A = abs(g1(:, :, :, tau1_index_CBV)) ./ ( abs(g1(:, :, :, tau1_index_CBV .* n_CBV)) ./ abs(g1(:, :, :, tau1_index_CBV))) .^ ( 1/(n_CBV^2 - 1) );
%         A = abs(g1(:, :, :, tau1_index_CBV)) ./ ( abs(g1(:, :, :, tau1_index_CBV .* n_CBV)) ./ abs(g1(:, :, :, tau1_index_CBV))) .^ ( 1/(n_CBV - 1) ); % not squaring tau1 model
        CBVi = A ./ (1 - A);
    end

    % ztest = 36;
    % xtest = 206;
    % CBF_1(ztest, xtest)

%     figure
%     imagesc(CBF_1 ./ max(CBF_1, [], 'all'))
%     colormap hot
    
%     figure
%     imagesc(CBV_test ./ max(CBV_test, [], 'all'))
%     colormap hot

end