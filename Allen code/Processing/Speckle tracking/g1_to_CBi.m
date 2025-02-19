%% to do
% add directional filtering? I don't think it's applicable
%%

function [CBFi, CBVi] = g1_to_CBi(g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
% use calculated g1 to get rCBF
    CBF = @(tau_2, tau_1, g1_tau_2, g1_tau_1) 1./sqrt(tau_2 ^2 - tau_1 ^2) .* sqrt(abs(log10(abs(g1_tau_1) ./ abs(g1_tau_2))));
    CBV = @(tau_1, g1_tau_1) abs(squeeze(g1_tau_1)) ./ (1 - abs(squeeze(g1_tau_1)));


    if numel(size(g1)) == 3                                           % linear array
        CBFi = CBF(tau(tau2_index_CBF), tau(tau1_index_CBF), g1(:, :, tau2_index_CBF), g1(:, :, tau1_index_CBF));
        % CBV calculation relies on noise...
        CBVi = CBV(tau(tau1_index_CBV), g1(:, :, tau1_index_CBV));
    else
        CBFi = CBF(tau(tau2_index_CBF), tau(tau1_index_CBF), g1(:, :, :, tau2_index_CBF), g1(:, :, :, tau1_index_CBF)); % row column array
        % CBV calculation relies on noise...
        CBVi = CBV(tau(tau1_index_CBV), g1(:, :, :, tau1_index_CBV));

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