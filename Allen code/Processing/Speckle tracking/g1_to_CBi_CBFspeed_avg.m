% Use calculated g1 to get CBF and CBV index with the time-lagged fUS
% analysis

%% to do
% add directional filtering? I don't think it's applicable
%%

function [CBFi, CBVi] = g1_to_CBi_CBFspeed_avg(g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
    % CBF = @(tau_2, tau_1, g1_tau_2, g1_tau_1) 1./sqrt(tau_2 ^2 - tau_1 ^2) .* sqrt(abs(log10(abs(g1_tau_1) ./ abs(g1_tau_2))));
    CBF = @(tau_2, tau_1, g1_tau_2, g1_tau_1) 1./(sqrt(tau_2 ^2 - tau_1 ^2) * 1) .* sqrt(abs(log(abs(g1_tau_1) ./ abs(g1_tau_2))));
    % CBF = @(tau_2, tau_1, g1_tau_2, g1_tau_1) 1./sqrt(tau_2 ^2 - tau_1 ^2) .* sqrt(abs(abs(g1_tau_1) ./ abs(g1_tau_2)));
    CBV = @(tau_1, g1_tau_1) abs(squeeze(g1_tau_1)) ./ (1 - abs(squeeze(g1_tau_1)));


    if numel(size(g1)) == 3     % 2D data (linear array)
        counter = 0;
        for tau1 = tau1_index_CBF:tau2_index_CBF - 1 % Go through each possible pair of tau1 and tau2 in the prescribed range
            for tau2 = tau1 + 1:tau2_index_CBF
                if exist('CBFi', 'var')
                    CBFi = CBFi + CBF(tau(tau2), tau(tau1), g1(:, :, tau2), g1(:, :, tau1));
                else
                    CBFi = CBF(tau(tau2), tau(tau1), g1(:, :, tau2), g1(:, :, tau1));
                end
                counter = counter + 1;
            end
        end
        CBFi = CBFi ./ counter; % Average

        % CBV calculation relies on noise...
        CBVi = CBV(tau(tau1_index_CBV), g1(:, :, tau1_index_CBV));
    else                        % 3D data
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