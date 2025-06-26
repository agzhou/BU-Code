% Take the SVD processed parameters and apply the thresholding on some
% lower and upper SV limit

% Output:
%           RFData_cf (clutter filtered RF data)

function [RFData_cf] = applySVsRF(RFData, PP, EVs, V_sort, sv_threshold_lower, sv_threshold_upper)

    [ns, nc, nf] = size(RFData); % The RF data should have the dimensions [# samples, # channels, # frames]

    EVs_f = EVs;
    
    if length(EVs_f) < sv_threshold_upper
        error('Upper threshold is larger than the number of eigenvalues (frames)')
    elseif length(EVs_f) == sv_threshold_upper
        EVs_f(1:sv_threshold_lower - 1) = 0; % get rid of the data for eigenvalues past a threshold
    else
        EVs_f([1:sv_threshold_lower - 1, sv_threshold_upper + 1:end]) = 0; % get rid of the data for eigenvalues past a threshold
    end

    I_f = eye(size(diag(EVs_f)));

    I_f(:, [1:sv_threshold_lower - 1, sv_threshold_upper + 1:end]) = 0;

    P_f = PP * V_sort * I_f * V_sort'; % filtered beamformed/reconstructed data

    % Unstack the spatial dimension
%     RFData_cf = zeros(ns, nc, nf);
    RFData_cf = reshape(P_f, [ns, nc, nf]);
%     for x = 1:xp
%         RFData_cf(:, x, :) = P_f( (x-1)*zp + 1:x*zp, :);
%     end

end