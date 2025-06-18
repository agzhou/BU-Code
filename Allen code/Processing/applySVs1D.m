% Take the SVD processed parameters and apply the thresholding on some
% lower and upper SV limit

function [IQ_f] = applySVs1D(IQ_coherent_sum, PP, EVs, V_sort, sv_threshold_lower, sv_threshold_upper)

    [zp, xp, nf] = size(IQ_coherent_sum);

    EVs_f = EVs;
    
%     if length(EVs) > sv_threshold_upper
%         EVs_f([1:sv_threshold_lower - 1, sv_threshold_upper + 1:end]) = 0; % get rid of the data for eigenvalues past a threshold
%     else
%         error('Upper threshold is larger than the number of eigenvalues (frames)')
%     end

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

    IQ_f = zeros(zp, xp, nf);
    % Unstack the spatial dimension
    for x = 1:xp
        IQ_f(:, x, :) = P_f( (x-1)*zp + 1:x*zp, :);
    end

end