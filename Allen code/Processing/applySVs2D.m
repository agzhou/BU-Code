% Take the SVD processed parameters and apply the thresholding on some
% lower and upper SV limit

function [IQ_f] = applySVs2D(IQ_coherent_sum, PP, EVs, V_sort, sv_threshold_lower, sv_threshold_upper)

    %% SVD processing test with the covariance method
    
    [xp, yp, zp, nf] = size(IQ_coherent_sum);

    EVs_f = EVs; % filtered eigenvalues

%     if length(EVs_f) > sv_threshold_upper
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
    
    P_f = PP * V_sort * I_f * V_sort'; % filtered beamformed/reconstructed data, x*y*z by # frames
    
    IQ_f = zeros(xp, yp, zp, nf);

    % Unstack the spatial dimension for the final filtered beamformed/reconstructed volumetric data
    for x = 1:xp
        for y = 1:yp
            IQ_f(x, y, :, :) = P_f( (x-1)*yp*zp + (y-1)*zp + 1 : (x-1)*yp*zp + y*zp, :);
        end
    end
    
end

