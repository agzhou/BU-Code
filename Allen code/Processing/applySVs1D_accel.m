% Take the SVD processed parameters and apply the thresholding on some
% lower and upper SV limit

function [IQ_f, varargout] = applySVs1D_accel(IQ_coherent_sum, PP, EVs, V_sort, sv_threshold_lower, sv_threshold_upper)

    [zp, xp, nf] = size(IQ_coherent_sum);

    EVs_f = EVs;
    
%     if length(EVs) > sv_threshold_upper
%         EVs_f([1:sv_threshold_lower - 1, sv_threshold_upper + 1:end]) = 0; % get rid of the data for eigenvalues past a threshold
%     else
%         error('Upper threshold is larger than the number of eigenvalues (frames)')
%     end

    if length(EVs_f) < sv_threshold_upper
        error('Upper threshold is larger than the number of eigenvalues (frames)')
    end

    V_sel = V_sort(:, sv_threshold_lower:sv_threshold_upper);

    P_f = (PP * V_sel) * V_sel';
%     end
    IQ_f = reshape(P_f, [zp, xp, nf]);

    %% Jianbo's Noise thing, adapted (added 8/7/25)
    % I believe this is looking at the last 50 singular subspaces and
    % turning them into an image, then using a smoothing filter and taking
    % the average across frames
    
    if nargout > 1
        V_noise = V_sort(:, end-50:end);
        Noise = reshape((PP * V_noise) * V_noise', [zp, xp, nf]);
        sNoiseMed = imgaussfilt3(abs(squeeze(mean(Noise, 3))), 25);
        Noise = sNoiseMed/min(sNoiseMed(:));
        varargout{1} = Noise; % Assign the noise term to the optional output
    end

end