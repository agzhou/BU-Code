% Take the SVD processed parameters and apply the thresholding on some
% lower and upper SV limit

function [IQ_f, varargout] = applySVs2D_accel(IQ_coherent_sum, PP, EVs, V_sort, sv_threshold_lower, sv_threshold_upper)

    %% SVD processing test with the covariance method
    
    [xp, yp, zp, nf] = size(IQ_coherent_sum);

    if sv_threshold_upper > length(EVs)
        error('Upper threshold is larger than the number of eigenvalues (frames)')
    end

    % Reconstructing with V_sort * I_f * V_sort' (I_f = identity with the
    % columns outside [lower, upper] zeroed) is algebraically identical to
    % V_sel * V_sel', where V_sel is just the kept columns of V_sort.
    % Slicing first means the two x*y*z-by-nf matrix products below only
    % cost ~x*y*z*nf*k flops (k = # of kept SVs) instead of ~x*y*z*nf^2,
    % and never forms the dense nf-by-nf projector at all. (Optimized 9/23/26)
    V_sel = V_sort(:, sv_threshold_lower:sv_threshold_upper);

    P_f = (PP * V_sel) * V_sel'; % filtered beamformed/reconstructed data, x*y*z by # frames

    IQ_f = reshape(P_f, [xp, yp, zp, nf]);


    %% Jianbo's Noise thing, adapted (added 8/7/25)
    % I believe this is looking at the last 50 singular subspaces and
    % turning them into an image, then using a smoothing filter and taking
    % the average across frames

    % Only compute this (another matrix product plus a 3D Gaussian filter
    % over the full volume) when the caller actually asks for it -- most
    % call sites only use IQ_f. (Optimized 9/23/26)
    if nargout > 1
        V_noise = V_sort(:, end-50:end);
        Noise = reshape((PP * V_noise) * V_noise', [xp, yp, zp, nf]);
        sNoiseMed = imgaussfilt3(abs(squeeze(mean(Noise, 4))), 25);
        Noise = sNoiseMed/min(sNoiseMed(:));
        varargout{1} = Noise; % Assign the noise term to the optional output
    end
end

