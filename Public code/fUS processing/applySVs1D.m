% Take the SVD processed parameters and apply the thresholding on some
% lower and upper SV limit

function [IQ_f, varargout] = applySVs1D(IQ_coherent_sum, PP, EVs, V_sort, sv_threshold_lower, sv_threshold_upper)

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

%     IQ_f = zeros(zp, xp, nf);
    % Unstack the spatial dimension
%     for x = 1:xp
%         IQ_f(:, x, :) = P_f( (x-1)*zp + 1:x*zp, :);
%     end
    IQ_f = reshape(P_f, [zp, xp, nf]); % Change 06/25/2025

    %% Jianbo's Noise thing, adapted (added 8/7/25)
    % I believe this is looking at the last 50 singular subspaces and
    % turning them into an image, then using a smoothing filter and taking
    % the average across frames
    
    UDelta = PP*V_sort;
    Vnoise = zeros(size(V_sort));
    Vnoise(:, end-50:end) = V_sort(:, end - 50:end);
    Noise = reshape(UDelta*Vnoise', [zp, xp, nf]);
    % sNoiseMed=medfilt2(abs(squeeze(mean(Noise,3))),[50 50],'symmetric');
    sNoiseMed = imgaussfilt(abs(squeeze(mean(Noise, 3))), 25);
    Noise = sNoiseMed/min(sNoiseMed(:));
    varargout{1} = Noise; % Assign the noise term to the optional output
end