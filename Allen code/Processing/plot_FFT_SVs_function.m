%%
% Plot the frequency spectrum of singular vectors
% component

%%
function [] = plot_FFT_SVs_function(V_sort, P)
    V_sort_temp = V_sort;
    V_sort_zeroed = V_sort_temp - mean(V_sort_temp, 1); %%%

    singVecFS = fft(V_sort_zeroed, [], 1);

    singVecFS_shifted = fftshift(singVecFS, 1);

%     subFrameRate = P.fps_target;
    frameRate = P.frameRate;
%     tstep = 1/subFrameRate; 
    faxis = linspace(-frameRate/2, frameRate/2, length(singVecFS_shifted))';
    
    % test = abs(singVecFS_shifted);
    % 1D plot
%     figure; plot(faxis, test(:, 1))
    % 2D
%     figure; imagesc(20.*log10(abs(singVecFS_shifted)))
    figure; imagesc(abs(singVecFS_shifted))
    colormap hot
    xlabel('Singular vector number')
    ylabel('Frequency (Hz)')
    title('Temporal singular vector power density spectrum')
%     temp_axis = axis;
%     axis([temp_axis(1:2), [min(faxis), max(faxis)]]);
    yticks([0; str2double(yticklabels)])
    yticklabels_temp = yticklabels;
    newyticklabels = linspace(min(faxis), max(faxis), numel(yticks));
    yticklabels(num2cell(newyticklabels));
    

%     figure; plot(tstep:tstep:tstep*nf, temp); xlabel('time (s)'); ylabel('real(IQ coherent sum)')
    
%     figure; plot(faxis, abs(fft_val)); ylabel('magnitude of FFT'); xlabel('frequency (Hz)')
%     fD_predicted = -2*Trans.frequency*1e6*flow_v_mm_s/1e3/Resource.Parameters.speedOfSound
%     xline(fD_predicted)
%     xline(-fD_predicted)

    % Look at the "mean frequency" of each subspace of the temporal singular vectors
    % Do a weighted sum..
    
    singVec_weighted_mean_freq = -sum(abs(singVecFS_shifted) .* repmat(abs(faxis), 1, size(singVecFS_shifted, 2))) ./ length(faxis);
    % Negative so the curve decreases

    figure; plot(singVec_weighted_mean_freq)
    % figure; plot(movmean(singVec_weighted_mean_freq, 3))
    title('Weighted mean frequency of the temporal singular vectors'' power spectrum')
    xlabel('Singular vector number')
    ylabel('Frequency [Hz]')
end