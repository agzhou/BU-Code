%%
% Plot the frequency spectrum of singular vectors
% component

%%
% function [] = plot_FFT_SVs(V_sort, P)
%     V_sort_temp = abs(V_sort);
    V_sort_temp = V_sort; % don't do abs!!!!

    V_sort_zeroed = V_sort_temp - mean(V_sort_temp, 1); %%%
%     V_sort_zeroed = V_sort_temp;
    singVecFS = fft(V_sort_zeroed, [], 1);

    singVecFS_shifted = fftshift(singVecFS, 1);
%     singVecFS_shifted = singVecFS;
    subFrameRate = fps_target;
%     tstep = 1/subFrameRate;
    faxis = linspace(-subFrameRate/2, subFrameRate/2, length(singVecFS_shifted));
    
    test = abs(singVecFS_shifted);
    % 1D plot
%     figure; plot(faxis, test(:, 1))
    % 2D plot
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
    
%%
%     figure; plot(tstep:tstep:tstep*nf, temp); xlabel('time (s)'); ylabel('real(IQ coherent sum)')
    
%     figure; plot(faxis, abs(fft_val)); ylabel('magnitude of FFT'); xlabel('frequency (Hz)')
%     fD_predicted = -2*Trans.frequency*1e6*flow_v_mm_s/1e3/Resource.Parameters.speedOfSound
%     xline(fD_predicted)
%     xline(-fD_predicted)

% end