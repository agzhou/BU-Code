% For use with raw RF data from ultrasound measurements
% Input:  RF Data (unstacked, if the acquisition stacks frames)
% Output: PP     (space vs. frame/time data matrix),
%         EVs    (sorted eigenvalues, square of singular values)
%         V_sort (sorted eigenvectors)

function [PP, EVs, V_sort] = getSVsRF(RFData)
    
    %% SVD processing to extract motion, with the covariance method
    RFData = double(RFData);
    [ns, nc, nf] = size(RFData); % The RF data should have the dimensions [# samples, # channels, # frames]
    
    % Main data matrix PP to be manipulated
%     PP = zeros(ns*nc, nf); % # samples * # channels by nf matrix
    
%     % stack all the data for each x value
%     for x = 1:xp
%         PP((x-1)*zp + 1:x*zp, :) = RFData(:, x, :);
%     end
    % Stack all the data for each channel
    PP = reshape(RFData, [ns*nc, nf]);

    CM_V = PP'*PP; % covariance matrix for V

    [V, D_V] = eig(CM_V); % get the angular eigenvectors that are also the right singular vector Vi of P
    % V: eigenvectors
    % D_V: (unsorted) eigenvalues

    [D_V_Sort, ind_D_V_Sort]= sort(diag(D_V), 'descend'); % get the eigenvalues in descending order

    V_sort = V(:, ind_D_V_Sort); % rearrange V to be with the descending eigenvalues

    EVs = D_V_Sort; % (# acquisitions, # frames)
    
    % Look at the singular vectors' power density spectrums
%     plot_FFT_SVs_function(V_sort, P)

end