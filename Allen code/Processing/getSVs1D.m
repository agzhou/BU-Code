% For use with a 1D array and a single 2D image plane
% Input:  IQ coherent sum across angles, lower and upper singular value #
%         bounds
% Output: PP     (space vs. frame/time data matrix),
%         EVs    (sorted eigenvalues, square of singular values)
%         V_sort (sorted eigenvectors),

function [PP, EVs, V_sort] = getSVs1D(IQ_coherent_sum)
    
    %% SVD processing to extract motion, with the covariance method
    
    [zp, xp, nf] = size(IQ_coherent_sum);
    
    % Main data matrix PP to be manipulated
    PP = zeros(xp*zp, nf); % x*z pixels by nf matrix
    
    % stack all the data for each x value
    for x = 1:xp
        PP((x-1)*zp + 1:x*zp, :) = IQ_coherent_sum(:, x, :);
    end

    CM_V = PP'*PP; % covariance matrix for V

    [V, D_V] = eig(CM_V); % get the angular eigenvectors that are also the right singular vector Vi of P
    % V: eigenvectors
    % D_V: eigenvalues

    [D_V_Sort, ind_D_V_Sort]= sort(diag(D_V), 'descend'); % get the eigenvalues in descending order

    V_sort = V(:, ind_D_V_Sort); % rearrange V to be with the descending eigenvalues

    EVs = D_V_Sort; % (# acquisitions, # frames)
    
%     %% unnecessary plotting
%     
%     abs_IQ_f = abs(IQ_f);
%     
%     figure;
%     imagesc(abs_IQ_f(:, :, 1))
%     title(strcat("SVD filtered with ", num2str(sv_threshold_lower), " to ", num2str(sv_threshold_upper), " ordered singular values kept (Subframe 1)"))
%     xlabel('x pixels')
%     ylabel('z pixels')
%     figure;
%     imagesc(abs_IQ_f(:, :, end))
%     title(strcat("SVD filtered with ", num2str(sv_threshold_lower), " to ", num2str(sv_threshold_upper), " ordered singular values kept (Last subframe)"))
%     xlabel('x pixels')
%     ylabel('z pixels')
end