%% SVD processing with the covariance method
% For use with a 2D array and a stack of volumetric data
% Input:  IQ coherent sum across angles
% Output: PP     (space vs. frame/time data matrix),
%         EVs    (sorted eigenvalues, square of singular values)
%         V_sort (sorted eigenvectors)

%% Changes
% 12/06/2024: uses an input of the coherent sum instead of IQ with the
%             angle pages
% 12/09/2024: removed abs from CM_V and added sqrt to D_V

%% Main function
function [PP, EVs, V_sort] = getSVs2D(IQ_coherent_sum)

    
    % IQ = IData{1, 1} + 1i.*QData{1, 1}; 
    % dimensions (# x pixels, # y pixels, # z pixels, # pages (# angles*2), # frames)
    
    [xp, yp, zp, nf] = size(IQ_coherent_sum);

    % Main data matrix PP to be manipulated
    PP = zeros(xp*yp*zp, nf); % x*y*z pixels by na by nf matrix

    % Reshape into a space x time (frame) matrix
    % For each x value, stack all the (z) data for each y value
    for x = 1:xp
        for y = 1:yp
            PP( (x-1)*yp*zp + (y-1)*zp + 1 : (x-1)*yp*zp + y*zp, :) = IQ_coherent_sum(x, y, :, :);
        end
    end

    CM_V = PP'*PP; % covariance matrix for V
    
    [V, D_V] = eig(CM_V); % get the angular eigenvectors that are also the right singular vector Vi of PP
    % V: eigenvectors
    % D_V: (unsorted) eigenvalues

    [D_V_Sort, ind_D_V_Sort]= sort(diag(D_V), 'descend'); % get the eigenvalues in descending order

    V_sort = V(:, ind_D_V_Sort); % rearrange V to be with the descending eigenvalues
    
    EVs = D_V_Sort; % sorted eigenvalues
    
end

