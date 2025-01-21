% For use with a 1D array and a single 2D image plane

function [IQ_f, SVs, V_sort] = svd_declutter_1D_acrossframes(IQ_coherent_sum, sv_threshold_lower, sv_threshold_upper)
    
    %% SVD processing to extract motion, with the covariance method
    
    [zp, xp, nf] = size(IQ_coherent_sum);
    
    % Main data matrix PP to be manipulated
    
    % %%
    PP = zeros(xp*zp, nf); % x*y*z pixels by nf matrix
    % CM_V = zeros(nf, nf);
    % V = zeros(nf, nf);
    % V_sort = zeros(nf, nf);
    % D_V = zeros(nf, nf);
    % D_V_Sort = zeros(1, nf);
    % D_V_Sort_f = zeros(1, nf);
    % D_V_Sort_f_keep1 = zeros(1, nf);
    % D_V_Sort_f_keep2 = zeros(1, nf);
    % ind_D_V_Sort = zeros(1, nf);
    
    
    % for f = 1:nf  % go through all (sub)frames
        %%% go through each angle and stack all the data for each x value
            for x = 1:xp
                PP((x-1)*zp + 1:x*zp, :) = IQ_coherent_sum(:, x, :);
            end
    
        CM_V = PP'*PP; % covariance matrix for V
        % CM_V = P'*P;
        % sv = diag(CM); % singular values are the diagonal of the covariance matrix??
        
        % [U, D_U] = eig(CM_V); % get the angular eigenvectors that are also the right singular vector Vi of P
        % [D_U_Sort, ind_D_U_Sort]= sort(diag(D_U),'descend'); % get the eigenvalues in descending order
        
        [V_temp, D_V_temp] = eig(CM_V); % get the angular eigenvectors that are also the right singular vector Vi of P
        V = V_temp;        % eigenvectors
        D_V = D_V_temp;    % eigenvalues
    
        [D_V_Sort_temp, ind_D_V_Sort_temp]= sort(diag(D_V), 'descend'); % get the eigenvalues in descending order
        D_V_Sort = D_V_Sort_temp;
        ind_D_V_Sort = ind_D_V_Sort_temp;
    
        V_sort = V(:, ind_D_V_Sort); % rearrange V to be with the descending eigenvalues
        
        D_V_Sort_f = D_V_Sort;
        if length(D_V_Sort_temp) > sv_threshold_upper
            D_V_Sort_f([1:sv_threshold_lower - 1, sv_threshold_upper + 1:end]) = 0; % get rid of the data for eigenvalues past a threshold
        end
    
        % I_f = sqrt(diag(D_V_Sort_f)); % prev lab code did sqrt
        I_f = eye(size(diag(D_V_Sort_f)));
    
        % I_f(:, sv_threshold+1:end) = 0; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        I_f(:, [1:sv_threshold_lower - 1, sv_threshold_upper + 1:end]) = 0;  
    
        P_f = PP * V_sort * I_f * V_sort'; % filtered beamformed/reconstructed data
    
    % end
    
    % P_f_coherent_sum_keep1 = squeeze(sum(P_f_keep1, 2)); % coherently sum across angles
    % P_f_coherent_sum_keep2 = squeeze(sum(P_f_keep2, 2)); % coherently sum across angles
    % 
    % IQ_f = zeros(zp, xp, nf);
    % IQ_f_keep1 = IQ_f;
    % IQ_f_keep2 = IQ_f;
    % for f = 1:nf
        for x = 1:xp
            IQ_f(:, x, :) = P_f( (x-1)*zp + 1:x*zp, :);
        end
    % end
    
    SVs = D_V_Sort; % (# acquisitions, # frames)
    
    %% unnecessary plotting
    
    abs_IQ_f = abs(IQ_f);
    
    figure;
    imagesc(abs_IQ_f(:, :, 1))
    title(strcat("SVD filtered with ", num2str(sv_threshold_lower), " to ", num2str(sv_threshold_upper), " ordered singular values kept (Subframe 1)"))
    xlabel('x pixels')
    ylabel('z pixels')
    figure;
    imagesc(abs_IQ_f(:, :, end))
    title(strcat("SVD filtered with ", num2str(sv_threshold_lower), " to ", num2str(sv_threshold_upper), " ordered singular values kept (Last subframe)"))
    xlabel('x pixels')
    ylabel('z pixels')
end