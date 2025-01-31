%% Changes
% 12/06/2024: uses an input of the coherent sum instead of IQ with the
%             angle pages
% 12/09/2024: removed abs from CM_V and added sqrt to D_V


%% For use with a 2D array and multiple 2D slices in a 3D image space
function [IQ_f, SVs, V_sort] = svd_declutter_2D_acrossframes_test_test(IQ_coherent_sum, sv_threshold_lower, sv_threshold_upper)
    % number of (ordered) singular values to use/keep
    tstart = clock;
    %% SVD processing test with the covariance method
    
    % IQ = IData{1, 1} + 1i.*QData{1, 1}; 
    % dimensions (# x pixels, # y pixels, # z pixels, # pages (# angles*2), # frames)
    
    [xp, yp, zp, nf] = size(IQ_coherent_sum);

    %%
    PP = zeros(xp*yp*zp, nf); % x*y*z pixels by na by nf matrix
    % CM_V = zeros(nf, nf);
    % V = zeros(np, np, nf);
    % V_sort = zeros(np, np, nf);
    % D_V = zeros(np, np, nf);
    % D_V_Sort = zeros(np, nf);
    % D_V_Sort_f = zeros(np, nf);
    % ind_D_V_Sort = zeros(np, nf);
    % I_f = zeros(np, np, nf);
    % P_f = zeros(xp*yp*zp, np, nf);
    IQ_f = zeros(xp, yp, zp, nf);

    
    %% test
    
    %%% for each x value, stack all the (z) data for each y value
            for x = 1:xp
                for y = 1:yp
        %             PP( (x-1)*yp + (y-1)*zp + 1 : (x-1)*yp + y*zp, p) = IQ(x, y, :, p);
                    PP( (x-1)*yp*zp + (y-1)*zp + 1 : (x-1)*yp*zp + y*zp, :) = IQ_coherent_sum(x, y, :, :);
                end
            end
    
%         figure; imagesc(abs(squeeze(PP(:, :, 10))))
    %%

        % if ~evalin('base', 'exist(''PP'', ''var'')')
            
        % end
        
        % CM_U = abs(PP*PP'); % covariance matrix for U

        CM_V = PP'*PP; % covariance matrix for V
        % CM_V = PP'*PP;
        % sv = diag(CM); % singular values are the diagonal of the covariance matrix??
        
        % [U, D_U] = eig(CM_V); % get the angular eigenvectors that are also the right singular vector Vi of PP
        % [D_U_Sort, ind_D_U_Sort]= sort(diag(D_U),'descend'); % get the eigenvalues in descending order
        
        [V_temp, D_V_temp] = eig(CM_V); % get the angular eigenvectors that are also the right singular vector Vi of PP
        V = V_temp;        % eigenvectors
        D_V = sqrt(D_V_temp);    % eigenvalues
    
        [D_V_Sort_temp, ind_D_V_Sort_temp]= sort(diag(D_V), 'descend'); % get the eigenvalues in descending order
        D_V_Sort = D_V_Sort_temp;
        ind_D_V_Sort = ind_D_V_Sort_temp;
    
        V_sort = V(:, ind_D_V_Sort); % rearrange V to be with the descending eigenvalues
        
        D_V_Sort_f = D_V_Sort;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if length(D_V_Sort_temp) > sv_threshold_upper
            D_V_Sort_f([1:sv_threshold_lower - 1, sv_threshold_upper + 1:end]) = 0; % get rid of the data for eigenvalues past a threshold
        end
    
%         D_V_Sort_f(:, f) = D_V_Sort_f_temp;
%         I_f = sqrt(diag(D_V_Sort_f)); % prev lab code did sqrt
        I_f = eye(size(diag(D_V_Sort_f)));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        I_f(:, [1:sv_threshold_lower - 1, sv_threshold_upper + 1:end]) = 0;  
        
        P_f = PP * V_sort * I_f * V_sort'; % filtered beamformed/reconstructed data
    
        % figure(2)
        % imagesc(abs(PP))
        % figure(3)
        % imagesc(abs(P_f))
        %
    
    % P_f_coherent_sum = sum(P_f, 2); % coherently sum across angles
    
    SVs = D_V_Sort; % (# acquisitions, # frames)
    
    % Unstack the spatial dimension
 
        for x = 1:xp
            for y = 1:yp
        %         IQ_f(x, y, :) = P_f_coherent_sum( (x-1)*yp + (y-1)*zp + 1 : (x-1)*yp + y*zp );
                IQ_f(x, y, :, :) = P_f( (x-1)*yp*zp + (y-1)*zp + 1 : (x-1)*yp*zp + y*zp, :);
            end
        end
    
    % for x = 1:xp
    %     img(:, x) = P_f_coherent_sum((x-1)*zp + 1:x*zp); %%%%%%%%%%%%%%%%%%
    % end
    
    % figure(3)
    %%
    % [X, Y, Z] = meshgrid(1:xp, 1:yp, 1:zp);
    % ptZ = unique(Media.MP(:, 3)); % all z points of scatterers from simulation
    % xslice = xp/2 + 1;   % middle
    % yslice = yp/2 + 1;   % middle
    % % zslice = ptZ*2;
    % zslice = ceil((ptZ ./ endDepth) .* zp);
    % 
    abs_IQ_f = abs(IQ_f);
    
    %% Plot processed planes
    fn = 1;
    figure
    imagesc(squeeze(abs_IQ_f(:, 40, :, fn))')
    title('SVD processed yz plane at x = 0 mm')
    ylabel('z pixels')
    xlabel('y pixels')
    
    figure
    imagesc(squeeze(abs_IQ_f(40, :, :, fn))')
    title('SVD processed xz plane at y = 0 mm')
    ylabel('z pixels')
    xlabel('x pixels')

    tend = clock;
    disp(strcat("SVD processing done, elapsed time is ", num2str(etime(tend, tstart)), "s"))
end

