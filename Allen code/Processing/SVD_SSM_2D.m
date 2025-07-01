% Determine the optimal SV thresholds with the spatial similarity matrix
    [zp, xp, nf] = size(IQ);
    PP = reshape(IQ, [zp*xp, nf]);
    [U, S, V] = svd(PP); % Already sorted in decreasing order
    
    SSM = zeros(nf, nf); % Initialize the spatial similarity matrix

    SSM_const = 1/(zp * xp); % constant in front of the summation term
%%
    tic
    for n = 1:nf
%     for n = 1:10
        abs_u_n = abs(U(:, n)); % The nth column vector from U
        mean_abs_u_n = sum(abs_u_n) / length(abs_u_n);
        stddev_abs_u_n = std(abs_u_n);
%         for m = 1:nf
        for m = 1:n % leverage the symmetry of the SSM
            abs_u_m = abs(U(:, m)); % The mth column vector from U
            mean_abs_u_m = sum(abs_u_m) / length(abs_u_m);
            SSM(n, m) = sum( ((abs_u_n - mean_abs_u_n) .* (abs_u_m - mean_abs_u_m)) ...
                        ./ stddev_abs_u_n ...
                        ./ std(abs_u_m) );
        end
    end
    SSM = SSM .* SSM_const;
    toc
    figure; imagesc(SSM)