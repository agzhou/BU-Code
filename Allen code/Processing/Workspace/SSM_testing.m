%%
    tic
%     for n = 1:nf
    for n = 1:10
        abs_u_n = abs(U(:, n)); % The nth column vector from U
        mean_abs_u_n = sum(abs_u_n) / length(abs_u_n);
        stddev_abs_u_n = std(abs_u_n);
        for m = 1:nf
            abs_u_m = abs(U(:, m)); % The mth column vector from U
            mean_abs_u_m = sum(abs_u_m) / length(abs_u_m);
            SSM(n, m) = sum( ((abs_u_n - mean_abs_u_n) .* (abs_u_m - mean_abs_u_m)) ...
                        ./ stddev_abs_u_n ...
                        ./ std(abs_u_m) );
        end
    end
    SSM = SSM .* SSM_const;
    toc

    %% Test "std" speed
    clearvars std_test std_func_test
    u_m_length = length(abs_u_m);
    fac = u_m_length - 1;
    tic
%     std_test = sqrt ( sum( abs((abs_u_m - mean_abs_u_m)) .^ 2) ./ (u_m_length - 1) );
%     std_test = sqrt ( sum( (abs_u_m - mean_abs_u_m) .^ 2) ./ (u_m_length - 1) );
    std_test = sqrt ( sum( (abs_u_m - mean_abs_u_m) .^ 2) ./ fac );
%     std_test = ( sum( (abs_u_m - mean_abs_u_m) .^ 2) ./ fac ) .^ 0.5;
    toc

    tic
    std_func_test = std(abs_u_m);
    toc