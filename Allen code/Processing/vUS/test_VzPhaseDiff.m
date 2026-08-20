    %% 
    x = [0.0, 0.05, 0.0, 1, 0];
    % [g1] = awgn(g1vUS2D_vec(x, tau, sigma, k0, true, false), 20).';
    % [g1] = awgn(g1vUS2D_vec(x, tau, sigma, k0, true, false), 200).';
    [g1] = g1vUS2D_vec(x, tau, sigma, k0, true, false).';
    figure; plot(g1)
    figure; plot(tau, abs(g1))

    %% Adjacent-lag phase differences and their reliability weights
    dphi = angle( g1(:, 2:end) .* conj(g1(:, 1:end-1)) );  % [nVoxels, nTau-1], each in (-pi, pi]
    w    = abs(g1(:, 2:end)) .* abs(g1(:, 1:end-1));       % down-weights already-decorrelated pairs

    % Magnitude-weighted circular mean of the phase-rotation rate
    wsum = sum(w, 2);
    wsum(wsum == 0) = eps; % guard fully decorrelated / no-signal voxels against 0/0
    rate = sum(w .* dphi, 2) ./ wsum / dt; % [nVoxels, 1], rad/s --> rate = (2*k0*Vz*dt)/dt

    Vz0 = rate ./ (2*k0); % Eq. 15: phase term is exp(i*2*k0*Vz0*tau)
