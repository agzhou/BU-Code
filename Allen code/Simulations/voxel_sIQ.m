function sIQ = voxel_sIQ(voxel, SP)
    % k = 1/SP.wl;
    k = 2*pi/SP.wl;
    sIQ = sum(  voxel.data(:, 4) ...
                .* exp( -(voxel.data(:, 1) - voxel.center(1)).^2./(2.*SP.sigma(1)^2) -(voxel.data(:, 2) - voxel.center(2)).^2./(2.*SP.sigma(2)^2) -(voxel.data(:, 3) - voxel.center(3)).^2./(2.*SP.sigma(3)^2) ) ...
                .* exp( 2.*1i.*k.*(voxel.data(:, 3) - voxel.center(3)) )    );
    
    sIQ = awgn(sIQ, SP.snr, 'measured', 'db'); % Add white Gaussian noise
end