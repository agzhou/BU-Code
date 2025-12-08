function sIQ = voxel_sIQ(voxel, SP)

    sIQ = sum(  voxel.data(:, 4) .* exp( -(voxel.data(1) - voxel.center(1)).^2./(2.*SP.sigma(1)^2) -(voxel.data(2) - voxel.center(2)).^2./(2.*SP.sigma(2)^2) -(voxel.data(1) - voxel.center(3)).^2./(2.*SP.sigma(3)^2) ) .* exp( 2.*1i./SP.wl.*(voxel.data(:, 3) - voxel.center(3)) )    );

end