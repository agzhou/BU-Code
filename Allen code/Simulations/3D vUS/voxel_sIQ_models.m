% Description: calculate the sIQ values for numerical simulations

% Inputs:
%   voxel: struct with a 'data' field: a (nVoxels, 4) matrix. Columns 1-3 are the x, y, z coordinates;
%          column 4 are their reflectivities
%   sigma: field-based sigma vector for the PSF.
%          1. For the Single-Gaussian model: [sigma_x, sigma_y, sigma_z]
%          2. For the overlaid-Gaussians model: [sigma_narrow, sigma_wide, sigma_axial]
%   model: 'SG' or 'RCA' to specify the PSF model
%   SP: simulation parameter struct (see main_3D_vUS_simulation.m)


function sIQ = voxel_sIQ_models(voxel, sigma, model, SP)
    k = 2*pi/SP.wl; % Angular wavenumber

    switch model
        case 'SG' % Single-Gaussian PSF model
            PSF = exp( -(voxel.data(:, 1) - voxel.center(1)).^2./(2.*sigma(1)^2) ...
                  - (voxel.data(:, 2) - voxel.center(2)).^2./(2.*sigma(2)^2) ...
                  - (voxel.data(:, 3) - voxel.center(3)).^2./(2.*sigma(3)^2) );
    
        case 'RCA' % Two overlaid Gaussians PSF model (RCA OPW)
            sigma_n = sigma(1); sigma_w = sigma(2); sigma_a = sigma(3); % parse sigma vector
            PSF = ( exp( -(voxel.data(:, 1) - voxel.center(1)).^2./(2.*sigma_n^2) -(voxel.data(:, 2) - voxel.center(2)).^2./(2.*sigma_w^2)) ...
                  + exp( -(voxel.data(:, 1) - voxel.center(1)).^2./(2.*sigma_w^2) -(voxel.data(:, 2) - voxel.center(2)).^2./(2.*sigma_n^2)) ) ...
                  .* (voxel.data(:, 3) - voxel.center(3)).^2./(2.*sigma_a^2);
    end

    % Calculate the sIQ using the specific PSF term
    sIQ = sum(  voxel.data(:, 4) ...
                        .* PSF ...
                        .* exp( 2.*1i.*k.*(voxel.data(:, 3) - voxel.center(3)) )    );
    
    % Add white Gaussian noise
    sIQ = awgn(sIQ, SP.snr, 'measured', 'db');
end