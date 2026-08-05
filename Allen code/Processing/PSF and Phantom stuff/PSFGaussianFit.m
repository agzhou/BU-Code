% expects pixels as a row vector
% Optional input for the fit type
function [psf_gFit, psf_gFit_values, pixels_finer] = PSFGaussianFit(pixels, psf, gfit_pixel_spacing, varargin)

    % gfit_type = 'gauss2';
    gfit_type = 'gauss1'; % Default fit type: single-term Gaussian
    if nargin > 2
        gfit_type = varargin{1};
    end

%     [lb, ub] = findLocalMinsOfPSF(psf);
    lb = 1; ub = length(psf); % Change 7/6/26
    psf_cut = psf(lb:ub);
    pixels_cut = pixels(lb:ub);

    if size(psf_cut, 1) == 1 % row vector
        psf_gFit = fit(pixels_cut', psf_cut', gfit_type); % get the fit object
    else
        psf_gFit = fit(pixels_cut', psf_cut, gfit_type); % get the fit object
    end

    pixels_finer = lb:gfit_pixel_spacing:ub; % define finer point spacing for more accurate FWHM
    psf_gFit_values = psf_gFit(pixels_finer);

end