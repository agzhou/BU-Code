% expects pixels as a row vector
function [psf_gFit, psf_gFit_values, pixels_finer] = PSFGaussianFit(pixels, psf, gfit_pixel_spacing)

    gfit_type = 'gauss2';
    [lb, ub] = findLocalMinsOfPSF(psf);
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