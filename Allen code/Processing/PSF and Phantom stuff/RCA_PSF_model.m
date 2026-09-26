

x0 = 0; y0 = 0; z0 = 0;
sigma_narrow = 50e-6;
sigma_wide = 150e-6;
sigma_axial = 40e-6;
x = linspace(x0 - sigma_wide*5, x0 + sigma_wide*5, 100);
y = linspace(y0 - sigma_wide*5, y0 + sigma_wide*5, 100);
z = linspace(z0 - sigma_wide*5, z0 + sigma_wide*5, 100);

[X, Y, Z] = meshgrid(x, y, z);


PSF = @(x, y, z, x0, y0, z0, sigma_narrow, sigma_wide, sigma_axial) ...
        ( exp(-(x - x0).^2 ./ (2*sigma_narrow^2) - (y - y0).^2 ./ (2*sigma_wide^2)) ...
        + exp(-(x - x0).^2 ./ (2*sigma_wide^2) - (y - y0).^2 ./ (2*sigma_narrow^2)) ) ...
        .* exp(-(z - z0).^2 ./ (2*sigma_axial^2));

test = PSF(X, Y, Z, x0, y0, z0, sigma_narrow, sigma_wide, sigma_axial);
volumeViewer(test)
