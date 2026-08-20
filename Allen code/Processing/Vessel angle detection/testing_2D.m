%%
sigmas = 1:5;
spacing = [PData.PDelta(1) .* P.wl, PData.PDelta(3) .* P.wl]; % x, z pixel size [m]

% Create vectors defining x and z coordinates in actual length units
x_mm = (1:PData.Size(2)) .* PData.PDelta(1) .* P.wl .* 1e3; % x [mm]
z_mm = (1:PData.Size(1)) .* PData.PDelta(3) .* P.wl .* 1e3; % z [mm]
%%
vesselness = vesselness2D(PDI./noise, sigmas, spacing, tau, brightondark);