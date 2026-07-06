%% Description: semi-automatically calculate the ultrasound PSF shape in 3D, taking in simulation data

%% Load simulation data and parameters (assuming Verasonics)
load('D:\Allen\Data\RC15gV PSF sim\params.mat')
load('D:\Allen\Data\RC15gV PSF sim\PSF.mat')

%% Get the position of the point scatterer - assuming there is only one scatterer
pos_wl = P.Media.MP(1, 1:3); % Position in wavelengths
pos_m = pos_wl .* P.wl; % [x, y, z] position in meters

%%