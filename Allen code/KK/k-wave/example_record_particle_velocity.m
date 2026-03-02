% Steer a linear array transducer in 2D
clearvars

%% create the computational grid
Nx = 128;           % number of grid points in the x (row) direction
Ny = 128;           % number of grid points in the y (column) direction
dx = 0.1e-3;        % grid point spacing in the x direction [m]
dy = 0.1e-3;        % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

%% define the properties of the propagation medium
medium.sound_speed = 1500;  % [m/s]
medium.alpha_coeff = 0.75;  % [dB/(MHz^power cm)]
medium.alpha_power = 1.5;

medium.density = 1000 * ones(Nx, Ny);       % [kg/m^3]
%% create initial pressure distribution in the shape of a disc, using makeDisc
disc_magnitude = 5; % [Pa]
disc_x_pos = 50;    % [grid points]
disc_y_pos = 50;    % [grid points]
disc_radius = 8;    % [grid points]
disc_1 = disc_magnitude * makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);

disc_magnitude = 3; % [Pa]
disc_x_pos = 80;    % [grid points]
disc_y_pos = 60;    % [grid points]
disc_radius = 5;    % [grid points]
disc_2 = disc_magnitude * makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);

source.p0 = disc_1 + disc_2;

%% define four sensor points centered about source.p0
sensor_radius = 40; % [grid points]
sensor.mask = zeros(Nx, Ny);
sensor.mask(Nx/2 + sensor_radius, Ny/2) = 1;
sensor.mask(Nx/2 - sensor_radius, Ny/2) = 1;
sensor.mask(Nx/2, Ny/2 + sensor_radius) = 1;
sensor.mask(Nx/2, Ny/2 - sensor_radius) = 1;

% set the acoustic variables that are recorded
sensor.record = {'p', 'u'}; % p = pressure, u = particle velocity
%% Run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor);
%% define source mask for a linear transducer with an odd number of elements   
num_elements = 21;      % [grid points]
x_offset = 25;          % [grid points]
source.p_mask = zeros(Nx, Ny);
start_index = Ny/2 - round(num_elements/2) + 1;
source.p_mask(x_offset, start_index:start_index + num_elements - 1) = 1;