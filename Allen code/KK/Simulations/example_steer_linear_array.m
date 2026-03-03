% Steer a linear array transducer in 2D
clearvars

%% create the computational grid
Nx = 128;           % number of grid points in the x (row) direction
Ny = 128;           % number of grid points in the y (column) direction
dx = 0.1e-3;        % grid point spacing in the x direction [m]
dy = 0.1e-3;        % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);
% kgrid.dt = 6e-6;
% kgrid.Nt = 1000;
%% define the properties of the propagation medium
medium.sound_speed = 1500;  % [m/s]
medium.alpha_coeff = 0.75;  % [dB/(MHz^power cm)]
medium.alpha_power = 1.5;

% medium.density = 1000 * ones(Nx, Ny);       % [kg/m^3]

kgrid.makeTime(medium.sound_speed);
%% define source mask for a linear transducer with an odd number of elements   
num_elements = 21;      % [grid points]
x_offset = 25;          % [grid points]
source.p_mask = zeros(Nx, Ny);
start_index = Ny/2 - round(num_elements/2) + 1;
source.p_mask(x_offset, start_index:start_index + num_elements - 1) = 1;

% define the properties of the tone burst used to drive the transducer
sampling_freq = 1/kgrid.dt;     % [Hz]
steering_angle = 30;            % [deg]
element_spacing = dx;           % [m]
tone_burst_freq = 1e6;          % [Hz]
tone_burst_cycles = 8;

% create an element index relative to the centre element of the transducer
element_index = -(num_elements - 1)/2:(num_elements - 1)/2;

% use geometric beam forming to calculate the tone burst offsets for each
% transducer element based on the element index
tone_burst_offset = 40 + element_spacing * element_index * ...
    sin(steering_angle * pi/180) / (medium.sound_speed * kgrid.dt);

% create the tone burst signals
source.p = toneBurst(sampling_freq, tone_burst_freq, tone_burst_cycles, ...
    'SignalOffset', tone_burst_offset);
%% % create a sensor mask covering the entire computational domain using the
% opposing corners of a rectangle
sensor.mask = [1, 1, Nx, Ny].';

% set the acoustic variables that are recorded
sensor.record = {'p'}; % p = pressure, u = particle velocity
%% Run the simulation
% create a display mask to display the transducer
display_mask = source.p_mask;
% assign the input options
input_args = {'DisplayMask', display_mask, 'PMLInside', false, 'PlotPML', false};

% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
