% Comparison of kWaveArray and kWaveTransducer Class

clearvars;
% close all

%% DEFINE LITERALS - Setting up parameters for the simulation

% Selection of K-Wave code execution model
model = 3;  % Options: 1 - MATLAB CPU, 2 - MATLAB GPU, 3 - C++ code, 4 - CUDA code
USE_STATISTICS = true;      % set to true to compute the rms or peak beam patterns, set to false to compute the harmonic beam patterns

% Medium parameters
c0 = 1540;        % Sound speed in the medium [m/s]
rho0 = 1020;      % Density of the medium [kg/m^3]

% Source parameters
source_f0 = (250/48)*1e6;         % Frequency of the ultrasound source [Hz]
source_amp = 1e6;          % Amplitude of the ultrasound source [Pa]
source_cycles = 3;         % Number of cycles in the tone burst signal
% source_focus = 5e-3;     % Focal length of the source [m]
numEl = 16;               % Number of elements in the transducer array
element_length = 2.3e-4;
element_width = 2.3e-4;    % Width of each transducer element [m]
element_pitch = 2.3e-4;    % Pitch - distance between the centers of adjacent elements [m]
elevation_length = 3e-3;   % Elevation Length - length along 3rd dimension of elements [m]
RF_fs = source_f0*4;       % Sampling Frequency of final RFData

% Define transmission angles for plane wave compounding
na = 1;  % Number of angles for transmission
if (na > 1)
    startAngle = -24*pi/180;
    thetaX = linspace(startAngle, -startAngle, na);
    thetaY = linspace(startAngle, -startAngle, na);
    [tX,tY] = meshgrid(thetaX,thetaY);
    TXangle = [tX(:),tY(:)];
else
    TXangle = [0*pi/180,0*pi/180];
end
na = size(TXangle,1);

% Transducer position parameters
translation = [0, 0];
rotation = 0;

% Grid parameters
grid_size_x = 5e-3;  % Grid size in x-direction [m]
grid_size_y = 5e-3;  % Grid size in y-direction [m]
grid_size_z = 5e-3;  % Grid size in z-direction [m]

% Computational parameters
ppw = 8;             % Points per wavelength
t_end = round(grid_size_z*2/c0,6);        % Simulation duration [s]
cfl = 0.5;            % Courant-Friedrichs-Lewy (CFL) number for stability

%% GRID - Creating the computational grid

% Calculate grid spacing based on PPW and source frequency
dx = c0 / (ppw * source_f0);  % Grid spacing [m]
dy = dx;
dz = dx;

% Compute grid size
Nx = roundEven(grid_size_x / dx);  % Number of grid points in x-direction
Ny = roundEven(grid_size_y / dy);  % Number of grid points in y-direction
Nz = roundEven(grid_size_z / dz);  % Number of grid points in z-direction

% Create the computational grid
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% Create the time array
kgrid.makeTime(c0, cfl, t_end);
dsFactor = (1/kgrid.dt)/RF_fs;

%% MEDIUM - Defining the medium properties 
medium.sound_speed = c0 * ones([Nx, Ny, Nz]);   % sound speed [m/s]
medium.density = rho0 * ones([Nx, Ny, Nz]);      % density [kg/m3]


%% SOURCE/SENSOR - KWaveArray

[karray, ElemPos] = initArray(kgrid, numEl, element_pitch, element_width, element_length);

% Plot Array
chkMask = karray.getArrayBinaryMask(kgrid);
[X,Y,Z] = meshgrid(kgrid.x_vec,kgrid.y_vec,kgrid.z_vec);
x = X(chkMask); y = Y(chkMask); z = Z(chkMask);
% Plot
figure
scatter3(x, y, z, 'SizeData', 1);
xlim([kgrid.x_vec(1) kgrid.x_vec(end)]);
ylim([kgrid.y_vec(1) kgrid.y_vec(end)]);
zlim([kgrid.z_vec(1) kgrid.z_vec(end)]);

arrayLen = element_length*numEl;
arrayWidth = element_width*numEl;
for i = 1:numEl
    line([ElemPos(i)+element_width/2, ElemPos(i)+element_width/2], [-arrayLen/2, arrayLen/2], [mean(z), mean(z)], 'Color', 'red', 'LineWidth', 2);    % Horizontal lines
    line([-arrayWidth/2, arrayWidth/2], [ElemPos(i)+element_length/2, ElemPos(i)+element_length/2], [mean(z), mean(z)], 'Color', 'green', 'LineWidth', 2);    % Vertical lines
end

xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('3D Scatter Plot of Logical Array');
grid on;
view(2)

% Create source signal using a tone burst
source_sig = source_amp .* toneBurst(1/kgrid.dt, source_f0, source_cycles);

% % Plotting the source signal
% figure;
% plot(kgrid.t_array(1:length(source_sig)) * 1e6, source_sig);
% xlabel('Microseconds (us)')
% title('Source Signal');

% Assign binary mask from karray to the sensor mask
sensor.mask = karray.getArrayBinaryMask(kgrid);

% set the record mode such that only the rms and peak values are stored
sensor.record = {'p'};

% Define frequency response of the sensor
sensor.frequency_response = [source_f0, 100];
%% SIMULATION - Running the simulation for different transmission angles

% Preallocate arrays for time delays and RF data
time_delays = zeros(numEl*numEl, na);

% Simulation input options
input_args = {'PMLSize', 'auto', 'PMLInside', false, 'PlotPML', false, 'DisplayMask', 'off','DeleteData',false};
RFData = zeros(numEl*numEl, kgrid.Nt, na);

% Loop over each angle for plane wave compounding
for i = 1:na
    % RFData based on kWaveArray
    [source, time_delays(:,i)] = genSource(kgrid, source_f0, source_cycles, source_amp, TXangle(i,:), karray, ElemPos, c0);
    sensor_data = runSim(kgrid, medium, source, sensor, input_args, model,source_amp);
    RFData(:, :, i) = karray.combineSensorData(kgrid, sensor_data.p);
end


% Rearrange RF data dimensions for further processing
RFData = downsample(flip(flip(reshape(permute(RFData, [2, 1, 3]),[kgrid.Nt,numEl,numEl,na]),2),3),dsFactor);

% figure; colormap gray
% imagesc(log10(abs(RFData)))
% xlabel('Horizontal Position [mm]')
% ylabel('Time [us]')
% title('RF data')

%% Beamforming  Parameter definition
% Define key parameter structure
param.fs = RF_fs;                           % [Hz]   sampling frequency
param.pitch = element_pitch;                % [m]
param.fc = source_f0;                       % [Hz]   center frequency
param.c = c0;                               % [m/s]  longitudinal sound speed
param.fnumber = 0.6;                        % [ul]   receive f-number

wavelength = param.c/param.fc;              % [m] convert from wavelength to meters
% samplesPerWave = param.fs/param.fc;     % the number of samples per wavelength
% note: this is off by a factor of two because you also account for
% roundtrip time. In otherwords, there are 4 samples per wavelength, but in
% practice that becomes 8 since you are also accounting for time to go to
% and from the transducer.

[~,I] = max(source_sig);
param.t0 = (kgrid.t_array(I))/param.fc; % Sequence start time (time offset)
param.TXdelay = time_delays;
param.DecimRate = 1;    % Decimation rate

xCoord = ((-numEl/2):0.25:(numEl/2)-1)*param.pitch;  % [m]   Beamformed points x coordinates
zCoord = (1:0.025:32)*wavelength;   % [m]    Beamformed points z coordinates
[X,Z] = meshgrid(xCoord,zCoord);

vsource = 10000*[tan(TXangle).',-ones(na,1)];  

%% Beamform
Recon = zeros(size(X,1),size(X,2),na);
for i = 1:na
    RFDataIQ = rf2iq(RFData(:,:,i),param);
    Recon(:,:,i) = ezdas(RFDataIQ,X,Z,vsource(i,:),param);
%     Recon(:,:,i) = das(RFData,X,Z,time_delays,param);
end

ReconC = abs(sum(Recon,3));
ReconC_log = 20*log10(ReconC/max(ReconC,[],'all'));

figure;
imagesc(xCoord*1e3,zCoord*1e3,ReconC_log);
axis image; colormap gray;colorbar
title('Beam formed image (dB)'); 
xlabel('Horizontal Position [mm]');
ylabel('Depth [mm]');

figure;
imagesc(xCoord*1e3,zCoord*1e3,ReconC);
axis image; colormap gray;
title('Beam formed image (linear)'); axis image
xlabel('Horizontal Position [mm]');
ylabel('Depth [mm]');

if na > 1
    genSliderV2(log10(abs(Recon)))
end


%% HELPER FUNCTIONS
function [karray, ElemPos] = initArray(kgrid, element_num, element_pitch, element_width, element_len)
    % Initializes the transducer array.
    % Args:
    %   kgrid: The k-Wave grid object.
    %   element_num: Number of elements in the array.
    %   element_pitch: Distance between the centers of adjacent elements.
    %   element_width: Width of each element.
    % Returns:
    %   karray: The k-Wave array object.
    %   ElemPos: The positions of the elements in the array.

    % Create empty kWaveArray object with specified BLI tolerance and upsampling rate
    karray = kWaveArray('BLITolerance', 0.05, 'UpsamplingRate', 10);

    % Calculate the center position for the first element
    L = element_num * element_pitch / 2;
    ElemPos = -(L - element_pitch / 2) + (0:element_num - 1) * element_pitch;
    [X,Y] = meshgrid(ElemPos,ElemPos);

    rotation = [0,0,0];
    % Add rectangular elements to the array
    for indY = 1:element_num
        for indX = 1:element_num
            % Set element position
            x_pos = X(indY,indX);
            y_pos = Y(indY,indX);

            % Define Rectangle dimensions
            position = [x_pos,y_pos,kgrid.z_vec(1)];
            Lx = element_width;
            Ly = element_len;

            % Add line element to the array
            karray.addRectElement(position, Lx, Ly, rotation);
        end
    end

%     karray.setArrayPosition([0,0,0], rotation)
end

function [source, time_delays] = genSource(kgrid, source_f0, source_cycles, source_amp, theta, karray, ElemPos, c0)
    % Generates the source signal with time delays for each transducer element.
    % Args:
    %   kgrid: The k-Wave grid object.
    %   source_f0: Frequency of the source.
    %   source_cycles: Number of cycles in the tone burst signal.
    %   source_amp: Amplitude of the source.
    %   theta: Steering angle of the plane wave.
    %   karray: The k-Wave array object.
    %   ElemPos: The positions of the elements in the array.
    %   c0: Speed of sound.
    % Returns:
    %   source: The source object containing the mask and signals.
    %   time_delays: Time delays applied to each element for focusing.

    % Calculate time delays for each element based on steering angle
    
    [X,Y] = meshgrid(ElemPos,ElemPos);
    
    time_delays0 = (X.*sin(theta(2))-Y.*sin(theta(1)))/c0;
    time_delays0 = time_delays0(:) - min(time_delays0(:));
    time_delays = time_delays0;
    
%     rng(10,'twister');
%     time_delays = 0.004*rand(length(ElemPos),length(ElemPos))/c0; %+ 

    
    % Create time-varying source signals for each physical element
    source_sig = source_amp .* toneBurst(1/kgrid.dt, source_f0, source_cycles, 'SignalOffset', round(time_delays / kgrid.dt));

    % Assign binary mask for the source
    source.p_mask = karray.getArrayBinaryMask(kgrid);

    % Generate per source signal from element signals. AKA Assign source signals to the source object
    source.p = karray.getDistributedSourceSignal(kgrid, source_sig);
end

function [sensor_data] = runSim(kgrid, medium, source, sensor, input_args, model, source_amp)
    % Runs the simulation based on the selected model (CPU or GPU).
    % Args:
    %   kgrid: The k-Wave grid object.
    %   medium: The medium in which waves propagate.
    %   source: The source object containing the ultrasound signal.
    %   sensor: The sensor object to record the pressure.
    %   input_args: Additional input arguments for the simulation.
    %   model: The selected model for running the simulation.
    % Returns:
    %   sensor_data: The recorded sensor data from the simulation.

    % Run the simulation based on the chosen model
    switch model
        case 1
            % MATLAB CPU
            sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, ...
                input_args{:}, ...
                'DataCast', 'single', ...
                'PlotScale', [-1, 1] * source_amp);

        case 2
            % MATLAB GPU
            sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, ...
                input_args{:}, ...
                'DataCast', 'gpuArray-single', ...
                'PlotScale', [-1, 1] * source_amp);

        case 3
            % C++ code
            sensor_data = kspaceFirstOrder3DC(kgrid, medium, source, sensor, input_args{:});

        case 4
            % C++/CUDA GPU
            sensor_data = kspaceFirstOrder3DG(kgrid, medium, source, sensor, input_args{:});
    end
end

