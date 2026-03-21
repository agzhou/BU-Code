

% Adapted from Nikunj Khetan's code
% Comparison of kWaveArray and kWaveTransducer Class

clearvars;
% close all

%% Add k-wave and MUST toolboxes to path
addKWavePath
addMUSTPath

%% DEFINE LITERALS - Setting up parameters for the simulation

% Selection of K-Wave code execution model
model = 3;  % Options: 1 - MATLAB CPU, 2 - MATLAB GPU, 3 - C++ code, 4 - CUDA code
USE_STATISTICS = true;      % set to true to compute the rms or peak beam patterns, set to false to compute the harmonic beam patterns
SUBTRACT_BASELINE = true; % Set this flag = true to run a simulation with a homogeneous medium, to see only the crosstalk effects from the probe elements/sources not having any directivity

% Medium parameters
c0 = 1540;        % Sound speed in the medium [m/s]
rho0 = 1020;      % Density of the medium [kg/m^3]

% Source parameters
source_f0 = (250/48)*1e6;  % Frequency of the ultrasound source [Hz]
% source_f0 = (15)*1e6;  % Frequency of the ultrasound source [Hz]
source_amp = 1e6;          % Amplitude of the ultrasound source [Pa]
source_cycles = 3;         % Number of cycles in the tone burst signal
% source_focus = 5e-3;     % Focal length of the source [m]
element.num = 16;                % Number of elements in the transducer array (in one dimension for a square matrix array)
element.length = 2.3e-4;   % Length of each transducer element [m]
element.width = 2.3e-4;    % Width of each transducer element [m]
element.pitch = 2.3e-4;    % Pitch - distance between the centers of adjacent elements [m]
% element.length = 2.3e-4 / 2;   % Length of each transducer element [m]
% element.width = 2.3e-4 / 2;    % Width of each transducer element [m]
% element.pitch = 2.3e-4 / 2;    % Pitch - distance between the centers of adjacent elements [m]
element.elevationlength = 3e-3;   % Elevation Length - length along 3rd dimension of elements [m]
RF_fs = source_f0*4;       % Sampling Frequency of final RFData

% Define transmission angles for plane wave compounding
naTX = 5;  % Number of angles for transmission (in one dimension)
maxAngle = 5; % [deg]
if (naTX > 1)
    startAngle = -maxAngle*pi/180;
    thetaX = linspace(startAngle, -startAngle, naTX);
    thetaY = linspace(startAngle, -startAngle, naTX);
    [tX, tY] = meshgrid(thetaX, thetaY);
    anglesTX = [tX(:), tY(:)];
    daTX = mean(diff(thetaX)); % TX angle ingrement [rad]
else
    anglesTX = [0*pi/180, 0*pi/180];
    daTX = 0;
end

ntaTX = size(anglesTX, 1);
% nta = length(TXangle(:)); % # of total transmit angles

% Transducer position parameters
Trans.translation = [0, 0, 0]; % [m]
% Trans.rotation = 0;
Trans.rotation = [0, 0, 0]; % [Degrees]

% Grid parameters
% note: grid size = the span of the grid, not the spacing
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
dsFactor = (1/kgrid.dt)/RF_fs; % Downsampling factor to turn the RF data on the simulation timestep to the RF_fs's corresponding timestep

%% MEDIUM - Defining the medium properties 
medium.sound_speed = c0 * ones([Nx, Ny, Nz]);   % sound speed [m/s]
medium.density = rho0 * ones([Nx, Ny, Nz]);      % density [kg/m3]

if SUBTRACT_BASELINE
    % Fully homogeneous medium for the baseline subtraction
    medium_baseline.sound_speed = c0 * ones([Nx, Ny, Nz]);   % sound speed [m/s]
    medium_baseline.density = rho0 * ones([Nx, Ny, Nz]);      % density [kg/m3]
end

% Add a ball target to the "real" medium
bc_mm = [0, 0, 2]; % Ball center coordinates in mm (x, y, z)
bc = bc_mm ./ 1e3 ./ [dx, dy, dz]; % Ball center coordinates in grid points
br_um = 50; % Ball radius in um
% br_um = 400; % Ball radius in um
br = br_um ./ 1e6 ./ dx; % Assumes dx = dy = dz
ball_mask = logical(makeBall(Nx, Ny, Nz, bc(1), bc(2), bc(3), br));

% Modify the properties at the locations of the ball target
% medium.sound_speed(ball_mask) = c0 * 1;
% medium.density(ball_mask) = rho0 * 1.05;      % density [kg/m3]
% medium.sound_speed(ball_mask) = c0 * 2;
% medium.density(ball_mask) = rho0 * 1;      % density [kg/m3]
% medium.sound_speed(ball_mask) = 2600;
medium.sound_speed(ball_mask) = c0;
medium.density(ball_mask) = 1120;      % density [kg/m3]

% Plot the ball
figure; imagesc(kgrid.y_vec * 1e3, kgrid.z_vec*1e3, squeeze(max(ball_mask, [], 1))'); xlabel('y [mm]'); ylabel('z [mm]'); colorbar; axis image

%% SOURCE/SENSOR - KWaveArray

[karray, ElemPos, element.coords] = initArray(kgrid, element, Trans);

% % Plot Array
% chkMask = karray.getArrayBinaryMask(kgrid);
% [X,Y,Z] = meshgrid(kgrid.x_vec, kgrid.y_vec, kgrid.z_vec);
% x = X(chkMask); y = Y(chkMask); z = Z(chkMask);
% % Plot
% figure
% scatter3(x, y, z, 'SizeData', 1);
% xlim([kgrid.x_vec(1) kgrid.x_vec(end)]);
% ylim([kgrid.y_vec(1) kgrid.y_vec(end)]);
% zlim([kgrid.z_vec(1) kgrid.z_vec(end)]);
% 
% arrayLen = element.length*element.num;
% arrayWidth = element.width*element.num;
% for i = 1:element.num
%     line([ElemPos(i)+element.width/2, ElemPos(i)+element.width/2], [-arrayLen/2, arrayLen/2], [mean(z), mean(z)], 'Color', 'red', 'LineWidth', 2);    % Horizontal lines
%     line([-arrayWidth/2, arrayWidth/2], [ElemPos(i)+element.length/2, ElemPos(i)+element.length/2], [mean(z), mean(z)], 'Color', 'green', 'LineWidth', 2);    % Vertical lines
% end
% 
% xlabel('X-axis');
% ylabel('Y-axis');
% zlabel('Z-axis');
% title('3D Scatter Plot of Logical Array');
% grid on;
% axis image
% view(2)

% Create source signal using a tone burst
source_sig = source_amp .* toneBurst(1/kgrid.dt, source_f0, source_cycles);

% % % Plotting the source signal
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
time_delays = zeros(element.num*element.num, ntaTX);

% Simulation input options
% input_args = {'PMLSize', 'auto', 'PMLInside', false, 'PlotPML', false, 'DisplayMask', 'off','DeleteData',false};
input_args = {'PMLSize', 'auto', 'PMLInside', false, 'PlotPML', false, 'DisplayMask', 'off'};

RFData = zeros(element.num*element.num, kgrid.Nt, ntaTX);
RFDataBL = zeros(element.num*element.num, kgrid.Nt, ntaTX);
RFDataTarget = zeros(element.num*element.num, kgrid.Nt, ntaTX);

tic
% Loop over each angle for plane wave compounding
for ai = 1:ntaTX
    % RFData based on kWaveArray
    disp("******** "+ num2str(ai) + " ********")
    % **** FIX THE BELOW TXangle(ai, :)!!!!!!! ****
    [source, time_delays(:, ai)] = genSource(kgrid, source_f0, source_cycles, source_amp, anglesTX(ai, :), karray, ElemPos, c0);
    sensor_data = runSim(kgrid, medium, source, sensor, input_args, model, source_amp); % run simulation
    p = sensor_data.p;
    RFDataTarget(:, :, ai) = karray.combineSensorData(kgrid, p);
    % Additional simulation for baseline crosstalk subtraction
    if SUBTRACT_BASELINE
        sensor_data_baseline = runSim(kgrid, medium_baseline, source, sensor, input_args, model, source_amp); % run "baseline" simulation: crosstalk only
        p = p - sensor_data_baseline.p; % Baseline subtraction
        RFDataBL(:, :, ai) = karray.combineSensorData(kgrid, sensor_data_baseline.p);
    end

    RFData(:, :, ai) = karray.combineSensorData(kgrid, p); % Data from each array element stored with dimensions [total # elements, kgrid.Nt]
end
toc

% Rearrange RF data dimensions for further processing
% RFData = downsample(flip(flip(reshape(permute(RFData, [2, 1, 3]), [kgrid.Nt, element.num, element.num, nta]), 2), 3), dsFactor);
RFData_raw = RFData; % Non-downsampled RFData
RFData = downsample(permute(RFData_raw, [2, 1, 3]), dsFactor);
% RFData = downsample(flip(flip(reshape(permute(RFData_raw, [2, 1, 3]), [kgrid.Nt, element.num, element.num, nta]), 2), 3), dsFactor);

% figure; colormap gray
% imagesc(log10(abs(RFData)))
% xlabel('Horizontal Position [mm]')
% ylabel('Time [us]')
% title('RF data')

%% KK parameters

naRX = 25; % # of RX angles in 1 dimension

o = fix(-naRX/2):1:fix(naRX/2); % Truncate towards zero
j = fix(naRX/2); % Shift parameter
anglesRXList = (sign(o) .* daTX .* (2.*abs(o)./naRX + j))'; % Receive angles [deg]
anglesRX = listToAngles(anglesRXList); % All the receive angles [theta_x, theta_y]
ntaRX = size(anglesRX, 1); % Total number of RX angles

%% KK compression
% s = 2 * element.pitch * RF_fs / c0; % aspect ratio...
ratio = RF_fs / c0;
RawDataKK = DataCompressKKMatrixArray(RFData, anglesRX, ratio, element.coords);
nSamples = size(RawDataKK,1);
% figure; imagesc(squeeze(RawDataKK(:, :, round(ntaRX/2))))
% figure; imagesc(squeeze(RawDataKK(:, round(ntaTX/2), :)))

figure; imagesc(reshape(RawDataKK,[nSamples,ntaTX*ntaRX]))

%% Beamforming  Parameter definition
% Define key parameter structure
param.fs = RF_fs;                           % [Hz]   sampling frequency
param.pitch = element.pitch;                % [m]
param.fc = source_f0;                       % [Hz]   center frequency
param.c = c0;                               % [m/s]  longitudinal sound speed
param.fnumber = [0.6, 0.6];                        % [ul]   receive f-number
param.elements = element.coords; % Element coordinates (x, y) [m]
wavelength = param.c/param.fc;              % [m] convert from wavelength to meters
% samplesPerWave = param.fs/param.fc;     % the number of samples per wavelength
% note: this is off by a factor of two because you also account for
% roundtrip time. In otherwords, there are 4 samples per wavelength, but in
% practice that becomes 8 since you are also accounting for time to go to
% and from the transducer.

[~,I] = max(source_sig); % Find the index where the source input signal is maximum (wrt kgrid.dt)
param.t0 = (kgrid.t_array(I))/param.fc; % Sequence start time (time offset) [s] --> convert the index wrt kgrid.dt to time wrt the RF sampling frequency
% param.TXdelay = time_delays;
param.DecimRate = 1;    % Decimation rate

xCoord = ((-element.num/2):0.25:(element.num/2))*param.pitch;  % [m]   Beamformed points x coordinates
% xCoord = ((-element.num/2):0.5:(element.num/2))*param.pitch;  % [m]   Beamformed points x coordinates
yCoord = xCoord;
zbounds_mm = [0, 5]; % Z bounds/extents [mm]
zbounds = zbounds_mm ./ 1e3; % Z bounds/extents [m]
zCoord = zbounds(1):0.25*wavelength:zbounds(2);   % [m]    Beamformed points z coordinates
% zCoord = (1:0.025:32)*wavelength;   % [m]    Beamformed points z coordinates
[X, Y, Z] = meshgrid(xCoord, yCoord, zCoord);

BFgrid = struct('X', X, 'Y', Y, 'Z', Z); % Struct for the beamforming grid
% vsource = 10000*[tan(TXangle).',-ones(na,1)];  

%% KK beamforming

tic;
RawDataKKHilb = hilbert(RawDataKK);
% [BFData, LUTTX, LUTRX] = BeamformKK_MatrixArray(RawDataKKHilb, anglesRX, anglesTX, BFgrid, param);
[BFData] = BeamformKK_MatrixArray(RawDataKKHilb, anglesRX, anglesTX, BFgrid, param);
toc

%% Look at LUTs
% [testBFData, testLUTTX, testLUTRX] = BeamformKK_MatrixArray(RawDataKKHilb, [0, 0], [0, 0], BFgrid, param);

%% Examine KK beamforming result
% figure; imagesc(xCoord * 1e3, zCoord*1e3, squeeze(max(abs(BFData), [], 1))'); xlabel('y [mm]'); ylabel('z [mm]'); colorbar; axis image
figure; imagesc(xCoord * 1e3, zCoord*1e3, squeeze(max(abs(sum(BFData,[4, 5])), [], 1))'); xlabel('y [mm]'); ylabel('z [mm]'); colorbar; axis image

% volumeViewer(abs(BFData).^ 0.5)
% volumeViewer(abs(sum(BFData,[4,5])).^ 0.5)

% volumeViewer(abs(BFData(:,:,:,13,1)).^ 0.5)

%% DAS beamforming
Recon = zeros(size(X, 1), size(X, 2), size(X, 3), ntaTX); % Initialize container for storing reconstructed data

tic
for ai = 1:ntaTX % Go through every angle
% for ai = 1
    disp("Reconstructing volume with angle # " + num2str(ai))
    RFDataIQ = rf2iq(RFData(:, :, ai), param);
    % Recon(:,:,ai) = ezdas(RFDataIQ,X,Z,vsource(i,:),param);
    % Recon(:, :, :, ai) = das3(squeeze(RFData(:, :, ai)), X, Y, Z, time_delays, param);
    Recon(:, :, :, ai) = das3(RFDataIQ, X, Y, Z, time_delays(:, ai), param);

end
toc

%% Testing the DAS result
vol_CPWC_DAS = squeeze(sum(Recon, 4));
figure; imagesc(xCoord * 1e3, zCoord*1e3, squeeze(max(abs(vol_CPWC_DAS), [], 1))'); xlabel('y [mm]'); ylabel('z [mm]'); colorbar; axis image

% volumeViewer(abs(vol_CPWC_DAS).^ 0.5)

%%
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

if naTX > 1
    genSliderV2(log10(abs(Recon)))
end

%% HELPER FUNCTIONS
function [karray, ElemPos, elementCoords] = initArray(kgrid, element, Trans)
    % Initializes the transducer array.
    % Args:
    %   kgrid: The k-Wave grid object.
    %   element: a struct for transducer element parameters, with fields:
    %       num: Number of elements in the array (in one dimension).
    %       pitch: Distance between the centers of adjacent elements.
    %       width: Width of each element.
    %       length: length of each element
    % Returns:
    %   karray: The k-Wave array object.
    %   ElemPos: The positions of the elements in the array.

    % Create empty kWaveArray object with specified BLI tolerance and upsampling rate
    %   - BLITolerance: Scalar value controlling where the spatial extent of the band-limited interpolant (BLI) at each point is trunctated as a portion of the maximum value.
    %   - UpsamplingRate: Oversampling used to distribute the off-grid points compared to the equivalent number of on-grid points.
    karray = kWaveArray('BLITolerance', 0.05, 'UpsamplingRate', 10);
    
    % Calculate the center position for the first element
    L = element.num * element.pitch / 2; % Half-length of the full array (at least, between the start and end element centers along one dimension)
    ElemPos = -(L - element.pitch / 2) + (0:element.num - 1) * element.pitch; % Element center positions in one dimension
    [X, Y] = meshgrid(ElemPos, ElemPos); % Meshgrid of array elements
    elementCoords = [X(:), Y(:)]'; % Save the coordinates of each element as x, y pairs
    
    % rotation = [0, 0, 0];
    rotation = Trans.rotation;

    % Add rectangular elements to the array
    for indY = 1:element.num
        for indX = 1:element.num
            % Extract/set element position from the meshgrid
            x_pos = X(indY, indX);
            y_pos = Y(indY, indX);

            % Define Rectangle dimensions
            position = [x_pos, y_pos, kgrid.z_vec(1)]; % Position of the center of the element. kgrid.z_vec(1) refers to the topmost/starting grid coordinate in z.
            Lx = element.width;
            Ly = element.length;

            % Add rectangular element to the array
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
    %   ElemPos: The positions of the elements in the array (1 dimension).
    %   c0: Speed of sound.
    % Returns:
    %   source: The source object containing the mask and signals.
    %   time_delays: Time delays applied to each element for focusing.

    % Calculate time delays for each element based on steering angle
    
    [X, Y] = meshgrid(ElemPos, ElemPos);

    % Create the time delays for each element in the matrix array
    % time_delays0 = ( X.*sin(theta(2)) - Y.*sin(theta(1)) )/c0; % Plane wave (old version from Nikunj)
    time_delays0 = ( X.*sin(theta(2)) - Y.*sin(theta(1)).*cos(theta(2)) )/c0; % Plane wave (new Allen version)
    time_delays0 = time_delays0(:) - min(time_delays0(:)); % Shift so the lowest time delay is zero
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

