clearvars
%%
addpath([cd, '\MUST'])
addpath([cd, '\..'])

%% Set up parameters for the simulated probe
% param.fc = 13.8889e6; % Center frequency [Hz]
param.fc = 3e6; % Center frequency [Hz]
% param.fc = 8e6;
param.fs = 4*param.fc; % Sampling frequency [Hz]f
param.bandwidth = 70; % Bandwidth [% of center frequency]
param.width = 250e-6; % Array element width (x) [m]
param.height = 250e-6; % Array element height (y) [m]
param.c = 1540; % Speed of sound [m/s]

% x/y coords of each element
ne = 32; % # of elements in one dimension, for ne x ne matrix
pitch = 300e-6; % pitch (in m)
% [xe, ye] = meshgrid(((1:32)-16.5)*pitch);
[xe, ye] = meshgrid(((1:ne) - (ne + 1)/2)*pitch);
param.elements = [xe(:).'; ye(:).'];

% Transmit apodization
h = cos(linspace(-pi/4, pi/4, ne));
h = h'*h;
param.TXapodization = h(:);

% Define volume grid
volume_grid.x_bounds = [-5e-3, 5e-3];
volume_grid.y_bounds = [-5e-3, 5e-3];
volume_grid.z_bounds = [0, 5e-3];

%% Choose tilt angles (about x and y)
% With a matrix array, we get na x na total transmits/plane waves
% TX.na = 5; % # of angles per axis
TX.na = 3; % # of angles per axis
% TX.na = 1; % # of angles per axis
TX.nta = TX.na^2; % # of total angles/transmissions
TX.ma = 5 * pi/180; % max angle [rad]
TX.angles = linspace(-TX.ma, TX.ma, TX.na); % Angles for one axis
[TX.tiltX, TX.tiltY] = meshgrid(TX.angles, TX.angles); % Get all the combinations of tilt angles
% Vectorize the tilt angles
TX.tiltX = TX.tiltX(:);
TX.tiltY = TX.tiltY(:);
% figure; scatter(TX.tiltX, TX.tiltY); axis square; xlabel('x angle [rad]'); ylabel('y angle [rad]') % Plot all the tilt angles used

% Get transmit time delays
TX.txdel = cell(TX.nta, 1);
for ti = 1:TX.nta % Go through each transmission index
    TX.txdel{ti} = txdelay3(param, TX.tiltX(ti), TX.tiltY(ti)); % [s]
end

%% Simulate the RMS pressure field

% n = 32;
% volume_grid.x = linspace(volume_grid.x_bounds(1), volume_grid.x_bounds(2), n); % (in m)
% volume_grid.y = linspace(volume_grid.y_bounds(1), volume_grid.y_bounds(2), n); % (in m)
% volume_grid.z = linspace(volume_grid.z_bounds(1), volume_grid.z_bounds(2), 4*n); % (in m)
% [x, y, z] = meshgrid(volume_grid.x, volume_grid.y, volume_grid.z);
% RP = pfield3(x,y,z,TX.txdel,param);

%% Display the acoustic pressure field and element centers
% RPdB = 20*log10(RP/max(RP(:))); % convert to dB
% slice(x*1e2,y*1e2,z*1e2,RPdB,0,0,1:5)
% shading flat
% colormap(hot), caxis([-6 0])
% set(gca,'zdir','reverse'), axis equal
% % Fine-tune the figure:
% alpha color % some transparency
% c = colorbar; c.YTickLabel{end} = '0 dB';
% zlabel('[cm]')
% title('Plane wave - RMS pressure field')
% hold on
% plot3(xe*1e2,ye*1e2,xe*0,'b.')
% hold off



%% Define coordinates of scatterers (stored in struct ss)
% Define Simulation Parameter struct
% SP.c = 1540; % Speed of sound [m/s]
SP.c = param.c; % Speed of sound [m/s]
% SP.f = 13.8889 * 1e6; % Ultrasound frequency [Hz]
SP.f = param.fc; % Ultrasound frequency [Hz]
SP.wl = SP.c / SP.f; % Wavelength [m]
SP.frameRate = 2500; % Frame rate [Hz]
SP.scatterReflectivity = 1.0;
SP.sigma = [300e-6, 300e-6, 150e-6]; %%%% PSF testing %%%%

SP.snr = 50; % Choose the SNR for the data vs. Gaussian white noise (5 is what Bingxue used)

% SP.vesselDiam = 50e-6; % Vessel diameter [m]
SP.vesselDiam = 100e-6; % Vessel diameter [m]

SP.vesselLength = 5/1e3;  % Vessel length [m]

% Define the center of the vessel [m]
SP.xstart = 0;
SP.ystart = 0;
SP.zstart = 3 * 1e-3;

% Get points in a cylindrical "vessel"
% cyl_vessel = genRandomPts3D_cyl(vesselDiam, vesselLength, startDepthMM/1e3, xstart, ystart, zstart);
[cyl_vessel, SP] = genRandomPts3D_cyl(SP);
% plotPoints(cyl_vessel, SP)
SP.dim = 3;

SP.flow_v_mm_s = 30;
SP.flow_dim = 3; %%%%%%%%

% ss.x = cyl_vessel(:, 1); % x coordinates [m]
% ss.y = cyl_vessel(:, 2); % y coordinates [m]
% ss.z = cyl_vessel(:, 3); % z coordinates [m]

% ---- Rotate the vessel ---- %

% Make the vessel horizontal
xa = 90;
ya = 0;
za = 0;
cyl_vessel_rot = rotateVessel(cyl_vessel, xa, ya, za, SP);




% plotPoints(cyl_vessel, SP)
% plotPoints(cyl_vessel_rot, SP)

% Set the scatterer points for SIMUS
ss.x = cyl_vessel_rot(:, 1); % x coordinates [m]
ss.y = cyl_vessel_rot(:, 2); % y coordinates [m]
ss.z = cyl_vessel_rot(:, 3); % z coordinates [m]
ss.Rc = cyl_vessel(:, 4); % Reflection coefficients

%% %-- Display the elements and the scatterers
  figure
  scatter3(ss.x*1e3, ss.y*1e3, ss.z*1e3, 30, 'filled')
  colormap(cool)
  hold on
  scatter3(xe*1e3,ye*1e3,0*xe,3,'b','filled')
  axis equal, box on
  set(gca,'zdir','reverse')
  zlabel('[mm]')
  title([int2str(size(cyl_vessel, 1)) ' scatterers'])

%% Simulate backscattered RF signals
RF = cell(TX.nta, 1);

sim_options.ParPool = true; % Enable parallel computing
for ti = 1:TX.nta % Go through each transmission index
    RF{ti} = simus3(ss.x, ss.y, ss.z, ss.Rc, TX.txdel{ti}, param, sim_options);
end

%% Look at 4 random RF signals (adapted from the example code)

% Use RF cell #:
test_cell = 1;

%-- Choose 4 elements randomly
  n = sort(randi(ne^2, 1, 4));
  plot(xe(n)*1e3,ye(n)*1e3,'ro','MarkerFaceColor','r')
  for k = 1:4
      text(xe(n(k))*1e3+0.3,ye(n(k))*1e3,int2str(n(k)),...
          'Color','r','BackgroundColor','w')
  end
  title('32{\times}32 elements')
  %-- Display their RF signals
  figure
  tl = tiledlayout(4,1);
  title(tl,'RF signals')
  RF_test = RF{test_cell}/max(RF{test_cell}(:, n),[], 'all');
  for k = 1:4
      nexttile
      plot((0:size(RF_test, 1) - 1)/param.fs*1e6, RF_test(:, n(k)))
      title(['Element #' int2str(n(k))])
      ylim([-1 1])
  end

%% IQ Demodulation and beamforming
IQ = cell(TX.nta, 1);
for ti = 1:TX.nta % Go through each transmission index
    IQ{ti} = rf2iq(RF{ti}, param); % in s
end

%% Look at the same random 4 elements' IQ signals
% Look at 4 random RF signals (adapted from the example code)
  %-- Display their IQ signals
  figure
  tl = tiledlayout(4,1);
  title(tl,'|IQ| signals')
  for k = 1:4
      nexttile
      plot((0:size(IQ{test_cell}, 1)-1)/param.fs*1e6, abs(IQ{test_cell}(:, n(k))))
      title(['Element #' int2str(n(k))])
      % ylim([-1 1])
  end

%% Beamforming with the matrix testing
lambda = param.c/param.fc;
% beamforming grid
bf.xvals = (volume_grid.x_bounds(1) : lambda/2 : volume_grid.x_bounds(2))';
bf.yvals = (volume_grid.y_bounds(1) : lambda/2 : volume_grid.y_bounds(2))';
bf.zvals = (volume_grid.z_bounds(1) : lambda/2 : volume_grid.z_bounds(2))';
[bf.x, bf.y, bf.z] = meshgrid(bf.xvals, ...
                              bf.yvals, ...
                              bf.zvals);
% figure; scatter3(bf.x, bf.y, bf.z, 20, 'filled')

% IQbf = cell(TX.nta, 1);

M = dasmtx3([size(IQ{1}, 1) size(IQ{1}, 2)], bf.x, bf.y, bf.z, TX.txdel{1}, param);
%% Beamforming

lambda = param.c/param.fc;
% beamforming grid
bf.xvals = (volume_grid.x_bounds(1) : lambda/2 : volume_grid.x_bounds(2))';
bf.yvals = (volume_grid.y_bounds(1) : lambda/2 : volume_grid.y_bounds(2))';
bf.zvals = (volume_grid.z_bounds(1) : lambda/2 : volume_grid.z_bounds(2))';
[bf.x, bf.y, bf.z] = meshgrid(bf.xvals, ...
                              bf.yvals, ...
                              bf.zvals);
% figure; scatter3(bf.x, bf.y, bf.z, 20, 'filled')

IQbf = cell(TX.nta, 1);
tic
for ti = 1:TX.nta % Go through each transmission index
    IQbf{ti} = das3(IQ{ti}, bf.x, bf.y, bf.z, TX.txdel{ti}, param);
end
toc
%% Combine all the transmissions into one coherently compounded volume
IQbf_cpwc = IQbf{1}; % Initialize the CPWC IQ volume
if TX.nta > 1 % If there is more than one transmission
    for ti = 2:TX.nta % Go through each transmission index
        IQbf_cpwc = IQbf_cpwc + IQbf{ti};
    end
end
% volshow(abs(IQbf_cpwc)); axis square
figure; imagesc(squeeze(max(abs(IQbf_cpwc), [], 1))'); title('xz MIP')
figure; imagesc(squeeze(max(abs(IQbf_cpwc), [], 2))'); title('yz MIP')
figure; imagesc(squeeze(max(abs(IQbf_cpwc), [], 3))); title('xy MIP')

% figure; imagesc(squeeze(max(abs(IQbf{1}), [], 1))')
% figure; imagesc(squeeze(max(abs(IQbf{1}), [], 2))')
% figure; imagesc(squeeze(max(abs(IQbf{1}), [], 3)))