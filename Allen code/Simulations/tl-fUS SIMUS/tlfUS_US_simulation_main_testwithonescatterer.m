clearvars
%%
addpath([cd, '\MUST'])

%% Set up parameters for the simulated probe
% param.fc = 13.8889e6; % Center frequency [Hz]
param.fc = 3e6; % Center frequency [Hz]
% param.fc = 1e6; % Center frequency [Hz]
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

% Choose angles about x and y
tiltX = 0; % (in rad)
tiltY = 0; % (in rad)

% Get transmit time delays
txdel = txdelay3(param, tiltX, tiltY); % in s

% Define volume grid
n = 32;
volume_grid.x_bounds = [-5e-3, 5e-3];
volume_grid.y_bounds = [-5e-3, 5e-3];
volume_grid.z_bounds = [0, 10e-3];
volume_grid.x = linspace(volume_grid.x_bounds(1), volume_grid.x_bounds(2), n); % (in m)
volume_grid.y = linspace(volume_grid.y_bounds(1), volume_grid.y_bounds(2), n); % (in m)
volume_grid.z = linspace(volume_grid.z_bounds(1), volume_grid.z_bounds(2), 4*n); % (in m)
[x, y, z] = meshgrid(volume_grid.x, volume_grid.y, volume_grid.z);

%% Simulate the RMS pressure field
% RP = pfield3(x,y,z,txdel,param);

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

%% Simulate backscattered RF signals
param.fs = 4*param.fc; % Sampling frequency [Hz]

% Define coordinates of scatterers (stored in struct ss)
ss.x = [0].* 1e-3; % x coordinates [m]
ss.y = [0].* 1e-3; % y coordinates [m]
ss.z = [5].* 1e-3; % z coordinates [m]
ss.Rc = ones(size(ss.x)); % Reflection coefficient = 1

RF = simus3(ss.x, ss.y, ss.z, ss.Rc, [txdel], param);

%% Look at 4 random RF signals (adapted from the example code)
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
  RF = RF/max(RF(:,n),[],'all');
  for k = 1:4
      nexttile
      plot((0:size(RF,1)-1)/param.fs*1e6,RF(:,n(k)))
      title(['Element #' int2str(n(k))])
      ylim([-1 1])
  end

%% IQ Demodulation and beamforming
IQ = rf2iq(RF, param);

%% Look at the same random 4 elements' IQ signals
%% Look at 4 random RF signals (adapted from the example code)
  %-- Display their IQ signals
  figure
  tl = tiledlayout(4,1);
  title(tl,'|IQ| signals')
  % RF = RF/max(RF(:,n),[],'all');
  for k = 1:4
      nexttile
      plot((0:size(IQ, 1)-1)/param.fs*1e6, abs(IQ(:, n(k))))
      title(['Element #' int2str(n(k))])
      % ylim([-1 1])
  end
%% Beamforming

lambda = param.c/param.fc;
% beamforming grid
bf.xvals = (volume_grid.x_bounds(1) : lambda : volume_grid.x_bounds(2))';
bf.yvals = (volume_grid.y_bounds(1) : lambda : volume_grid.y_bounds(2))';
bf.zvals = (volume_grid.z_bounds(1) : lambda : volume_grid.z_bounds(2))';
[bf.x, bf.y, bf.z] = meshgrid(bf.xvals, ...
                              bf.yvals, ...
                              bf.zvals);
figure; scatter3(bf.x, bf.y, bf.z, 20, 'filled')

IQbf = das3(IQ, bf.x, bf.y, bf.z, [txdel], param);
% IQbf = das3(IQ, bf.x(:), bf.y(:), bf.z(:), [txdel], param);

%%
% volshow(abs(IQbf)); axis square
figure; imagesc(squeeze(max(abs(IQbf), [], 1))')
figure; imagesc(squeeze(max(abs(IQbf), [], 2))')
figure; imagesc(squeeze(max(abs(IQbf), [], 3)))