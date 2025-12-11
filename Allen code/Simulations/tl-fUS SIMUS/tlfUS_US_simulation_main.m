clearvars
%%
addpath([cd, '\MUST'])

%% Set up parameters for the simulated probe
% param.fc = 13.8889e6; % Center frequency [Hz]
% param.fc = 3e6; % Center frequency [Hz]
param.fc = 100; % Center frequency [Hz]
param.bandwidth = 70; % Bandwidth [% of center frequency]
param.width = 250e-6; % Array element width (x) [m]
param.height = 250e-6; % Array element height (y) [m]

% x/y coords of each element
pitch = 300e-6; % pitch (in m)
[xe,ye] = meshgrid(((1:32)-16.5)*pitch);
param.elements = [xe(:).'; ye(:).'];

% Transmit apodization
h = cos(linspace(-pi/4,pi/4,32));
h = h'*h;
param.TXapodization = h(:);

% Choose angles about x and y
tiltX = 0; % (in rad)
tiltY = 0; % (in rad)

% Get transmit time delays
txdel = txdelay3(param, tiltX, tiltY); % in s

% Define volume grid
n = 32;
volume_grid.x = linspace(-5e-3, 5e-3, n); % (in m)
volume_grid.y = linspace(-5e-3, 5e-3, n); % (in m)
volume_grid.z = linspace(0, 10e-3, 4*n); % (in m)
[x, y, z] = meshgrid(volume_grid.x, volume_grid.y, volume_grid.z);

%%
% Simulate the RMS pressure field
RP = pfield3(x,y,z,txdel,param);

%% Display the acoustic pressure field and element centers
RPdB = 20*log10(RP/max(RP(:))); % convert to dB
slice(x*1e2,y*1e2,z*1e2,RPdB,0,0,1:5)
shading flat
colormap(hot), caxis([-6 0])
set(gca,'zdir','reverse'), axis equal
% Fine-tune the figure:
alpha color % some transparency
c = colorbar; c.YTickLabel{end} = '0 dB';
zlabel('[cm]')
title('Plane wave - RMS pressure field')
hold on
plot3(xe*1e2,ye*1e2,xe*0,'b.')
hold off

%% Simulate backscattered RF signals

% Define coordinates of scatterers (stored in struct ss)
ss.x = [0].* 1e3; % x coordinates [m]
ss.y = [0].* 1e3; % y coordinates [m]
ss.z = [5].* 1e3; % z coordinates [m]
ss.Rc = ones(size(ss.x)); % Reflection coefficient = 1

RF = simus3(ss.x, ss.y, ss.z, ss.Rc, [txdel], param);
