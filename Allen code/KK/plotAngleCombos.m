% Plot histograms of k = k_out - k_in for different receive angle schemes
% (in 2D case: k_x = k_x,out - k_x,in)
% And here, plotting delta_theta = theta_out - theta_in is equivalent for
% small angles.

%% Global parameters
f = 13.6e6; % Plane wave temporal frequency [Hz]
c = 1540; % Speed of sound (in tissue) [m/s]
%% Define TX angles ("in")

naTX = 21; % # of angles for TX
maTX = 10; % Max angle for TX [deg]

anglesTX = linspace(-maTX, maTX, naTX)'; % Transmit angles [deg]
daTX = mean(diff(anglesTX)); % TX angle increment [deg]

%% Case 1: Use the same TX angles for RX
anglesRX = linspace(-maTX, maTX, naTX)'; % Receive angles [deg]
delta_theta = repmat(anglesRX, 1, naTX) - repmat(anglesTX', length(anglesRX), 1);
k_x = f*2*pi .* sind(delta_theta) ./ c;
% figure; imagesc(delta_theta); colorbar; axis image; xlabel('TX angle'); ylabel('RX angle'); title('Delta theta [deg]')
% figure; imagesc(k_x); colorbar; axis image; xlabel('TX angle'); ylabel('RX angle'); title('k_x [radians/m]')

% Plot histograms
figure; histogram(delta_theta, BinMethod="integers"); title('Delta theta counts'); xlabel('Delta theta [deg]'); ylabel('Counts')
figure; histogram(k_x); title('Delta k_x'); xlabel('Delta k_x [radians/m]'); ylabel('Counts')

%% Case 2: Use the shifted angles for RX
naRX = 21;
o = fix(-naRX/2):1:fix(naRX/2); % Truncate towards zero
j = 7; % Shift parameter
anglesRX = (sign(o) .* daTX .* (2.*abs(o)./naRX + j))'; % Receive angles [deg]
% anglesRX = (sign(o) .* daTX .* (2.*abs(o) + j))'; % Receive angles [deg]
delta_theta = repmat(anglesRX, 1, naTX) - repmat(anglesTX', naRX, 1);
k_x = f*2*pi .* sind(delta_theta) ./ c;
% figure; imagesc(delta_theta); colorbar; axis image; xlabel('TX angle'); ylabel('RX angle'); title('Delta theta [deg]')
% figure; imagesc(k_x); colorbar; axis image; xlabel('TX angle'); ylabel('RX angle'); title('k_x [radians/m]')

% Plot histograms
figure; histogram(delta_theta, BinMethod="integers"); title("Delta theta counts, with j = " + num2str(j)); xlabel('Delta theta [deg]'); ylabel('Counts')
% figure; histogram(k_x); title("Delta k_x counts, with j = " + num2str(j)); xlabel('Delta k_x [radians/m]'); ylabel('Counts')
