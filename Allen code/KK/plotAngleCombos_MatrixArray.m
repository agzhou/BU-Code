% Plot histograms of k = k_out - k_in for different receive angle schemes
% (in 2D case: k_x = k_x,out - k_x,in)
% And here, plotting delta_theta = theta_out - theta_in is equivalent for
% small angles.

clearvars
%% Global parameters
f = 13.6e6; % Plane wave temporal frequency [Hz]
c = 1540; % Speed of sound (in tissue) [m/s]

binlimits = [-19.5, 19.5]; % Test histogram bin limits
binwidth = 1;
%% Define TX angles ("in")

naTX = 25; % # of angles for TX
maTX = 5; % Max angle for TX [deg]

anglesTXList = linspace(-maTX, maTX, naTX)'; % Transmit angles [deg]
[anglesTXX, anglesTXY] = meshgrid(anglesTXList, anglesTXList);
anglesTX = [anglesTXX(:), anglesTXY(:)];
daTX = maTX*2/(naTX - 1); % TX angle increment [deg]

% Plot the transmit angles
figure
scatter(anglesTX(:, 1), anglesTX(:, 2), 'o', 'filled'); axis image

%% Case 1: Use the same TX angles for RX
anglesRX = anglesTX; % Receive angles [deg]
delta_theta = calcDeltaTheta(anglesTX, anglesRX);

% k_x = f*2*pi .* sind(delta_theta) ./ c;
% figure; imagesc(delta_theta); colorbar; axis image; xlabel('TX angle'); ylabel('RX angle'); title('Delta theta [deg]')
% figure; imagesc(k_x); colorbar; axis image; xlabel('TX angle'); ylabel('RX angle'); title('k_x [radians/m]')

% Plot histograms
% figure; histogram2(delta_theta(:, 1), delta_theta(:, 2), BinMethod="integers"); title('Delta theta counts'); xlabel('Delta theta [deg]'); ylabel('Counts')
figure; histogram2(delta_theta(:, 1), delta_theta(:, 2), XBinLimits = binlimits, YBinLimits = binlimits, BinWidth = binwidth); title('Delta theta counts'); xlabel('Delta theta [deg]'); ylabel('Delta theta y [deg]'); zlabel('Counts')
% figure; histogram(k_x); title('Delta k_x'); xlabel('Delta k_x [radians/m]'); ylabel('Counts')

%% Case 2: Use the shifted angles for RX
% naRX = 5;
naRX = naTX; % Compare by using the same # of transmit angles
o = fix(-naRX/2):1:fix(naRX/2); % Truncate towards zero
j = 12; % Shift parameter
anglesRXList = (sign(o) .* daTX .* (2.*abs(o)./naRX + j))'; % Receive angles [deg]
anglesRX = listToAngles(anglesRXList);

delta_theta = calcDeltaTheta(anglesTX, anglesRX);

figure; histogram2(delta_theta(:, 1), delta_theta(:, 2), XBinLimits = binlimits, YBinLimits = binlimits, BinWidth = binwidth); title('Delta theta counts'); xlabel('Delta theta x [deg]'); ylabel('Delta theta y [deg]'); zlabel('Counts')

% figure; histogram(delta_theta(:, 1), BinMethod="integers")
% figure; histogram(delta_theta(:, 2), BinMethod="integers")

% delta_theta = repmat(anglesRX, 1, naTX) - repmat(anglesTX', naRX, 1);
% k_x = f*2*pi .* sind(delta_theta) ./ c;
% % figure; imagesc(delta_theta); colorbar; axis image; xlabel('TX angle'); ylabel('RX angle'); title('Delta theta [deg]')
% % figure; imagesc(k_x); colorbar; axis image; xlabel('TX angle'); ylabel('RX angle'); title('k_x [radians/m]')
% 
% % Plot histograms
% figure; histogram(delta_theta, BinMethod="integers"); title("Delta theta counts, with j = " + num2str(j)); xlabel('Delta theta [deg]'); ylabel('Counts')
% % figure; histogram(k_x); title("Delta k_x counts, with j = " + num2str(j)); xlabel('Delta k_x [radians/m]'); ylabel('Counts')

%% Helper functions
function angles = listToAngles(anglesList) % Turn the single-dimension angle list into a full list of x and y angles
    [anglesX, anglesY] = meshgrid(anglesList, anglesList);
    angles = [anglesX(:), anglesY(:)];
end

function delta_theta = calcDeltaTheta(anglesTX, anglesRX)
    delta_theta = [];
    for ind = 1:size(anglesTX, 1)
        delta_theta = [delta_theta; anglesRX - anglesTX(ind, :)];
    
    end
end