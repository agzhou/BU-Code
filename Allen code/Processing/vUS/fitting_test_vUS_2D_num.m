%% Parameters for testing
v_xgp = 0.005; % m/s
v_zgp = 0.01; % m/s
tau = (0:1/5000:20e-3).'; % s
nTau = length(tau);
c0 = 1540; % m/s
fc = 15.625e6; % Hz
lambda0 = c0./fc; % m
k0 = 2*pi/lambda0; % m^-1
sigma = [58.9110, 73.9967].*1e-6; % Field-based 1/e PSF values (x, z) [m] for the L22-14v probe at 15.625 MHz and 17 angles from -10 to 10 deg. The y component is set to some arbitrary positive number but it won't really be used. (G:\My Drive\Data\PSF Simulations\L22-14v PSF sim - 17 angles from -10 to 10 deg)

F = 1; DC = 0;

%% Set the flow profile parameters and calculate the g1 model (single curve)
a = 1;
k = 2;
fmin = 0;
fmax = a^2 * (k + 2) ./ k;

% Calculate g1 model
g1 = zeros(length(tau), 1); % Initialize g1 vector
x = [v_xgp, v_zgp, F, DC, k, a]; % Parameter vector
% for ti = 1:length(tau)
%     t = tau(ti);
%     g1(ti) = integral(vUS_2D_num_vec(x, t, k0, sigma), fmin, fmax);
% end
g1(ti) = integral(vUS_2D_num_vec(x, tau, k0, sigma), fmin, fmax, 'ArrayValued', true);

% plot result
figure; plot(tau, abs(g1), 'LineWidth', 2), xlabel('tau'); ylabel('|g1|')

%% Sweep over the flow profile parameters (only k for now)
a = 1;
k = [2:0.2:3];
fmin = zeros(length(k), 1);
fmax = a^2 * (k + 2) ./ k;

% Calculate g1 model
g1 = cell(length(k), 1);

for ki = 1:length(k)
    x = [v_xgp, v_zgp, F, DC, k(ki), a]; % Parameter vector
    % temp_g1 = zeros(length(tau), 1); % Initialize g1 vector
    % for ti = 1:length(tau)
    %     t = tau(ti);
    %     temp_g1(ti) = integral(vUS_2D_num_vec(x, t, k0, sigma), fmin(ki), fmax(ki));
    % end
    % g1{ki} = temp_g1;
    g1{ki} = integral(vUS_2D_num_vec(x, tau, k0, sigma), fmin(ki), fmax(ki), 'ArrayValued', true);
end

% Plot
figure; hold on
for ki = 1:length(k)
    plot(tau.*1e3, abs(g1{ki}), 'LineWidth', 2), xlabel('tau [ms]'); ylabel('|g1|'); title('|g1| at different values of k')
end
legend(num2str(k.'))

figure; hold on
for ki = 1:length(k)
    plot(g1{ki}, 'LineWidth', 2), title('g1 at different values of k')
end
axis equal; xlim([-1, 1]); ylim([-1, 1]); 
legend(num2str(k.'), 'Location', 'southeast')