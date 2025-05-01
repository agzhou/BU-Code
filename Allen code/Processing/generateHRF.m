
pps = P.daqrate; % Points per second
% All below values are in [s] and converted according some sample rate
delay_response = 1.5 * pps;
delay_undershoot = 10 * pps;
dispersion_response = 0.5 * pps;
dispersion_undershoot = 1 * pps;
ratio_response_to_undershoot = 6 * pps;
onset = 0 * pps;
kernel_length = 16 * pps;

duration = 30 * pps;

%%
% temp to figure out what it does
b = 1/10;
c = 1;

t = 0:duration;
% HRF = 1 .* (t - delay_response).^2 .* exp(-b .* (t - delay_response).^c); % .* heaviside(t - delay_response)
% HRF = gampdf(t, 3, 5);

alpha = 5;
lambda = 0.001;
% gamma = t .^ (alpha-1) .* exp(-lambda .* t) .* lambda .^ alpha ./ factorial(alpha-1);
% HRF = (t - delay_response) .^ (alpha-1) .* exp(-lambda .* (t - delay_response)) .* lambda .^ alpha ./ factorial(alpha-1);

% figure; plot(gamma)

%%
GF = @(t, alpha, lambda) t .^ (alpha-1) .* exp(-lambda .* t) .* lambda .^ alpha ./ factorial(alpha-1);

alpha = 4;
lambda = 0.005;

stim_1trial = TD.airPuffOutput(trial_windows{1});
gamma = GF(1:length(stim_1trial), alpha, lambda);
% HRF = conv(stim_1trial, gamma, 'same');
HRF = conv(stim_1trial, gamma);
%
figure
plot(stim_1trial)
hold on
plot(gamma)
plot(HRF)
hold off
