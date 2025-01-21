% ztest = 114;
ztest = 561;
xtest = 64;

% figure; scatter(Media.MP(:, 1), Media.MP(:, 3), '.')
figure; scatter3(Media.MP(:, 1), Media.MP(:, 2), Media.MP(:, 3), '.')

figure; scatter3(replacePoints(:, 1), replacePoints(:, 2), replacePoints(:, 3), '.')
% figure; scatter(temp_Media_sorted(:, 1), temp_Media_sorted(:, 3), '.')
% figure; scatter(replacePoints(:, 1), replacePoints(:, 3), '.')
% figure; scatter(tt(:, 1), tt(:, 3), '.')

% 3D pt plot
% pts = randomPts3D_func(30e-6, 30e-6, .02, .0001);
% figure; scatter3(pts(:, 1), pts(:, 2), pts(:, 3), '.')


xlabel('x (wl)')
ylabel('z (wl)')
% temp = abs(squeeze(IQ_coherent_sum(ztest, xtest, :)));

temp = real(squeeze(IQ_coherent_sum(ztest, xtest, :))); %%%%%%%

temp = temp - mean(temp); %%%
tstep = 1/fps_target; 
figure; plot(tstep:tstep:tstep*nf, temp); xlabel('time (s)'); ylabel('real(IQ coherent sum)')
fft_val = fft(temp);
fft_val = fftshift(fft_val);

faxis = linspace(-fps_target/2, fps_target/2, length(temp));
figure; plot(faxis, abs(fft_val)); ylabel('magnitude of FFT'); xlabel('frequency (Hz)')
fD_predicted = -2*Trans.frequency*1e6*flow_v_mm_s/1e3/Resource.Parameters.speedOfSound
xline(fD_predicted)
xline(-fD_predicted)
% figure; plot(fftest)