ztest = 140;
xtest = 145;

temp = real(squeeze(IQ_f(ztest, xtest, :))); %%%%%%%

temp = temp - mean(temp); %%%
tstep = 1/fps_target; 
figure; plot(tstep:tstep:tstep*nf, temp); xlabel('time (s)'); ylabel('real(SVD filtered IQ coherent sum)')
fft_val = fft(temp);
fft_val = fftshift(fft_val);

faxis = linspace(-fps_target/2, fps_target/2, length(temp));
figure; plot(faxis, abs(fft_val)); ylabel('magnitude of FFT'); xlabel('frequency (Hz)')
fD_predicted = -2*Trans.frequency*1e6*flow_v_mm_s/1e3/Resource.Parameters.speedOfSound
xline(fD_predicted)
xline(-fD_predicted)
% figure; plot(fftest)