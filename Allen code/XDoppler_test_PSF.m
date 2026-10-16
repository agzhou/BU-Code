IQ_OPW = IQ_RC + IQ_CR;

% XD = 2.*real(IQ_RC .* conj(IQ_CR));
XD = IQ_RC .* conj(IQ_CR) + IQ_CR .* conj(IQ_RC);

figure; imagesc(squeeze(max(abs(IQ_OPW), [], 1))')
figure; imagesc(squeeze(max(abs(XD), [], 1))')

volumeViewer(abs(IQ_OPW));
volumeViewer(XD);