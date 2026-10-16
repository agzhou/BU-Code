IQ_OPW = IQ_RC + IQ_CR;

XD = IQ_RC .* conj(IQ_CR);

volumeViewer(abs(IQ_OPW));
volumeViewer(2.*real(XD));