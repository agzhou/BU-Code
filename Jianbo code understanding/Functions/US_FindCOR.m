function ret = US_FindCOR(RR)

	a = real(RR);
	b = imag(RR);

	A1 = std(a,1,3).^2;
	A2 = mean( (a - repmat(mean(a,3),[1 1 size(a,3) 1])) .* b ,3);
	A3 = mean( (b - repmat(mean(b,3),[1 1 size(b,3) 1])) .* a ,3);
	A4 = std(b,1,3).^2;
	
	B1 = 1/2 * ( mean(a.^3,3) - mean(a,3).*mean(a.^2,3) + mean( (a - repmat(mean(a,3),[1 1 size(a,3) 1])) .* b.^2 ,3) );
	B2 = 1/2 * ( mean(b.^3,3) - mean(b,3).*mean(b.^2,3) + mean( (b - repmat(mean(b,3),[1 1 size(b,3) 1])) .* a.^2 ,3) );

	ret = ( (A4.*B1 - A2.*B2) + 1i*(A1.*B2 - A3.*B1) ) ./ (A1.*A4 - A2.*A3);
    