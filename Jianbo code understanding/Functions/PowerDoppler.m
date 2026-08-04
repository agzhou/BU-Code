%% Power Doppler data processing 
function [PDI,PDINeg,PDIPos]=PowerDoppler(sIQ,Noise)
[nx,nz,nt]=size(sIQ);
sBlood=sIQ./repmat(Noise,[1,1,nt]);
PDI=mean(abs(sBlood).^2,3); 
%% directional PDI
nf= 2^nextpow2(2*nt);           % Fourier transform points
SpecBlood=fftshift(fft(sBlood,nf,3),3);
PDINeg=squeeze(sum(abs(SpecBlood(:,:,1:floor(nf/2)-1)).^2,3));
PDIPos=squeeze(sum(abs(SpecBlood(:,:,floor(nf/2)+1:nf)).^2,3));
