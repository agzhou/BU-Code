%% g1 based Vz calcualtion, CPU
% formula: 2pi/T=2n*2pi/Lambda*Vz => Vz=Lambda/(2nT), n is the optical refractive index
% input: 
    % GG, 3D array, (nz,nx,nTau)
    % PRSSinfo: data processing parameters, including 
        % PRSSinfo.rFrame: sIQ frame rate, Hz
        % PRSSinfo.f0: Transducer center frequency, Hz
        % PRSSinfo.C: Sound speed in the sample, m/s
        % PRSSinfo.Dim: GG original dimension, [nz,nx,nTau]
        % PRSSinfo.rfnScale: spatial refind scale
% output: 
    % Vz, 2D, [nz,nx], m/s
%%%%%%%%%%%%%%% EXAMPLE %%%%%%%%%%%%%%%%%%%%%%

function [Vz]=GG2Vz_FREQ(GG, PRSSinfo)

[nz,nx,nt]=size(GG);
fCoor=linspace(-PRSSinfo.rFrame/2,PRSSinfo.rFrame/2,nt);
PN=sign(mean(mean(mean(imag(GG(:,:,1:3)),3))));
if PN==1
    fCoor(fCoor<0)=0;
else
    fCoor(fCoor>0)=0;
end

fBlood=fftshift(fft(GG,nt,3),3);
% fBlood(abs(fBlood)<1*std(abs(fBlood(abs(fBlood)>0))))=0; % thresholding
% fBlood(abs(fBlood)<4*mean(abs(fBlood(abs(fCoor)>1000))))=0; % thresholding
fD=sum(repmat(permute(fCoor,[1,3,2]),[nz,nx,1]).*abs(fBlood).^2,3)./sum(abs(fBlood).^2,3);
Vz0=fD.*PRSSinfo.C/(2*PRSSinfo.f0); % axial speed obtained with color Doppler, mm/s
gR=mean(abs(real(GG(:,:,1:2))),3);
gRm=mean(gR(:));
gRstd=std(gR(:));
gCR=(gR>(max((gRm-0.1*gRstd),0.12))*1)+((gR<=(max((gRm-0.1*gRstd),0.12))).*gR*3);
% gCR=imgaussfilt(gCR,0.8);
gCR=medfilt2(gCR,[3, 3]);
Vz=Vz0.*gCR;
Vz(isnan(abs(Vz)))=0;

