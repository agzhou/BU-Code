%% US g1 fit, fit all frequency signal, for single flow direction
% For algorithm development, data analysis
% input: IQ(SVD filtered) 
% g1Info (g1 calculation parameters, including g1Info.tStart, g1Info.nt, g1Info.nTau)
% A_FWHM (X, Y, Z spatial resolution, um)
% rFrame (IQ frame rate, Hz)
function [Ms, Mf, Vx, Vy, Vz,Pvz,Vcz,R,GG, PDI]=sIQ2vUS_SglFD(sIQ, g1Info, A_FWHM, rFrame, rfnScale, nItpT)
if nargin<5
    rfnScale=1;
    nItpT=1;
end
%% I. refine sIQ
[nz0,nx0,nt]=size(sIQ);
cIQ=zeros(nz0*rfnScale,nx0*rfnScale,nt);
if rfnScale>1
    for it=1:nt
        cIQ(:,:,it)=imresize(sIQ(:,:,it),[nz0,nx0]*rfnScale); % spatial interpolation
    end
else
    cIQ=sIQ;
end
clear sIQ
% Power Doppler
[PDI0]=sIQ2PDI(cIQ);
PDI=squeeze(PDI0(:,:,3));
%% II. constant
[nz,nx,nt]=size(cIQ);
Ms=zeros(nz,nx); % 1-positive frequency, 2-negative frequency
Vx=Ms; Vy=Ms; Vz=Ms;  Pvz=Vz; Vcz=zeros(nz,nx);
Mf=zeros(nz,nx,2);
R=zeros(nz,nx);
GG=zeros(nz,nx,g1Info.nTau*nItpT,2); %  dimension 5-[gg,ggfit]
rFrameItp=rFrame*nItpT;
f0=16.625e6; % ultrasound transducer center frequency, Hz
%% III. determine the signal (|f|<800Hz) to noise ratio
fCoor=linspace(-rFrame/2,rFrame/2,nt)';
fCoorSig=zeros(size(fCoor));
fCoorSig(abs(fCoor)<800)=1; % signal frequency range
fCoorSig=circshift(fCoorSig,nt/2);
fIQ=(fft(cIQ,nt,3)); % no fft shift
SNR=squeeze(sum(abs(fIQ.*repmat(permute(fCoorSig,[3 2 1]),[nz nx 1])),3))./squeeze(sum(abs(fIQ),3)); % SNR of oringla data

%% IV. Fitting criteria
warning off;
for iz=1:nz
    for ix=1:nx
        %% IV.1 fitting criteria - SNR>0.45-0.001*iz
        if SNR(iz,ix)>mean(SNR(:))+1.5*std(SNR(:))
            GG(iz,ix,:,1)=interp(squeeze(IQ2g1(cIQ(iz,ix,:),g1Info.tStart,g1Info.nt,g1Info.nTau)),nItpT); 
            [Ms(iz,ix), Mf(iz,ix,:), Vx(iz,ix), Vy(iz,ix), Vz(iz,ix),Pvz(iz,ix),R(iz,ix),GG(iz,ix,:,2)]=USg1fit_SV(GG(iz,ix,:,1), A_FWHM, rFrameItp);
            %% IV.2 color Doppler
            [Vcz(iz,ix)]=ColorDoppler(cIQ(iz,ix,:),f0,rFrame); % color Doppler, all frequency
        end
    end
end
Vz=-1*Vz;
Vcz=-1*Vcz;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
