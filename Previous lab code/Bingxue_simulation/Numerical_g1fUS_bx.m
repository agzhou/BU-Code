function [IQ, FIQ, g1,Numer, tCoor, fCoor, tauCoor, Noise, nparticle] = Numerical_g1fUS_bx(Vtot, WidthVes, Noise,agl)
% snr = 20 dB
% Vtot = -10e-3; %m/s
% nScatter = 10
% WidthVes = 20

% rng(2,'twister');

%% ultrasound numerical simulation - Multiple vessels
% ith Particle’s initial position: (xi0,yi0,zi0)
% Particle moving speed: (vxi,vyi,vzi )
% Central position of the resolution voxel: (x0,y0,z0)
% time dependent position of the scatter: (xit, yit, zit)
% spatial resolution: (Rx,Ry,Rz)
% Ns: number of static particle
% Nf1: number of dynamic particle flowing in direction 1
% clear all;
% clc;
%addpath('D:\OneDrive\Work\PROJ - FUS\CODE\Functions') % Path on JTOPTICS
%% 0. system parameter
pixsize = 0.1; % in um
f0=15e6; % center frequency, Hz
C=1540; % sound speed, m/s
k0=2*pi*f0/C/1; 
z0=0e-6;  x0= 0e-6;  y0=0;
A_FWHM=[110 300 85]*1*1e-6; % [x,y,z]spatial resolution, FWHM, m 110 [110 300 85]
% Sigma=R_FWHM/(2*sqrt(2*log(2)));
% Sigma=R_FWHM/(exp(1));
Sigma=A_FWHM*1/(2*sqrt(2*log(2)));
fRate=5000; % frame rate, HZ
dt=1/fRate; 
nlength=400;% length for vessel in vetical, um 
nt0 = 200; % 150 for 10mm/s 5000Hz run 300 um: 150/5000 *10*1e3 = 300 um 
nt=20000+nt0; %  #of time points 4s
nTau=100;
nTauFit=nTau;
nNewlength = abs(Vtot)*10*dt*1e6; % vessel length for concantented parts: use 10 
%% 1. imaging field
pixZ=pixsize; % pixel size in Z, um
pixX=pixsize; % pixel size in X, um
pixY=pixsize; % pixel size in Y, um
zCoor=[-500:pixZ:500]*1e-3; 
xCoor=[-500:pixX:500]*1e-3;  
yCoor=[-0:pixY:0]*1e-3; 
nz=length(zCoor);
nx=length(xCoor);
ny=length(yCoor);
%% 2. Vessel parameters
% 2.1 parameters for vessel 1, flowing down
AngleVes = agl; % vessel 1 orientation, angle degree
WidthVes0 = WidthVes;   % WidthVes um vessel width, in pixel, real vessel width=WidthVes*pixX
density= 1/20; % 1/um^2 number of RBCs per cross section for each vessel
vzVes=Vtot*cos(AngleVes/180*pi); % axial flow speed in Vessel 1, m/s
vxVes=Vtot*sin(AngleVes/180*pi); % Lateral-X flow speed in Vessel 1, m/s
V1z0=0; V1x0=0; V1y0=0;  % vessel 1 center location, um

%% 3. vessel location
[xP0,zP0] = GenRdPtcl(WidthVes0, nlength, pixsize, density, 0);
   
figure,plot(xP0,zP0,'.'); axis equal tight; grid on
 xlim([-200 200]); ylim([-200 200]);
xlabel('x [um]'); ylabel('z [um]');

%% 4. time series signal 
% 4.1 RBC flow speed
IQ0=zeros(1,1, nt);
xP = xP0; zP = zP0;
[xPr, zPr] = RotatVes(xP,zP,AngleVes);
data = exp(-(xPr*1e-6-x0).^2/(2*Sigma(1)^2)).*exp(-(zPr*1e-6-z0).^2/(2*Sigma(3)^2)).*exp(2i*k0*(zPr*1e-6-z0)); 
IQ0(1,1,1)=sum((data(:,:)));%...1:round(end/2)
%     figure; plot(abs(data));
nparticle(:,1) = sum(abs(data)>1e-4);
 tic
for it=2:nt
    %% generate new RBC every v*dt um for Ves 
%      [xNewP, zNewP] = GenRdPtcl(WidthVes0, nlength, pixsize, density, AngleVes);  
     [xNewP, zNewP] = GenRdPtcl(WidthVes0, nNewlength, pixsize, density, 0);
%      rP = sqrt(xP.^2+zP.^2).*sign(zP); %sign(zP)
%      rNewP = sqrt(xNewP.^2+zNewP.^2).*sign(zNewP);%sign(zNewP)
     indKept = find(zP<=nlength/2-abs(Vtot)*dt*1e6);
     indIn = find(-zNewP>nNewlength/2-abs(Vtot)*dt*1e6);
     xP = [xP(indKept),xNewP(indIn)];%
     zP = [zP(indKept)-Vtot*dt*1e6, zNewP(indIn)-(nlength-nNewlength)/2];%
     [xPr,zPr] = RotatVes(xP,zP,AngleVes);
    %% IQ data for the resolution voxel/pixel
    data = exp(-(xPr*1e-6-x0).^2/(2*Sigma(1)^2)).*exp(-(zPr*1e-6-z0).^2/(2*Sigma(3)^2)).*(exp(2i*k0*(zPr*1e-6-z0))); 
     IQ0(1,1,it)=sum((data(:,:)));%...1:round(end/2)
%     figure; plot(abs(data));
    nparticle(:,it) = sum(abs(data)>1e-4);
    nall(:,it) = length(data);
end
 toc
nparticle = mean(nparticle((nt0+1):end));
% nparticle = mean(nall((nt0+1):end));
nt = nt-nt0;
IQ(1,1,:) = IQ0(1,1,(nt0+1):end); 

% figure; plot(nall(2:end))
% figure; plot(nparticle)
% IQ = IQ.*conj(IQ); xlim = 200;
[vx,vz] = meshgrid([-xlim*1e-6:pixsize*1e-6:xlim*1e-6]+x0,[-xlim*1e-6:pixsize*1e-6:xlim*1e-6]+z0); 
val0 = exp(-(vx-x0).^2/(2*Sigma(1)^2)).*exp(-(vz-z0).^2/(2*Sigma(3)^2)).*exp(2i*k0*(vz-z0));    
% figure; imagesc([-xlim:pixsize:xlim],[-xlim:pixsize:xlim],abs(val0));axis image; xlabel('x[um]');ylabel('z[um]');colorbar;
% figure; imagesc([-xlim:pixsize:xlim],[-xlim:pixsize:xlim],real(val0));axis image; xlabel('x[um]');ylabel('z[um]');colorbar;set(gca,'YTick',[-100,-75,-50,-25,0,25,50,75,100]);alpha(h1, abs(real(val0)));
% figure; imagesc([-xlim:pixsize:xlim],[-xlim:pixsize:xlim],imag(val0));axis image; xlabel('x[um]');ylabel('z[um]');colorbar;
% figure; imagesc([-xlim:pixsize:xlim],[-xlim:pixsize:xlim],angle(val0));axis image; xlabel('x[um]');ylabel('z[um]');colorbar;
% figure;subplot(141); plot([-xlim:pixsize:xlim],abs(val0(:,xlim/pixsize)));subplot(142); plot([-xlim:pixsize:xlim],real(val0(:,xlim/pixsize)));subplot(143); plot([-xlim:pixsize:xlim],imag(val0(:,xlim/pixsize)));subplot(144); plot([-xlim:pixsize:xlim],angle(val0(:,xlim/pixsize)));
% figure;subplot(141); plot([-xlim:pixsize:xlim],abs(val0(xlim/pixsize,:)));subplot(142); plot([-xlim:pixsize:xlim],real(val0(xlim/pixsize,:)));subplot(143); plot([-xlim:pixsize:xlim],imag(val0(xlim/pixsize,:)));subplot(144); plot([-xlim:pixsize:xlim],angle(val0(xlim/pixsize,:)));
% figure;subplot(141); plot([-xlim:pixsize:xlim],abs(diag(val0)));subplot(142); plot([-xlim:pixsize:xlim],real(diag(val0)));subplot(143); plot([-xlim:pixsize:xlim],imag(diag(val0)));subplot(144); plot([-xlim:pixsize:xlim],angle(diag(val0)));
% 5. plot and save particle flow animation
% nFrame=100;
% fig=figure;
% for it=1:nFrame %floor(nFrame/20):nFrame
%     plot(V1Ptcl_xP(:,it)*1e6,V1Ptcl_zP(:,it)*1e6,'b.'); 
%     hold on, plot(V2Ptcl_xP(:,it)*1e6,V2Ptcl_zP(:,it)*1e6,'r.'); 
%     hold off
%     axis equal;
%     xlim([-200 200]); ylim([-200 200]);
%     xlabel('x [um]')
%     ylabel('z [um]')
%     title (['t=',num2str(it*dt*1e3),' ms'])
%     grid on
%     drawnow; pause(0.1);
%     frames(it)=getframe(gcf);
% end
% FilePath= 'D:\g1_based_fUS\NumericalSimulation\';
% FileName='FlowParticle';
% outfileGIF=[FilePath,FileName,'1.gif'];
% iFile=1;
% while exist(outfileGIF)==2
%     iFile=iFile+1;
%     outfileGIF=[FilePath,FileName,num2str(iFile),'.gif'];
% end
% close(fig)  
% for it=1:nFrame
%     % save GIF 
%     im = frame2im(frames(it));
%     [imind,cm] = rgb2ind(im,256);  
%     if it==1
%         imwrite(imind,cm,outfileGIF,'gif','DelayTime',0.1,'loopcount',inf);
%     else
%         imwrite(imind,cm,outfileGIF,'gif','DelayTime',0.1,'writemode','append');
%     end
% end
%% add noise
if nargin < 3 
    IQ0 = squeeze(IQ);
    snr = 1;
    IQn = awgn(squeeze(IQ0),snr,'measured','db');
    IQ(1,1,:) = IQn;
    Noise = IQn - IQ0;
else
    if max(size(Noise)) ~= 1
    IQ0 = squeeze(IQ);
    IQ(1,1,:) = squeeze(IQ0)+Noise; 
    else
    IQ0 = squeeze(IQ);
    snr = Noise;
    IQn = awgn(squeeze(IQ0),snr,'measured','db');
    IQ(1,1,:) = IQn;
    Noise = IQn - IQ0;
    end
end
% check real snr
Signal_rms = sqrt(mean(abs(IQ0(:)).^2));
Noise_rms = sqrt(mean(abs(Noise(:)).^2));
snr_real = 20 .* log10(Signal_rms ./ Noise_rms)
% print({'SNR is', num2str(snr_real),'dB'})
% %% SVD or HP filter
% SignalRank = [2:nt];
% [sIQ, NoiseIQ]=SVDfilter(IQ,SignalRank);
% %%
% IQ = sIQ;
% %%
% IQ = IQ-mean(IQ,"all");% matlab R2021b
IQ = IQ-sum(IQ(:))/length(squeeze(IQ));% matlab R2017a
%% 6. raw data and frequency spectrum
Lfft=2^nextpow2(nt);
FIQ=fftshift(fft((squeeze(IQ)),Lfft));

tCoor=linspace(dt,nt*dt,nt)*1e3;
fCoor=linspace(-fRate/2, fRate/2,Lfft);
figure,
subplot(2,1,1),plot(abs(squeeze(IQ)));
xlabel('Time[ms]')
ylabel('sIQ')
subplot(2,1,2),plot(fCoor,abs(FIQ));
%xlim([-500 500])
xlabel('Frequency [Hz]')
ylabel('Power');
%% 7. simulated and fitted g1
% 7.1 g1 
[g1, Numer] =IQ2g1(IQ,1,nt,nTau);
tauCoor=linspace(dt,nTau*dt,nTau)*1e3;
% %plot g1
% fig=figure;
% set(fig,'Position',[400 400 800 650])
% subplot(2,2,1);hold on, plot(tauCoor,real(squeeze(g1)));title('real(g1), RAW'); xlabel('time lag, [ms]')
% subplot(2,2,2);hold on, plot(tauCoor,imag(squeeze(g1)));title('imag(g1), RAW');xlabel('time lag, [ms]')
% subplot(2,2,3);hold on, plot(squeeze(g1).');title('cplx(g1), RAW');%ylim([-1,1]); xlim([-1,1])
% subplot(2,2,4);hold on, plot(tauCoor,abs(squeeze(g1)));title('abs(g1), RAW'); 
% xlabel('time lag, [ms]')
end