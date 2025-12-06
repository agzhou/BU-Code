function [IQ, FIQ, g1,Numer, tCoor, fCoor, tauCoor, Noise, nparticle] = Numerical_g1fUS(Vtot, WidthVes, Noise)
% snr = 20 dB
% Vtot = -10e-3 m/s
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
pixsize = 0.1;
f0=15e6; % center frequency, Hz
C=1540; % sound speed, m/s
k0=2*pi*f0/C/1; 
z0=0.5e-6;  x0=-300e-6;  y0=0;
A_FWHM=[110 300 110]*1*1e-6; % [x,y,z]spatial resolution, FWHM, m 110 [110 300 85]
% Sigma=R_FWHM/(2*sqrt(2*log(2)));
% Sigma=R_FWHM/(exp(1));
Sigma=A_FWHM*1/(2*sqrt(2*log(2)));
fRate=5000;max(50000,5000*round(2/(pixsize))); % % frame rate, HZ
dt=1/fRate; 
nLineBase=500; 
nt=10000; % #of time points 2s
nTau=100;
nTauFit=nTau;
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
Vtot_1=Vtot; % total speed, m/s
AngleVes1= 0; % vessel 1 orientation, angle degree
WidthVes1= WidthVes;  % vessel width, in pixel, real vessel width=WidthVes*pixX
nPtcPerLineV1= 1/540*WidthVes;%WidthVes;%round(0.2*WidthVes);%WidthVes;%round(0.2*WidthVes);%0.2*WidthVes;%round(0.2*WidthVes); % number of RBCs per cross section for each vessel
vzVes1=Vtot_1*sin(AngleVes1/180*pi); % axial flow speed in Vessel 1, m/s
vxVes1=Vtot_1*cos(AngleVes1/180*pi); % Lateral-X flow speed in Vessel 1, m/s
V1z0=0; V1x0=0; V1y0=0;  % vessel 1 center location, um

% 2.2 parameters for vessel 2, 
Vtot_2=-5e-3; % total speed, m/s
AngleVes2=70; % vessel 1 orientation, angle degree
WidthVes2=20;  % vessel width, in pixel, real vessel width=WidthVes*pixX
nPtcPerLineV2=0; % number of RBCs per cross section for each vessel
vzVes2=Vtot_2*sin(AngleVes2/180*pi); % axial flow speed in Vessel 1, m/s
vxVes2=Vtot_2*cos(AngleVes2/180*pi); % Lateral-X flow speed in Vessel 1, m/s
V2z0=50; V2x0=0; V2y0=0;  % vessel 1 center location, um

%% 3. vessel location
% Vessel 1
Ves10=ones(WidthVes1,200); % vessel
if sign(Vtot_1)>0
    Ves10(:,1)=2; % source area to generate RBC
else
    Ves10(:,end)=2; % source area to generate RBC
end
Ves1=imrotate(Ves10,-AngleVes1,'nearest'); % rotate vessel 
%figure,imagesc(Ves1);axis equal tight
[nzRot,nxRot]=size(Ves1);
[I1S, J1S]=find(Ves1==2);
WidthVes1Rot=length(I1S);
V1SzP=((floor(nzRot/2)-I1S)*pixZ*1+V1z0)*1e-6; % rotated vessel source line position-Z, m
V1SxP=((J1S-floor(nxRot/2))*pixX*1+V1x0)*1e-6; % rotated vessel source line position-X, m
% figure,plot(V1SxP,V1SzP)

% Vessel 2
Ves20=ones(WidthVes2,200); % vessel
if sign(Vtot_2)>0
    Ves20(:,1)=2; % source area to generate RBC
else
    Ves20(:,end)=2; % source area to generate RBC
end
Ves2=imrotate(Ves20,-AngleVes2,'nearest'); % rotate vessel 
% figure,imagesc(Ves2);axis equal tight
[nzRot,nxRot]=size(Ves2);
[I2S, J2S]=find(Ves2==2);
WidthVes2Rot=length(I2S);
V2SzP=((floor(nzRot/2)-I2S)*pixZ+V2z0)*1e-6; % rotated vessel source line position-Z, m
V2SxP=((J2S-floor(nxRot/2))*pixX+V2x0)*1e-6; % rotated vessel source line position-X, m
% figure,plot(V2SxP,V2SzP)


% 3.4 generate RBC initial position
V1zP=zeros(nLineBase*nPtcPerLineV1,1);
V1xP=zeros(nLineBase*nPtcPerLineV1,1);

V2zP=zeros(nLineBase*nPtcPerLineV2,1);
V2xP=zeros(nLineBase*nPtcPerLineV2,1);
for iLine=1:nLineBase
    % randomly generate RBC position in the source area
    Temp1=rand(WidthVes1Rot,1);
    dTemp1=sort(Temp1(:),'descend');
    ind1=find(Temp1>dTemp1(nPtcPerLineV1+1)); 
    V1NewPtcl_zP{iLine}=V1SzP(ind1); % new RBC position-Z, m
    V1NewPtcl_xP{iLine}=V1SxP(ind1); % new RBC position-X, m    
    
    Temp2=rand(WidthVes2Rot,1);
    dTemp2=sort(Temp2(:),'descend');
    ind2=find(Temp2>dTemp2(nPtcPerLineV2+1)); 
    V2NewPtcl_zP{iLine}=V2SzP(ind2); % new RBC position-Z, m
    V2NewPtcl_xP{iLine}=V2SxP(ind2); % new RBC position-X, m    
  
    if iLine==1
        V1Ptcl_zP(:,iLine)=[V1NewPtcl_zP{iLine};V1zP]; % particle Z position for vessel 1
        V1Ptcl_xP(:,iLine)=[V1NewPtcl_xP{iLine};V1xP]; % particle X position for vessel 1
        
        V2Ptcl_zP(:,iLine)=[V2NewPtcl_zP{iLine};V2zP]; % particle Z position for vessel 2
        V2Ptcl_xP(:,iLine)=[V2NewPtcl_xP{iLine};V2xP]; % particle X position for vessel 2
    else
        V1Ptcl_zP(:,iLine)=circshift(V1Ptcl_zP(:,iLine-1),nPtcPerLineV1);
        V1Ptcl_zP(:,iLine)=[V1NewPtcl_zP{iLine};V1Ptcl_zP(nPtcPerLineV1+1:end,iLine)+sign(-vzVes1)*2*1e-6*sin(AngleVes1/180*pi)];
        V1Ptcl_xP(:,iLine)=circshift(V1Ptcl_xP(:,iLine-1),nPtcPerLineV1);
        V1Ptcl_xP(:,iLine)=[V1NewPtcl_xP{iLine};V1Ptcl_xP(nPtcPerLineV1+1:end,iLine)+sign(vxVes1)*2*1e-6*cos(AngleVes1/180*pi)];
              
        V2Ptcl_zP(:,iLine)=circshift(V2Ptcl_zP(:,iLine-1),nPtcPerLineV2);
        V2Ptcl_zP(:,iLine)=[V2NewPtcl_zP{iLine};V2Ptcl_zP(nPtcPerLineV2+1:end,iLine)+sign(-vzVes2)*5e-6*abs(sin(AngleVes2/180*pi))];
        V2Ptcl_xP(:,iLine)=circshift(V2Ptcl_xP(:,iLine-1),nPtcPerLineV2);
        V2Ptcl_xP(:,iLine)=[V2NewPtcl_xP{iLine};V2Ptcl_xP(nPtcPerLineV2+1:end,iLine)+sign(vxVes2)*5e-6*abs(cos(AngleVes2/180*pi))];
               
%         V2Ptcl_zP(:,iLine)=circshift(V2Ptcl_zP(:,iLine-1),nPtcPerLineV2);
%         V2Ptcl_zP(:,iLine)=[V2Ptcl_zP(1:end-nPtcPerLineV2,iLine)+sign(-vzVes2)*5e-6*sin(AngleVes2/180*pi);V2NewPtcl_zP{iLine}];
%         V2Ptcl_xP(:,iLine)=circshift(V2Ptcl_xP(:,iLine-1),nPtcPerLineV2);
%         V2Ptcl_xP(:,iLine)=[V2Ptcl_xP(1:end-nPtcPerLineV2,iLine)+sign(vxVes2)*5e-6*cos(AngleVes2/180*pi);V2NewPtcl_xP{iLine}];
    end
end
V1zP=V1Ptcl_zP(1:end-nPtcPerLineV1,nLineBase);
V1xP=V1Ptcl_xP(1:end-nPtcPerLineV1,nLineBase);
%figure,plot(V1xP*1e6,V1zP*1e6,'.'); axis equal tight; grid on
%  xlim([-200 200]); ylim([-200 200]);
 xlabel('x [um]'); ylabel('z [um]');
clear V1Ptcl_xP V1Ptcl_zP 

V2zP=V2Ptcl_zP(1:end-nPtcPerLineV2,nLineBase);
V2xP=V2Ptcl_xP(1:end-nPtcPerLineV2,nLineBase);
% figure,plot(V2xP,V2zP,'.'); axis equal tight; grid on
clear V2Ptcl_xP V2Ptcl_zP 
%% 4. time series signal 
% 4.1 RBC flow speed
IQ=zeros(1,1, nt);
iS1=1; % source index for vesel 1
iS2=1; % source index for vesel 2
for it=1:nt
    %% generate new RBC every 5 um for Ves1 
    if floor(it*dt*abs(Vtot_1)*1e6/pixZ)>=iS1
        Temp1=rand(WidthVes1Rot,1);
        dTemp1=sort(Temp1(:),'descend');
        ind1=find(Temp1>dTemp1(nPtcPerLineV1+1));
        V1NewPtcl_zP{it}=V1SzP(ind1);
        V1NewPtcl_xP{it}=V1SxP(ind1);
        iS1=floor(it*dt*abs(Vtot_1)*1e6/pixZ)+1;
        %% append new RBC position to existing position matrix
        if it==1
            V1Ptcl_zP(:,it)=[V1NewPtcl_zP{it};V1zP];
            V1Ptcl_xP(:,it)=[V1NewPtcl_xP{it};V1xP];
        else
            V1Ptcl_zP(:,it)=circshift(V1Ptcl_zP(:,it-1),nPtcPerLineV1);
            V1Ptcl_zP(:,it)=[V1NewPtcl_zP{it};V1Ptcl_zP(nPtcPerLineV1+1:end,it)-vzVes1*dt];
            V1Ptcl_xP(:,it)=circshift(V1Ptcl_xP(:,it-1),nPtcPerLineV1);
            V1Ptcl_xP(:,it)=[V1NewPtcl_xP{it};V1Ptcl_xP(nPtcPerLineV1+1:end,it)+vxVes1*dt];
        end
    else
        if it==1
            V1Ptcl_zP(:,it)=V1zP;
            V1Ptcl_xP(:,it)=V1xP;
        else
            V1Ptcl_zP(:,it)=V1Ptcl_zP(1:end,it-1)-vzVes1*dt;
            V1Ptcl_xP(:,it)=V1Ptcl_xP(1:end,it-1)+vxVes1*dt;
        end
    end
    
    if floor(abs(it*dt*Vtot_2*1e6/pixZ))>=iS2
        Temp2=rand(WidthVes2Rot,1);
        dTemp2=sort(Temp2(:),'descend');
        ind2=find(Temp2>dTemp2(nPtcPerLineV2+1));
        V2NewPtcl_zP{it}=V2SzP(ind2);
        V2NewPtcl_xP{it}=V2SxP(ind2);
        iS2=floor(it*dt*Vtot_2*1e6/pixZ)+1;
        %% append new RBC position to existing position matrix
        if it==1
            V2Ptcl_zP(:,it)=[V2NewPtcl_zP{it};V2zP];
            V2Ptcl_xP(:,it)=[V2NewPtcl_xP{it};V2xP];
        else
            
            V2Ptcl_zP(:,it)=circshift(V2Ptcl_zP(:,it-1),nPtcPerLineV2);
            V2Ptcl_zP(:,it)=[V2NewPtcl_zP{it};V2Ptcl_zP(nPtcPerLineV2+1:end,it)-vzVes2*dt];
            V2Ptcl_xP(:,it)=circshift(V2Ptcl_xP(:,it-1),nPtcPerLineV2);
            V2Ptcl_xP(:,it)=[V2NewPtcl_xP{it};V2Ptcl_xP(nPtcPerLineV2+1:end,it)+vxVes2*dt];
            
        end
    else
        if it==1
            V2Ptcl_zP(:,it)=V2zP;
            V2Ptcl_xP(:,it)=V2xP;
        else
            V2Ptcl_zP(:,it)=V2Ptcl_zP(1:end,it-1)-vzVes2*dt;
            V2Ptcl_xP(:,it)=V2Ptcl_xP(1:end,it-1)+vxVes2*dt;
        end
    end
    
%     hold on,plot(V1Ptcl_xP(:,it),V1Ptcl_zP(:,it),'.'); axis equal tight; grid on
    %% IQ data for the resolution voxel/pixel
%     IQ(1,1,it)=sum((exp(-(V1Ptcl_xP(:,it)-x0).^2/(2*Sigma(1)^2)).*exp(-(V1Ptcl_zP(:,it)-z0).^2/(2*Sigma(3)^2)).*exp(2i*k0*(V1Ptcl_zP(:,it)-z0))));%...
    data = exp(-(V1Ptcl_xP(:,it)-x0).^2/(2*Sigma(1)^2)).*exp(-(V1Ptcl_zP(:,it)-z0).^2/(2*Sigma(3)^2)).*exp(2i*k0*(V1Ptcl_zP(:,it)-z0)); 
     IQ(1,1,it)=sum((data(:,:)));%...1:round(end/2)
%     figure; plot(abs(data));
    data0(:,it) = data;
    nparticle0(:,it) = abs(data)>1e-4;%mean(abs(data(end-500:end)));   
    %+exp(-(V2Ptcl_xP(:,it)-x0).^2/(2*Sigma(1)^2)).*exp(-(V2Ptcl_zP(:,it)-z0).^2/(2*Sigma(3)^2)).*exp(2i*k0*(V2Ptcl_zP(:,it)-z0)));
end
nparticle = mean(sum(nparticle0,1),2);
% IQ = IQ.*conj(IQ); xlim = 200;
[vx,vz] = meshgrid([-200*1e-6:pixsize*1e-6:200*1e-6]+x0,[-200*1e-6:pixsize*1e-6:200*1e-6]+z0); 
val0 = exp(-(vx-x0).^2/(2*Sigma(1)^2)).*exp(-(vz-z0).^2/(2*Sigma(3)^2)).*exp(2i*k0*(vz-z0));    
% figure; imagesc([-xlim:2:xlim],[-xlim:2:xlim],abs(val0));axis image; xlabel('x[um]');ylabel('z[um]');colorbar;
% figure; imagesc([-xlim:2:xlim],[-xlim:2:xlim],real(val0));axis image; xlabel('x[um]');ylabel('z[um]');colorbar;
% figure; imagesc([-xlim:2:xlim],[-xlim:2:xlim],imag(val0));axis image; xlabel('x[um]');ylabel('z[um]');colorbar;
% figure; imagesc([-xlim:2:xlim],[-xlim:2:xlim],angle(val0));axis image; xlabel('x[um]');ylabel('z[um]');colorbar;
% % % figure;subplot(141); plot([-xlim:pixsize:xlim],abs(val0(:,xlim/pixsize)));subplot(142); plot([-xlim:pixsize:xlim],real(val0(:,xlim/pixsize)));subplot(143); plot([-xlim:pixsize:xlim],imag(val0(:,xlim/pixsize)));subplot(144); plot([-xlim:pixsize:xlim],angle(val0(:,xlim/pixsize)));
% figure;subplot(141); plot([-xlim:2:xlim],abs(val0(xlim/2,:)));subplot(142); plot([-xlim:2:xlim],real(val0(xlim/2,:)));subplot(143); plot([-xlim:2:xlim],imag(val0(xlim/2,:)));subplot(144); plot([-xlim:2:xlim],angle(val0(xlim/2,:)));
% % 5. plot and save particle flow animation
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
snr_real = 20 .* log10(Signal_rms ./ Noise_rms);
% print({'SNR is', num2str(snr_real),'dB'})
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
%% plot g1
fig=figure;
set(fig,'Position',[400 400 800 650])
subplot(2,2,1);hold on, plot(tauCoor,real(squeeze(g1)));title('real(g1), RAW'); xlabel('time lag, [ms]')
subplot(2,2,2);hold on, plot(tauCoor,imag(squeeze(g1)));title('imag(g1), RAW');xlabel('time lag, [ms]')
subplot(2,2,3);hold on, plot(squeeze(g1).');title('cplx(g1), RAW');%ylim([-1,1]); xlim([-1,1])
subplot(2,2,4);hold on, plot(tauCoor,abs(squeeze(g1)));title('abs(g1), RAW'); 
xlabel('time lag, [ms]')
end
