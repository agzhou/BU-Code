% input: 
    % sIQ: bulk motion removed data
    % PRSSinfo: data acquistion information, including
        % PRSSinfo.FWHM: (X, Y, Z) spatial resolution, Full Width at Half Maximum of point spread function, m
        % PRSSinfo.rFrame: sIQ frame rate, Hz
        % PRSSinfo.f0: Transducer center frequency, Hz
        % PRSSinfo.C: Sound speed in the sample, m/s
        % PRSSinfo.g1nT: g1 calculation sample number
        % PRSSinfo.g1nTau: maximum number of time lag
        % PRSSinfo.SVDrank: SVD rank [low high]
        % PRSSinfo.HPfC:  High pass filtering cutoff frequency, Hz
        % PRSSinfo.NEQ: do noise equalization? 0: no noise equalization; 1: apply noise equalization
        % PRSSinfo.rfnScale: spatial refind scale
 % output:
    % Ms: static component fraction, [nz,nx]
    % Mf: dynamic component fraction, [nz,nx,2], 2: [real,imag]
    % Vx: x-direction velocity component, [nz,nx], mm/s
    % Vz: axial-direction velocity component, [nz,nx], mm/s
    % V=sqrt(Vx.^2+Vz.^2), [nz,nx], mm/s
    % pVz: Vz distribution (sigma-Vz), [nz,nx]
    % R: fitting accuracy, [nz,nx]
    % Vcz: Color Doppler axial velocity, [nz,nx], mm/s
    % GGf: gg fitting results, [nz,nx, nTau]
 % Bingxue Liu, 20220125
function [Vz,Vx,pVz,Ms,Mf,R,GGf, GGfid]=GG2vUS_fit(GG, g1Vz0, Ms0, MfR0, PRSSinfo)
%% O. constant
lambda0=PRSSinfo.C/PRSSinfo.f0;        % wavlength
k0 = 2*pi/lambda0;   % wave number
PRSSinfo.FWHM=[PRSSinfo.FWHM(1) 1000e-6 PRSSinfo.FWHM(2)]; % just to put the FWHM_y, any number other than 0
Sigma=PRSSinfo.FWHM*0.7/(2*sqrt(2*log(2))); % intensity-based sigma
Sigma2=2*Sigma;
dt = 1/PRSSinfo.rFrame;      % frame interval, s
tau = [1:PRSSinfo.g1nTau]*dt; % time lag, s
tn = tau / tau(end);
nz=PRSSinfo.Dim(1);nx = PRSSinfo.Dim(2); nTau = PRSSinfo.Dim(3);
%% IV. vUS initial
% [Vx0,Vz0,PVz0,MfI0,R0]=iniVx0Vz0Pvz0(GG, g1Vz0, Ms0, MfR0, PRSSinfo);
[Vz0, Vx0, PVz0, Ms0, Mf0, R0, GGf0]=GG2vUS(GG, reshape(g1Vz0,[nz*nx,1]), Ms0, MfR0, PRSSinfo);
MfI0 = reshape(Mf0(:,:,2),[nz*nx,1]);PVz0 = reshape(PVz0,[nz*nx,1]);
Vz0 = reshape(Vz0/1e3,[nz*nx,1]);Vx0 = reshape(Vx0/1e3,[nz*nx,1]);Ms0 = reshape(Ms0,[nz*nx,1]);
fitC0(:,:,1) = double(Ms0); % initials
fitC0(:,:,2) = double(MfR0); % initials
fitC0(:,:,3) = double(MfI0); % initials
fitC0(:,:,4) = double(Vx0); % initials *tau(end)/(Sigma2(1))
fitC0(:,:,5) = double(Vz0); % initials
fitC0(:,:,6) = double(PVz0); % initials
fitC0=double((fitC0));
%% V. vUS fitting
% V.1. Fitting constraint
Fmin_cstrn(:,:,1)=[Ms0-0.1 Ms0+0.1];   % Ms constrain
Fmin_cstrn(:,:,2)=[max(MfR0-0.15, 0) min(MfR0+0.1,1-Ms0)];   % MfR constrain
Fmin_cstrn(:,:,3)=[max(MfI0-0.15, 0) min(MfI0+0.2,1)];   % MfI constrain
Fmin_cstrn(:,:,4)=[0.5*Vx0 1.3*Vx0];%*tau(end)/(Sigma2(1)); % Vx constrain
Fmin_cstrn(:,:,5)=sign(Vz0).*[0.6*abs(Vz0) 1.2*abs(Vz0)];  % Vz constrain, mm/s
Fmin_cstrn(:,:,6)=[0.8*PVz0 min(PVz0*1.3,0.7)]; % PVz constrain
Fmin_cstrn0=double((Fmin_cstrn));
Index_cstrn = find(Fmin_cstrn0(:,1,:)>Fmin_cstrn0(:,2,:));
Fmin_cstrn_lb = Fmin_cstrn0(:,1,:);
Fmin_cstrn_ub = Fmin_cstrn0(:,2,:);
Fmin_cstrn_ub(Index_cstrn) = Fmin_cstrn_lb(Index_cstrn);
Fmin_cstrn(:,1,:) = Fmin_cstrn_lb;
Fmin_cstrn(:,2,:) = Fmin_cstrn_ub;
%% V.2 fit complex (g1)
fitE = @(c) double(sum(abs(c(:,1,1) + c(:,1,2).*exp( -(c(:,1,4).*tau).^2/(Sigma2(1))^2-(c(:,1,5).*tau).^2/(Sigma2(3))^2).*exp(-(k0*c(:,1,5).*c(:,1,6).*tau).^2).*cos(2*k0*c(:,1,5).*tau)+...
    1i*c(:,1,3).*exp( -(c(:,1,4).*tau).^2/(Sigma2(1))^2-(c(:,1,5).*tau).^2/(Sigma2(3))^2).*exp(-(k0*c(:,1,5).*c(:,1,6).*tau).^2).*sin(2*k0*c(:,1,5).*tau)- (GG) ).^2 ,2));
[fitC, fval] = fmincon(fitE, fitC0,[],[],[],[], ...
    [Fmin_cstrn(:,1,1) Fmin_cstrn(:,1,2) Fmin_cstrn(:,1,3) Fmin_cstrn(:,1,4) Fmin_cstrn(:,1,5) Fmin_cstrn(:,1,6)], ...
    [Fmin_cstrn(:,2,1) Fmin_cstrn(:,2,2) Fmin_cstrn(:,2,3) Fmin_cstrn(:,2,4) Fmin_cstrn(:,2,5) Fmin_cstrn(:,2,6)], ...
    [], optimset('Display','notify','TolFun',1e-6,'TolX',1e-6));%,'ScaleProblem','obj-and-constr'

Ms=reshape(fitC(:,1,1),[nz,nx]); 
Mf(:,:,1)=reshape(fitC(:,1,2),[nz,nx]);  % MfR
Mf(:,:,2)=reshape(fitC(:,1,3),[nz,nx]);  % MfI
Vx=reshape(fitC(:,1,4),[nz,nx]); % m/s  /(tau(end)/(Sigma2(1)))
Vz=reshape(fitC(:,1,5),[nz,nx]);  % m/s
V=sqrt(Vz.^2+Vx.^2);
pVz=reshape(fitC(:,1,6),[nz,nx]); % 
GGf0=fitC(:,:,1) + fitC(:,:,2).*exp( -(fitC(:,:,4).*tau).^2/(Sigma2(1))^2-(fitC(:,:,5).*tau).^2/(Sigma2(3))^2).*exp(-(k0*tau.*fitC(:,:,5).*fitC(:,:,6)).^2).*cos(2*k0*fitC(:,:,5).*tau)+...
    1i*fitC(:,:,3).*exp( -(fitC(:,:,4).*tau).^2/(Sigma2(1))^2-(fitC(:,:,5).*tau).^2/(Sigma2(3))^2).*exp(-(k0*tau.*fitC(:,:,5).*fitC(:,:,6)).^2).*sin(2*k0*fitC(:,:,5).*tau);
R=((reshape((1-sum(abs(GG-GGf0).^2,2)./sum(abs((GG)-mean(GG,2)).^2,2)),[nz,nx])));  
GGf=(reshape(GGf0,[nz,nx,nTau]));
GGf0id = exp( -(fitC(:,:,4).*tau).^2/(Sigma2(1))^2-(fitC(:,:,5).*tau).^2/(Sigma2(3))^2).*exp(-(k0*tau.*fitC(:,:,5).*fitC(:,:,6)).^2).*exp(2*1i*k0*fitC(:,:,5).*tau);
GGfid = reshape(GGf0id,[nz,nx,nTau]);
V=(V*1e3); Vz=(Vz*1e3); Vx=(Vx*1e3) ;% mm/s          
