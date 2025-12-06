%% DLS-OCT process, GPU-based
% input: % GG: [nVox,nTau], nVox=nz*nx*ny
% PRSinfo: processing information
% PRSinfo.FWHM: (transverse, axial), m
% PRSinfo.fAline: DAQ Aline rate, Hz
% PRSinfo.Lam: [light source center, wavelength bandwidth], m
% output:
    % Ms, Mf, R
    % Vt, mm/s
    % Vz, mm/s
    % D, um^2/s
function [Ms, Mf, Vt, Vz, D, R, GGf]=GG2VDR_NoFit_GPU(GG, PRSinfo)
%% constant %%%%%%%%%%%%
Sigma=PRSinfo.FWHM*0.7/(2*sqrt(2*log(2))); % intensity-based sigma
Sigma2=2*Sigma;
dt = 1/PRSinfo.fAline; 
k0 = 2*pi/PRSinfo.Lam(1); % /m
% dk = ( 2*pi/(PRSinfo.Lam(1)-PRSinfo.Lam(2)/2)-2*pi/(PRSinfo.Lam(1)+PRSinfo.Lam(2)/2) )/2*sqrt(2*log(2));
n = 1.35;  q = 2*n*k0; 
[nVox,nTau]=size(GG); % nVox=nz*nx*ny
nz=PRSinfo.Dim(1); nx=PRSinfo.Dim(2); ny=PRSinfo.Dim(3); 
tau = [1:nTau]*dt; % time lag, s
t = tau; tn = t / tau(end);
if nz*nx*ny>100*400*5
    if rem(ny,5)==0
        nyPchk=5;
    else
        nyPchk=1;
    end
    nyChk=ny/nyPchk;
else
    nyPchk=ny;
    nyChk=ny/nyPchk;
end
nVoxPchk=nz*nx*nyPchk;
PRSinfo.Dim=[nz,nx,nyPchk,nTau];
for iChk=1:nyChk
    iVoxStart=(iChk-1)*nVoxPchk+1;
    iVoxEnd=iChk*nVoxPchk;
    iGG=(GG(iVoxStart:iVoxEnd,:));
    % [vz0]=GG2Vz(GG, PRSinfo, 10); % m/s
    [vz0]=GG2Vz_GPU(iGG, PRSinfo, 10); % m/s
    iMs0 = min(max(real(FindCOR(iGG)),0),max(mean(real(iGG(:,floor(end*2/3):end)),2),0));
    iMe0=1-abs(iGG(:,1));
    iMf0=max(1-iMs0-iMe0,0);
    [Vt0(iVoxStart:iVoxEnd),Vz0(iVoxStart:iVoxEnd),D0(iVoxStart:iVoxEnd),R0(iVoxStart:iVoxEnd)]=iniDLSOCT_GPU(iGG, vz0, iMs0, iMf0, PRSinfo);
    Ms0(iVoxStart:iVoxEnd)=iMs0;
    Mf0(iVoxStart:iVoxEnd)=iMf0;
end
Ms=reshape(Ms0,[nz,nx,ny]);
Mf=reshape(Mf0,[nz,nx,ny]);
Vt=reshape(Vt0,[nz,nx,ny])*1e3; % mm/s
Vz=reshape(Vz0,[nz,nx,ny])*1e3;  % mm/s
D=reshape(D0,[nz,nx,ny])*1e12; % um^2/s
R=reshape(R0,[nz,nx,ny]);
% %% 1, determine the initial guess of vz0, Ms0, Me0, and Mf0
% % [vz0]=GG2Vz(GG, PRSinfo, 10); % m/s
% [vz0]=GG2Vz_GPU(GG, PRSinfo, 10); % m/s
% Ms0 = min(max(real(FindCOR(GG)),0),max(mean(real(GG(:,floor(end*2/3):end)),2),0));
% Me0=1-abs(GG(:,1));
% Mf0=max(1-Ms0-Me0,0);
% CR=1;%(Mf0>0.05);
% % [Vt0,Vz0,D0,R0]=iniDLSOCT(GG, vz0, Ms0, Mf0, PRSinfo);
% [Vt0,Vz0,D0,R0]=iniDLSOCT_GPU(GG, vz0, Ms0, Mf0, PRSinfo);
% % 2. Fitting constraint
% Fmin_cstrn(:,:,1)=[Ms0-0.05 Ms0+0.05];   % Ms constrain
% Fmin_cstrn(:,:,2)=[max(Mf0-0.00, 0) min(Mf0+0.1,1)];   % MfR constrain
% Fmin_cstrn(:,:,3)=[Vt0-4e-3 Vt0+4e-3]*tau(end)/(Sigma2(1)); % Vt constrain
% Fmin_cstrn(:,:,4)=sign(Vz0).*[0.5*abs(Vz0) 1.3*abs(Vz0)]*tau(end);  % Vz constrain
% Fmin_cstrn(:,:,5)=[0.9*D0 min(D0,100*1e-12)]*q^2*tau(end); % D constrain
% %% 3. non-linear least square fitting
% fitC0(:,:,1) = double(Ms0); % initials
% fitC0(:,:,2) = double(Mf0); % initials
% fitC0(:,:,3) = double(Vt0*tau(end)/(Sigma2(1))); % initials
% fitC0(:,:,4) = double(Vz0*tau(end)); % initials
% fitC0(:,:,5) = double(D0*q^2*tau(end)); % initials
% Tn=(tn);
% fitC0=gather(fitC0);
% Fmin_cstrn=(double(gather(Fmin_cstrn)));
% warning('off');
% fit = @(c) sum( abs(c(:,1,1)+ c(:,1,2).*exp( -(c(:,1,3).*Tn).^2-c(:,1,4).^2/(Sigma2(2).^2).*Tn.^2 -c(:,1,5).*Tn ).*exp(1i*q*c(:,1,4).*Tn) - GG ).^2 ,2);
% [fitC, fval] = fmincon(fit, fitC0, [],[],[],[], ...
%     [Fmin_cstrn(:,1,1) Fmin_cstrn(:,1,2) Fmin_cstrn(:,1,3) Fmin_cstrn(:,1,4) Fmin_cstrn(:,1,5)], ...
%     [Fmin_cstrn(:,2,1) Fmin_cstrn(:,2,2) Fmin_cstrn(:,2,3) Fmin_cstrn(:,2,4) Fmin_cstrn(:,2,5)], ...
%     [], optimset('Display','off','TolFun',1e-6,'TolX',1e-6));%
% Ms=reshape(fitC(:,1,1),[nz,nx,ny]); 
% Mf=reshape(fitC(:,1,2),[nz,nx,ny]);  
% Vt=reshape(fitC(:,1,3).*CR,[nz,nx,ny])/(tau(end)/(Sigma2(1)))*1e3; % mm/s
% Vz=reshape(fitC(:,1,4).*CR,[nz,nx,ny])/tau(end)*1e3;  % mm/s
% D=reshape(fitC(:,1,5).*CR,[nz,nx,ny])/(q^2*tau(end))*1e12; % um^2/s
% GGf=fitC(:,:,1)+fitC(:,:,2).*exp( -(fitC(:,:,3).*Tn).^2-fitC(:,:,4).^2/(Sigma2(2).^2).*Tn.^2 -fitC(:,:,5).*Tn).*exp(1i*q*fitC(:,:,4).*Tn);
% R=gather(reshape((1-sum(abs(GG-GGf).^2,2)./sum(abs(abs(GG)-mean(GG,2)).^2,2)).*CR,[nz,nx,ny]));  
% %% 3. non-linear least square fitting, grid-based
% [C1, ~]=meshgrid(Ms0, tn);
% [C2, ~]=meshgrid(Mf0, tn);
% [C3, ~]=meshgrid(Vt0*tau(end)/(Sigma2(1)), tn);
% [C4, ~]=meshgrid(Vz0*tau(end), tn);
% [C5, Tn]=meshgrid(D0*q^2*tau(end), tn);
% fitC0(:,:,1) = double(C1.'); % initials
% fitC0(:,:,2) = double(C2.'); % initials
% fitC0(:,:,3) = double(C3.'); % initials
% fitC0(:,:,4) = double(C4.'); % initials
% fitC0(:,:,5) = double(C5.'); % initials
% Tn=gpuArray(Tn.');
% fitC0=gather(fitC0);
% Fmin_cstrn=(double(gather(Fmin_cstrn)));
% warning('off');
% fit = @(c) sum( abs(c(:,:,1)+ c(:,:,2).*exp( -(c(:,:,3).*Tn).^2-c(:,:,4).^2/(Sigma2(2).^2).*Tn.^2 -c(:,:,5).*Tn ).*exp(1i*q*c(:,:,4).*Tn) - GG ).^2 ,2);
% [fitC, fval] = fmincon(fit, fitC0, [],[],[],[], ...
%     [Fmin_cstrn(:,1,1) Fmin_cstrn(:,1,2) Fmin_cstrn(:,1,3) Fmin_cstrn(:,1,4) Fmin_cstrn(:,1,5)], ...
%     [Fmin_cstrn(:,2,1) Fmin_cstrn(:,2,2) Fmin_cstrn(:,2,3) Fmin_cstrn(:,2,4) Fmin_cstrn(:,2,5)], ...
%     [], optimset('Display','off','TolFun',1e-6,'TolX',1e-6));%
% Ms=reshape(fitC(:,1,1),[nz,nx,ny]); 
% Mf=reshape(fitC(:,1,2),[nz,nx,ny]);  
% Vt=reshape(fitC(:,1,3).*CR,[nz,nx,ny])/(tau(end)/(Sigma2(1)))*1e3; % mm/s
% Vz=reshape(fitC(:,1,4).*CR,[nz,nx,ny])/tau(end)*1e3;  % mm/s
% D=reshape(fitC(:,1,5).*CR,[nz,nx,ny])/(q^2*tau(end))*1e12; % um^2/s
% GGf=fitC(:,:,1)+fitC(:,:,2).*exp( -(fitC(:,:,3).*Tn).^2-fitC(:,:,4).^2/(Sigma2(2).^2).*Tn.^2 -fitC(:,:,5).*Tn).*exp(1i*q*fitC(:,:,4).*Tn);
% R=gather(reshape(1-sum(abs(GG-GGf).^2,2)./sum(abs(abs(GG)-mean(GG,2)).^2,2),[nz,nx,ny]));    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% no non-linear least squared fitting 
% %% 1, determine the initial guess of vz0, Ms0, Me0, and Mf0
