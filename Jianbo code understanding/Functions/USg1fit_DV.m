%% fit g1, complex, in vivo, Gaussian distributied velocity model
function [Ms, Mf, Vx, Vy, Vz,Vz0,Pvz,R,GGf]=USg1fit_DV(GG, A_FWHM, rFrame)
% Sigma_vz=Pvz*Vz
% GG(nz,nx,nTau)
% Res: X, Y, Z spatial resolution, FWHM
% rFrame: frame rate, Hz
%% I. constant
[nz,nx,nTau]=size(GG);
f0=16.625e6;            % ultrasound transducer center frequency, Hz
C=1540;             % ultrasound speed in biological tissue, m/s
lambda0=C/f0;        % wavlength
k0 = 2*pi/lambda0;   % wave number
Sigma=A_FWHM*0.7/(2*sqrt(2*log(2))); % intensity-based sigma
Sigma2=2*Sigma;
nItpVz0=10;          % for Vz0 determination
dt = 1/rFrame;      % frame interval, s
tau = [1:nTau]*dt; % time lag, s
t = tau.';
tn = t / tau(end);
%% III. fitting procedure
Vz00=zeros(nz,nx); Vx=Vz00; Vy=Vz00; Vz=Vz00; GGf=zeros(nz,nx,nTau);
warning off;

for iz=1:nz
    for ix=1:nx
        gg=squeeze(GG(iz,ix,:));
        if abs(gg(1))>0.55 && abs(gg(2))<0.25 && abs(gg(2))<abs(gg(3))
            gg(1)=gg(2)+(abs(real(gg(2)-gg(3)))+1i*abs(imag(gg(2)-gg(3))))*1.5;
             %% III.1 determine Vz0
            [Vz00]=DtmVz0(gg, rFrame, lambda0, nItpVz0); % m/s
            Ms0=min(max(real(US_FindCOR(permute(gg(floor(end*1/2):end),[3,2,1]))),0),max(mean(real(gg(floor(end*2/3):end))),0));
            Me0 =1-abs(gg(1));
            MfR0 = max(1-Ms0-Me0,0);
            %% III.2 determine Vt0, PVz0
            [vx0,vz0,PVz0,MfI0, R0,GGf0]=DtmVx0nPVz0(gg, Vz00, Ms0, MfR0, k0,A_FWHM, rFrame);
            Vz00=Vz00*1e3;
            vz0=vz0*1e3;
            Vz0(iz,ix,:)=[Vz00, vz0];
            
            Ms(iz,ix)=0;
            R(iz,ix) = R0;
            Mf(iz,ix,:)=[MfI0 MfR0];
            GGf(iz,ix,:)=GGf0;
            Vx(iz,ix)=abs(vz0)*MfR0; Vy(iz,ix)=0; Vz(iz,ix)=vz0; Pvz(iz,ix)=PVz0; 
        else
            if abs(abs(gg(1))-abs(gg(2)))>2*abs((abs(gg(2))-abs(gg(3)))) % && abs(gg(2))<abs(gg(3))
                gg(1)=gg(2)+(abs(real(gg(2)-gg(3)))+1i*abs(imag(gg(2)-gg(3))))*1.5;
            end
            %% III.1 determine Vz0
            [Vz00]=DtmVz0(gg, rFrame, lambda0, nItpVz0); % m/s
            Ms0=min(max(real(US_FindCOR(permute(gg(floor(end*1/2):end),[3,2,1]))),0),max(mean(real(gg(floor(end*2/3):end))),0));
            Me0 =1-abs(gg(1));
            MfR0 = max(1-Ms0-Me0,0);
            %% III.2 determine Vt0, PVz0
            [vx0,vz0,PVz0,MfI0, R0]=DtmVx0nPVz0(gg, Vz00, Ms0, MfR0, k0,A_FWHM, rFrame);
            Vz0(iz,ix,:)=[Vz00, vz0]*1e3;
            %% III.3 fitting constrain
            Vy0=0;
            if R0<=0
                R0=0.3;
            end
            Fmin_cstrn(:,1)=[Ms0-0.1 Ms0+0.05];   % Ms constrain
            Fmin_cstrn(:,2)=[max(MfR0-0.00, 0) min(MfR0+0.05,1)];   % MfR constrain
            Fmin_cstrn(:,3)=[MfI0-0.2 MfI0+0.2];   % MfI constrain
            Fmin_cstrn(:,4)=(vx0+[-10*min(max(0.5*(0.9/R0)^4,0),1),min(5*max(0.5*(0.9/R0)^4,0),0.5*abs(vx0*1e3))]*1e-3)*tau(end)/(Sigma2(1));  % Vx constrain
            Fmin_cstrn(:,5)=[0 0]*1e-3*tau(end)/(Sigma2(2));  % Vy constrain
            %                 Fmin_cstrn(:,5)=(vx0+[-10*min(max(0.5*(0.9/R1)^4,0),1),2*min(max(0.5*(0.9/R1)^4,0),1)]*1e-3)*tau(end)/(Sigma2(2));  % Vy constrain
            Fmin_cstrn(:,6)=sign(vz0)*(abs(vz0)+[-1*abs(vz0) max(0.5*abs(vz0),1e-3)]);  % Vz constrain
            Fmin_cstrn(:,7)=PVz0+[-0.1 0.1];  % Pv constrain
            %         Fmin_cstrn(:,7)=[0.5 0.5]; % temp, constrain Pv to 0.5
            %% III.4 fit complex (g1)
            fitE = @(c) sum( abs(c(1) + c(2).*exp( -(c(4)*tn).^2-(c(5)*tn).^2-(c(6)*t).^2/(Sigma2(3))^2).*exp(-(k0*t*c(7)*c(6)).^2).*cos(2*k0*c(6)*t)+...
                1i*c(3).*exp( -(c(4)*tn).^2-(c(5)*tn).^2-(c(6)*t).^2/(Sigma2(3))^2).*exp(-(k0*t*c(7)*c(6)).^2).*sin(2*k0*c(6)*t)- gg ).^2 );
            fitC0 = [Ms0,MfR0,MfI0, vx0*tau(end)/(Sigma2(1)), Vy0*tau(end)/(Sigma2(2)), vz0,PVz0];
            [fitC, fval] = fmincon(fitE, fitC0, [],[],[],[], ...
                [Fmin_cstrn(1,1) Fmin_cstrn(1,2) Fmin_cstrn(1,3) Fmin_cstrn(1,4) Fmin_cstrn(1,5) Fmin_cstrn(1,6) Fmin_cstrn(1,7)], ...
                [Fmin_cstrn(2,1) Fmin_cstrn(2,2) Fmin_cstrn(2,3) Fmin_cstrn(2,4) Fmin_cstrn(2,5) Fmin_cstrn(2,6) Fmin_cstrn(2,7)], ...
                [], optimset('Display','off','TolFun',1e-6,'TolX',1e-6));%
            ms=fitC(1);
            mfR=fitC(2);
            mfI=fitC(3);
            vx = abs(fitC(4))*Sigma2(1)/tau(end);
            vy = abs(fitC(5))*Sigma2(2)/tau(end);
            vz = sign(vz0)*abs(fitC(6));
            Pv = fitC(7);
            ggFit=ms + mfR*exp(-vx^2*t.^2/(Sigma2(1))^2-vy^2*t.^2/(Sigma2(2))^2-vz^2*t.^2/(Sigma2(3))^2).*exp(-(k0*t*Pv*vz).^2).*cos(2*k0*vz*t)...
                +1i*mfI*exp(-vx^2*t.^2/(Sigma2(1))^2-vy^2*t.^2/(Sigma2(2))^2-vz^2*t.^2/(Sigma2(3))^2).*exp(-(k0*t*Pv*vz).^2).*sin(2*k0*vz*t);
            R(iz,ix) = 1 - sum( abs((gg)-ggFit).^2 ) / sum( abs((gg)-mean((gg))).^2);
            Ms(iz,ix)=ms;
            Mf(iz,ix,:)=[mfI mfR];
            GGf(iz,ix,:)=ggFit;
            Vx(iz,ix)=vx*1e3; Vy(iz,ix)=vy*1e3; Vz(iz,ix)=vz*1e3; Pvz(iz,ix)=Pv;
        end
    end
end


% %% fit g1, complex, in vivo, Gaussian distributied velocity model
% function [Ms, Mf, Vx, Vy, Vz,Pvz,R,GGf]=USg1fit_DV(GG, A_FWHM, rFrame)
% % Sigma_vz=Pvz*Vz
% % GG(nz,nx,nTau)
% % Res: X, Y, Z spatial resolution, FWHM
% % rFrame: frame rate, Hz
% %% I. constant
% [nz,nx,nTau]=size(GG);
% f0=16.625e6;            % ultrasound transducer center frequency, Hz
% C=1540;             % ultrasound speed in biological tissue, m/s
% lambda0=C/f0;        % wavlength
% k0 = 2*pi/lambda0;   % wave number
% Sigma=A_FWHM*0.7/(2*sqrt(2*log(2))); % intensity-based sigma
% Sigma2=2*Sigma;
% nItpVz0=10;          % for Vz0 determination
% dt = 1/rFrame;      % frame interval, s
% tau = [1:nTau]*dt; % time lag, s
% t = tau.';
% tn = t / tau(end);
% %% III. fitting procedure
% Vz0=zeros(nz,nx); Vx=Vz0; Vy=Vz0; Vz=Vz0; GGf=zeros(nz,nx,nTau);
% warning off;
% 
% for iz=1:nz
%     for ix=1:nx
%         gg=squeeze(GG(iz,ix,:));
%         if abs(gg(1))>0.55 && abs(gg(2))<0.25 && abs(gg(2))<abs(gg(3))
%             Ms(iz,ix)=0;
%             mfR=0;
%             mfI=0;
%             R(iz,ix) = 0;
%             Mf(iz,ix,:)=[mfI mfR];
%             GGf(iz,ix,:)=gg;
%             Vx(iz,ix)=0; Vy(iz,ix)=0; Vz(iz,ix)=0; Pvz(iz,ix)=0;
%         else
%             if abs(abs(gg(1))-abs(gg(2)))>2*abs((abs(gg(2))-abs(gg(3)))) % && abs(gg(2))<abs(gg(3))
%                 gg(1)=gg(2)+(abs(real(gg(2)-gg(3)))+1i*abs(imag(gg(2)-gg(3))))*1.2;
%             end
%             %% III.1 determine Vz0
%             [Vz00]=DtmVz0(gg, rFrame, lambda0, nItpVz0); % m/s
%             Ms0=min(max(real(US_FindCOR(permute(gg(floor(end*1/2):end),[3,2,1]))),0),max(mean(real(gg(floor(end*2/3):end))),0));
%             Me0 =1-abs(gg(1));
%             MfR0 = max(1-Ms0-Me0,0);
%             %% III.2 determine Vt0, PVz0
%             [vx0,vz0,PVz0,MfI0, R0]=DtmVx0nPVz0(gg, Vz00, Ms0, MfR0, k0,A_FWHM, rFrame);
%             %% III.3 fitting constrain
%             Vy0=0;
%             if R0<=0
%                 R0=0.3;
%             end
%             Fmin_cstrn(:,1)=[Ms0-0.1 Ms0+0.05];   % Ms constrain
%             Fmin_cstrn(:,2)=[max(MfR0-0.00, 0) min(MfR0+0.05,1)];   % MfR constrain
%             Fmin_cstrn(:,3)=[MfI0-0.2 MfI0+0.2];   % MfI constrain
%             Fmin_cstrn(:,4)=(vx0+[-10*min(max(0.5*(0.9/R0)^4,0),1),min(5*max(0.5*(0.9/R0)^4,0),0.5*abs(vx0*1e3))]*1e-3)*tau(end)/(Sigma2(1));  % Vx constrain
%             Fmin_cstrn(:,5)=[0 0]*1e-3*tau(end)/(Sigma2(2));  % Vy constrain
%             %                 Fmin_cstrn(:,5)=(vx0+[-10*min(max(0.5*(0.9/R1)^4,0),1),2*min(max(0.5*(0.9/R1)^4,0),1)]*1e-3)*tau(end)/(Sigma2(2));  % Vy constrain
%             Fmin_cstrn(:,6)=sign(vz0)*(abs(vz0)+[-1*abs(vz0) max(0.5*abs(vz0),1e-3)]);  % Vz constrain
%             Fmin_cstrn(:,7)=PVz0+[-0.1 0.1];  % Pv constrain
%             %         Fmin_cstrn(:,7)=[0.5 0.5]; % temp, constrain Pv to 0.5
%             %% III.4 fit complex (g1)
%             fitE = @(c) sum( abs(c(1) + c(2).*exp( -(c(4)*tn).^2-(c(5)*tn).^2-(c(6)*t).^2/(Sigma2(3))^2).*exp(-(k0*t*c(7)*c(6)).^2).*cos(2*k0*c(6)*t)+...
%                 1i*c(3).*exp( -(c(4)*tn).^2-(c(5)*tn).^2-(c(6)*t).^2/(Sigma2(3))^2).*exp(-(k0*t*c(7)*c(6)).^2).*sin(2*k0*c(6)*t)- gg ).^2 );
%             fitC0 = [Ms0,MfR0,MfI0, vx0*tau(end)/(Sigma2(1)), Vy0*tau(end)/(Sigma2(2)), vz0,PVz0];
%             [fitC, fval] = fmincon(fitE, fitC0, [],[],[],[], ...
%                 [Fmin_cstrn(1,1) Fmin_cstrn(1,2) Fmin_cstrn(1,3) Fmin_cstrn(1,4) Fmin_cstrn(1,5) Fmin_cstrn(1,6) Fmin_cstrn(1,7)], ...
%                 [Fmin_cstrn(2,1) Fmin_cstrn(2,2) Fmin_cstrn(2,3) Fmin_cstrn(2,4) Fmin_cstrn(2,5) Fmin_cstrn(2,6) Fmin_cstrn(2,7)], ...
%                 [], optimset('Display','off','TolFun',1e-6,'TolX',1e-6));%
%             ms=fitC(1);
%             mfR=fitC(2);
%             mfI=fitC(3);
%             vx = abs(fitC(4))*Sigma2(1)/tau(end);
%             vy = abs(fitC(5))*Sigma2(2)/tau(end);
%             vz = sign(vz0)*abs(fitC(6));
%             Pv = fitC(7);
%             ggFit=ms + mfR*exp(-vx^2*t.^2/(Sigma2(1))^2-vy^2*t.^2/(Sigma2(2))^2-vz^2*t.^2/(Sigma2(3))^2).*exp(-(k0*t*Pv*vz).^2).*cos(2*k0*vz*t)...
%                 +1i*mfI*exp(-vx^2*t.^2/(Sigma2(1))^2-vy^2*t.^2/(Sigma2(2))^2-vz^2*t.^2/(Sigma2(3))^2).*exp(-(k0*t*Pv*vz).^2).*sin(2*k0*vz*t);
%             R(iz,ix) = 1 - sum( abs((gg)-ggFit).^2 ) / sum( abs((gg)-mean((gg))).^2);
%             Ms(iz,ix)=ms;
%             Mf(iz,ix,:)=[mfI mfR];
%             GGf(iz,ix,:)=ggFit;
%             Vx(iz,ix)=vx*1e3; Vy(iz,ix)=vy*1e3; Vz(iz,ix)=vz*1e3; Pvz(iz,ix)=Pv;
%         end
%     end
% end

