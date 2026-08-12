%% calculate ACF (GG) of selected point (z,x)
function [Ms, Mf, Vx, Vy, Vz,R,GGf]=US_g1fit(GG, A_FWHM, rFrame)
% GG(nz,nx,nTau)
% Res: X, Y, Z spatial resolution
% rFrame: frame rate, Hz
%% 0. constant
[nz,nx,nTau]=size(GG);
f0=18e6;%15.625e6; % ultrasound transducer center frequency, Hz
C=1540;  % ultrasound speed in biological tissue, m/s
k0 = 2*pi*f0/C; % wave number
q=2*k0;
dt = 1/rFrame;  % frame interval, s
tau0 = [1:nTau]*dt; % time lag, s
Sigma=A_FWHM*0.7/(2*sqrt(2*log(2)));
Sigma2=2*Sigma;
% Sigma=R_FWHM/(sqrt(2)*exp(1));
%% 1. determine Ms, Me, and Mf
% Ms = min(max( abs( US_FindCOR(GG(:,:,floor(end*2/3):end)) ) ,0),1);
Ms = real(US_FindCOR(GG(:,:,floor(end*2/3):end)));
Me =1-abs(GG(:,:,1));
MfR = max(1-Ms-Me,0);
%% 2. determine Vz0 and fitting
Vz0=zeros(nz,nx); Vx=Vz0; Vy=Vz0; Vz=Vz0; R=Vz0;
GGf=zeros(nz,nx,nTau);
for iz=1:nz
    for ix=1:nx
        gg0=squeeze(GG(iz,ix,:));
        gg0_real=movmean(real(gg0),5);
        gg0_imag=movmean(imag(gg0),5);
        if abs(gg0(1))>0.2
             ms=Ms(iz,ix);mf0=MfR(iz,ix);me=Me(iz,ix);
             %% 2.1. determine Vz0 using the real part
             [r, lags]=xcorr(gg0_real);
             L_r=length(r);
             rDiffSign=[sign(diff(r(ceil(end/2):end)));-1*sign(diff(r(end-1:end)))];
             indSignChange=3+find(rDiffSign(4:end)>0);
             FirstVally_r_index=indSignChange(1);
             Tvz=lags(FirstVally_r_index-2+ceil(L_r/2))/rFrame*4/2;
             %                 Vz0(iz,ix)=2*pi/q/Tvz*sign(mean(diff(gg0_img(1:floor(lags(min(FirstVally_r_index+ceil(L_r/2),floor(L_r)))/2)))));
             Vz0(iz,ix)=2*pi/q/Tvz*sign(mean((gg0_imag(1:floor(lags(FirstVally_r_index-2+ceil(L_r/2))/2)))));
             if abs(Vz0(iz,ix))<0.0005
                 Vz0(iz,ix)=0.0005;  % m/s
                 Tvz=0;
             end
             vz0=Vz0(iz,ix);
             %% 2.2. Vt0 initial
             ggVt=gg0(1:10);tauVt=tau0(1:10);
             gvt=abs(ggVt-ms)/mf0.*exp((vz0*tauVt').^2/(Sigma2(3))^2);
             vt0=min(sqrt(abs(sum(Sigma2(1)^2*log(gvt))/sum(tauVt.^2))),10);
%              vt0=vz0;
             vx0=abs(vt0); vy0=abs(vt0);
             %% 2.3 constrain
             Fmin_cstrn(:,1)=[mf0-0.01 mf0+0.02];   % Mf constrain
             Fmin_cstrn(:,2)=(vt0+[-15 3]*1e-3)*tau0(end)/(Sigma2(1));  % Vx constrain
             Fmin_cstrn(:,3)=[0 0]*1e-3*tau0(end)/(Sigma2(2));  % Vy constrain
             Fmin_cstrn(:,4)=vz0+[-3 3]*1e-3;  % Vz constrain
             %% 2.4 fitting
             L_ggFit=floor(nTau);
             gg=gg0(1:L_ggFit);
             tau=tau0(1:L_ggFit);
             t = tau.';
             tn = t / tau(end);
             %% fit real part
             fitE = @(c) sum( abs(ms + c(1).*exp( -(c(2)*tn).^2-(c(3)*tn).^2-(c(4)*t).^2/(Sigma2(3))^2).*cos(2*k0*c(4)*t) - real(gg) ).^2 );
             fitC0 = [mf0, vx0*tau0(end)/(Sigma2(1)), vy0*tau0(end)/(Sigma2(2)), vz0];
             [fitC, fval] = fmincon(fitE, fitC0, [],[],[],[], ...
                 [Fmin_cstrn(1,1) Fmin_cstrn(1,2) Fmin_cstrn(1,3) Fmin_cstrn(1,4)], ...
                 [Fmin_cstrn(2,1) Fmin_cstrn(2,2) Fmin_cstrn(2,3) Fmin_cstrn(2,4)], ...
                 [], optimset('Display','off','TolFun',1e-6,'TolX',1e-6));%
             %% fit complex
             %         fitE = @(c) sum( abs(ms + c(1).*exp( -(c(2)*tn).^2-(c(3)*tn).^2-(c(4)*t).^2/(Sigma(3))^2).*exp(1i*2*k0*c(4)*t) - (gg) ).^2 );
             %         fitC0 = [mf0, vx0*tau0(end)/(Sigma(1)), vy0*tau0(end)/(Sigma(2)), vz0];
             %         [fitC, fval] = fmincon(fitE, fitC0, [],[],[],[], ...
             %             [Fmin_cstrn(1,1) Fmin_cstrn(1,2) Fmin_cstrn(1,3) Fmin_cstrn(1,4)], ...
             %             [Fmin_cstrn(2,1) Fmin_cstrn(2,2) Fmin_cstrn(2,3) Fmin_cstrn(2,4)], ...
             %             [], optimset('Display','off','TolFun',1e-6,'TolX',1e-6));%
             mf=fitC(1);
             vx = abs(fitC(2))*Sigma2(1)/tau0(end);
             vy = abs(fitC(3))*Sigma2(1)/tau0(end);
             vz = sign(vz0)*abs(fitC(4));
             %% fit imaginary part
             vz0=vz; vx0=vx; vy0=vy;
             ms=Ms(iz,ix);mf0=max(gg0_imag);me=Me(iz,ix);
             Fmin_cstrn(:,1)=[mf0-0.05 1];   % Mf constrain
             Fmin_cstrn(:,2)=(vx0+[-1 1]*1e-3)*tau0(end)/(Sigma2(1));  % Vx constrain
             Fmin_cstrn(:,3)=[0 0]*1e-3*tau0(end)/(Sigma2(2));  % Vy constrain
             Fmin_cstrn(:,4)=(vz0+[-0.5 0.5]*1e-3);  % Vz constrain
             fitE = @(c) sum( abs( c(1).*exp( -(c(2)*tn).^2-(c(3)*tn).^2-(c(4)*t).^2/(Sigma2(3))^2).*sin(2*k0*c(4)*t) - imag(gg) ).^2 );
             fitC0 = [mf0, vx0*tau0(end)/(Sigma2(1)), vy0*tau0(end)/(Sigma2(2)), vz0];
             [fitC, fval] = fmincon(fitE, fitC0, [],[],[],[], ...
                 [Fmin_cstrn(1,1) Fmin_cstrn(1,2) Fmin_cstrn(1,3) Fmin_cstrn(1,4)], ...
                 [Fmin_cstrn(2,1) Fmin_cstrn(2,2) Fmin_cstrn(2,3) Fmin_cstrn(2,4)], ...
                 [], optimset('Display','off','TolFun',1e-6,'TolX',1e-6));%
             mfI=fitC(1);
             vx = abs(fitC(2))*Sigma2(1)/tau0(end);
             vy = abs(fitC(3))*Sigma2(1)/tau0(end);
             vz = sign(vz0)*abs(fitC(4));
             %% %%%%%%%%%%%
             ggFit=ms + mf*exp(-vx^2*t.^2/(Sigma2(1))^2-vy^2*t.^2/(Sigma2(2))^2-vz^2*t.^2/(Sigma2(3))^2).*cos(2*k0*vz*t)+...
                 1i*mfI*exp(-vx^2*t.^2/(Sigma2(1))^2-vy^2*t.^2/(Sigma2(2))^2-vz^2*t.^2/(Sigma2(3))^2).*sin(2*k0*vz*t);
             Vx(iz,ix)=vx*1e3; Vy(iz,ix)=vy*1e3; Vz(iz,ix)=vz*1e3;
             Mf(:,iz,ix)=[mf mfI];
             R(iz,ix) = 1 - sum( abs(gg-ggFit).^2 ) / sum( abs(gg-mean(gg)).^2);
             GGf(iz,ix,1:L_ggFit)=ggFit;
        else
            Vx(iz,ix)=0; Vy(iz,ix)=0; Vz(iz,ix)=0;
            Mf(:,iz,ix)=[0 0];
            R(iz,ix) = 1;
            GGf(iz,ix,:)=gg0;
        end
        %% check fitted gg
        %         tCoor=tau*1e3;
        %         fig=figure;
        %         set(fig,'Position',[300 100 800 600])
        %         subplot(2,2,1), plot(tCoor,real(squeeze(gg)),'*'); hold on; plot(tCoor,real(squeeze(ggFit)));xlabel('time lag, [ms]');
        %         title(['Real of g1, Vx=',num2str(vx*1e3),' mm/s, Vy=', num2str(vy*1e3),' mm/s, Vz=', num2str(vz*1e3), ' mm/s'])
        %         subplot(2,2,2), plot(tCoor,imag(squeeze(gg)),'*'); hold on; plot(tCoor,imag(squeeze(ggFit))); xlabel('time lag, [ms]');
        %         title(['Imaginary of g1'])
        %         subplot(2,2,3), plot((squeeze(gg)),'*'); hold on; plot((squeeze(ggFit)));
        %         title(['Complex g1'])
        %         subplot(2,2,4), plot(tCoor,abs(squeeze(gg)),'*'); hold on; plot(tCoor,abs(squeeze(ggFit)));xlabel('time lag, [ms]');
        %         title(['Magnitude of g1'])
    end
end
%%%% Jonghwan's method
% GGa = unwrap(angle(GG-repmat(Ms,[1 1 nTau])),[],3);
% T = repmat(reshape(tau,[1 1 nTau]),[nz nx 1]);
% Vz0 = sum(GGa.*T,3) ./ sum(T.^2,3) /q;
%% 3. 

%% 2.1. determine Vz0 using the imaginary part
%             if max(abs(gg0_img))>0.1
%                 [r, lags]=xcorr(gg0_img(1:min(40*rFrame/1000,nTau)));
%                 L_r=length(r);
%                 rDiffSign=[sign(diff(r(ceil(end/2):end)));-1*sign(diff(r(end-1:end)))];
%                 indSignChange=3+find(rDiffSign(4:end)>0);
%                 FirstVally_r_index=indSignChange(1);
%                 Tvz=lags(FirstVally_r_index-2+ceil(L_r/2))/rFrame*4/2;
% %                 Vz0(iz,ix)=2*pi/q/Tvz*sign(mean(diff(gg0_img(1:floor(lags(min(FirstVally_r_index+ceil(L_r/2),floor(L_r)))/2)))));
%                 Vz0(iz,ix)=2*pi/q/Tvz*sign(mean((gg0_img(1:5))));
%                 if abs(Vz0(iz,ix))<0.0005
%                     Vz0(iz,ix)=0.0005;  % m/s
%                     Tvz=0;
%                 end
%                 vz0=Vz0(iz,ix); vx0=abs(vz0); vy0=abs(vz0);
%                 ms=Ms(iz,ix);mf0=Mf(iz,ix);me=Me(iz,ix);
%                 Fmin_cstrn(:,1)=[mf0-0.01 mf0+0.02];   % Mf constrain
%                 Fmin_cstrn(:,2)=[abs(vz0*0.5e3) 40]*1e-3*tau0(end)/(Sigma2(1));  % Vx constrain
%                 Fmin_cstrn(:,3)=[0 0]*1e-3*tau0(end)/(Sigma2(2));  % Vy constrain
%                 Fmin_cstrn(:,4)=[vz0-3e-3 vz0+3e-3];  % Vz constrain
%                 
%             else
%                 Vz0(iz,ix)=0.0005;  % m/s
%                 Tvz=0;
%                 vz0=Vz0(iz,ix); vx0=abs(vz0); vy0=abs(vz0);
%                 ms=Ms(iz,ix);mf0=Mf(iz,ix);me=Me(iz,ix);
%                 Fmin_cstrn(:,1)=[mf0-0.01 mf0+0.02];   % Mf constrain
%                 Fmin_cstrn(:,2)=[0 50]*1e-3*tau0(end)/(Sigma2(1));  % Vx constrain
%                 Fmin_cstrn(:,3)=[0 0]*1e-3*tau0(end)/(Sigma2(2));  % Vy constrain
%                 Fmin_cstrn(:,4)=[vz0-0.5e-3 vz0+0.5e-3];  % Vz constrain
%             end