%% fit g1, complex
function [Ms, Mf, Vx, Vy, Vz,R,GGf]=US_g1fit_INVIVO(GG, A_FWHM, rFrame)
% GG(nz,nx,nTau)
% Res: X, Y, Z spatial resolution
% rFrame: frame rate, Hz
%% 0. constant
[nz,nx,nTau]=size(GG);
f0=18e6; % ultrasound transducer center frequency, Hz
C=1540;  % ultrasound speed in biological tissue, m/s
k0 = 2*pi*f0/C; % wave number
q=2*k0;
dt = 1/rFrame;  % frame interval, s
tau0 = [1:nTau]*dt; % time lag, s
Sigma=A_FWHM*0.7/(2*sqrt(2*log(2)));
Sigma2=2*Sigma;
Nint_ggR=5;
% Sigma2(3)=Sigma2(3)*0.7;
% Sigma=R_FWHM/(sqrt(2)*exp(1));
%% 1. determine Ms, Me, and Mf
% Ms = min(max( abs( US_FindCOR(GG(:,:,floor(end*2/3):end)) ) ,0),1);
Ms = min(max(real(US_FindCOR(GG(:,:,floor(end*1/2):end))),0),max(mean(real(GG(:,:,floor(end*2/3):end)),3),0));
Me =1-abs(GG(:,:,1));
MfR = max(1-Ms-Me,0);
%% 2. determine Vz0 and fitting
Vz0=zeros(nz,nx); Vx=Vz0; Vy=Vz0; Vz=Vz0; 
GGf=zeros(nz,nx,nTau);
warning off;
for iz=1:nz
    for ix=1:nx
        gg0=squeeze(GG(iz,ix,:));
        gg0_real=movmean(real(gg0),5);
        gg0_imag=movmean(imag(gg0),5);
%         [peakValues0, valleyInd0] = findpeaks(gg0_real(1:end));
        if abs(gg0(1))>0.2
            ms0=Ms(iz,ix);
            mfR0=MfR(iz,ix);
            mfI0=max(abs(imag(gg0)));
            me=Me(iz,ix);
            %% nTau
            L_ggFit=floor(nTau);
            gg=gg0(1:L_ggFit);
            tau=tau0(1:L_ggFit);
            t = tau.';
            tn = t / tau(end);
            RI=0;
            %% 2.1. determine Vz0 using the real part
            gg0R=interp(gg0_real,Nint_ggR);
            [r, lags]=xcorr(gg0R);
            % figure,plot(lags, r)
            L_r=length(r);
            rDiffSign=[sign(diff(r(ceil(end/2):end)));-1*sign(diff(r(end-1:end)))];
            indSignChange=2+find(rDiffSign(3:end)>0)-1;
            FirstVally_r_index=indSignChange(1);
            gg0RIvs=-1*gg0R;
            [peakValues, valleyInd] = findpeaks(gg0RIvs(11:end));
            if ~isempty(valleyInd) && numel(peakValues)<5
                HalfCyc=valleyInd(1)+10;
            else
                HalfCyc=FirstVally_r_index;
            end
            Tvz=HalfCyc/rFrame*4/2/Nint_ggR;
            SignImag=sign(imag(gg));
            SignImagChange=find(SignImag~=SignImag(1));
            if ~isempty(SignImagChange) && SignImagChange(1)>5
                VzSign=SignImag(1);
            else
                VzSign=sign(mean((gg0_imag(1:floor(HalfCyc/Nint_ggR/2)))));
            end
            
            if abs(VzSign)~=1
                VzSign=1;
            end
            %                 Vz0(iz,ix)=2*pi/q/Tvz*sign(mean(diff(gg0_img(1:floor(lags(min(FirstVally_r_index+ceil(L_r/2),floor(L_r)))/2)))));
            Vz0(iz,ix)=2*pi/q/Tvz*VzSign;
            if max(abs(imag(gg(1:floor(10e-3*rFrame)))))<0.1 && abs(Vz0(iz,ix))>10e-3
                Vz0(iz,ix)=0.0002*sign(Vz0(iz,ix));  % m/s
                Tvz=0;
            end
            vz0=Vz0(iz,ix);
            vz0=sign(vz0)*min(abs(vz0),25e-3);
            %% 2.2. Vt0 initial
            vxmx0=[0:0.5:30]*1e-3;
            vy0=0;
            clear Rvt0;
            for ivx0=1:length(vxmx0)
%                 vy0=vxmx0(ivx0);
                ggvx0=ms0 + min(mfR0,1)*exp(-vxmx0(ivx0)^2*t.^2/(Sigma2(1))^2-vy0^2*t.^2/(Sigma2(2))^2-vz0^2*t.^2/(Sigma2(3))^2).*cos(2*k0*vz0*t);
                Rvt0(ivx0)=1 - sum( abs(real(gg)-ggvx0).^2 ) / sum( abs(real(gg)-mean(real(gg))).^2);
            end
            MInd=find(Rvt0==max(Rvt0));
            if isempty(MInd)
                MInd=1;
            end
            MInd=max(MInd,1);
            DR=max(Rvt0)-min(Rvt0(1:MInd));
            % test
%             saveFileBase='D:\OneDrive\Work\PROJ - FUS\PROJ - US velocimetry\Numerical Simulation\FlowAngleEffect\V5A0-10-90-1each\Rvt0-';
%             iFile=1;
%             saveFile=[saveFileBase,num2str(iFile),'.mat'];
%             while exist(saveFile,'file')==2
%                 iFile=iFile+1;
%                 saveFile=[saveFileBase,num2str(iFile),'.mat'];
%             end
%             save(saveFile,'Rvt0');
%             
%             figure,plot(vxmx0*1e3,Rvt0)
%             xlabel('Vx [mm/s]');
%             ylabel('R')
%             ylim([0.5 1])
%             grid on
            % determine vx0 -0
            if VzSign~=sign(mean((gg0_imag(1:floor(HalfCyc/Nint_ggR/2)))))
                vx0=min(min(vxmx0(MInd(1)),max((DR/0.02)^1.5*abs(vz0),3e-3)),5e-3*mfR0);
            else
                vx0=min(min(vxmx0(MInd(1)),max((DR/0.02)^1.5*abs(vz0),3e-3)),20e-3*mfR0);
            end
            ggvx0=ms0 + min(mfR0,1)*exp(-vx0^2*t.^2/(Sigma2(1))^2-vy0^2*t.^2/(Sigma2(2))^2-vz0^2*t.^2/(Sigma2(3))^2).*cos(2*k0*vz0*t);
            R0=1 - sum( abs(real(gg)-ggvx0).^2 ) / sum( abs(real(gg)-mean(real(gg))).^2);
            %% fit imaginary g1
            if  max(abs(imag(gg(1:floor(10e-3*rFrame)))))>0.08
%                 ggI=imag(gg(1:min(floor(HalfCyc/Nint_ggR*2.5),nTau)));
                ggI=gg0_imag(1:min(floor(HalfCyc/Nint_ggR*2.5),nTau));
                tauI=[1:length(ggI)]*dt; % time lag, s
                tI = tauI.';
                tnI = tI / tauI(end);
                %% 2.3 Imag-g1 constrain
                Fmin_cstrn(:,1)=[mfI0*1.2-0.05 min((mfI0+1),1)];   % Mf constrain
                Fmin_cstrn(:,2)=(vx0+[-1*min((2*(0.9/Rvt0(MInd))^4),vx0*1e3) min(2*(0.9/Rvt0(MInd))^4,3)]*1e-3)*tauI(end)/(Sigma2(1));  % Vx constrain
                Fmin_cstrn(:,3)=[0 0]*1e-3*tauI(end)/(Sigma2(2));  % Vy constrain
%                 Fmin_cstrn(:,3)=(vx0+[-1*min((2*(0.9/Rvt0(MInd))^4),vx0) min(2*(0.9/Rvt0(MInd))^4,5)]*1e-3)*tauI(end)/(Sigma2(2));  % Vy constrain
                Fmin_cstrn(:,4)=vz0+[-1 1]*1e-3;  % Vz constrain
                %% 2.4 Imag-g1 fitting
                fitE = @(c) sum( abs( c(1).*exp( -(c(2)*tnI).^2-(c(3)*tnI).^2-(c(4)*tI).^2/(Sigma2(3))^2).*sin(2*k0*c(4)*tI) - ggI ).^2 );
                fitC0 = [mfI0*1.2, vx0*tauI(end)/(Sigma2(1)), vy0*tauI(end)/(Sigma2(2)), vz0];
                [fitC, fval] = fmincon(fitE, fitC0, [],[],[],[], ...
                    [Fmin_cstrn(1,1) Fmin_cstrn(1,2) Fmin_cstrn(1,3) Fmin_cstrn(1,4)], ...
                    [Fmin_cstrn(2,1) Fmin_cstrn(2,2) Fmin_cstrn(2,3) Fmin_cstrn(2,4)], ...
                    [], optimset('Display','off','TolFun',1e-6,'TolX',1e-6));%
                mfI=fitC(1);
                vxI = abs(fitC(2))*Sigma2(1)/tauI(end);
                vyI = abs(fitC(3))*Sigma2(2)/tauI(end);
                vzI = sign(vz0)*abs(fitC(4));
                %% %%%%%%%%%%%
                ggFitR=ms0 + min(mfR0,1)*exp(-vxI^2*t.^2/(Sigma2(1))^2-vyI^2*t.^2/(Sigma2(2))^2-vzI^2*t.^2/(Sigma2(3))^2).*cos(2*k0*vzI*t);
                RI = 1 - sum( abs(real(gg)-ggFitR).^2 ) / sum( abs(real(gg)-mean(real(gg))).^2);
                % figure,plot(ggI,'.'); hold on, plot(ggFitI);
            end
            if RI>R0 && RI>0
                mfR0=MfR(iz,ix); vz0=vzI;
                %% 2.5. Vt0 initial
%                 vy0=vxI;
                vy0=0;
%                 vxmx1=vxI+[-5*min(max(1*(0.9/RI)^4,0.2),5):0.1:min(max(1*(0.9/RI)^4,0.2),2)]*1e-3;
                vxmx1=vxI+[-1*abs(vxI*1e3):0.5:5]*1e-3;
                clear Rvt1;
                for ivx0=1:length(vxmx1)
                    ggvx0=ms0 + min(mfR0,1)*exp(-vxmx1(ivx0)^2*t.^2/(Sigma2(1))^2-vy0^2*t.^2/(Sigma2(2))^2-vz0^2*t.^2/(Sigma2(3))^2).*cos(2*k0*vz0*t);
                    Rvt1(ivx0)=1 - sum( abs(real(gg)-ggvx0).^2 ) / sum( abs(real(gg)-mean(real(gg))).^2);
                end
            
%                 for ivx0=1:length(vxmx1)
%                     ggvx0= mfI*exp(-vxmx1(ivx0)^2*t.^2/(Sigma2(1))^2-vy0^2*t.^2/(Sigma2(2))^2-vz0^2*t.^2/(Sigma2(3))^2).*sin(2*k0*vz0*t);
%                     Rvt1(ivx0)=1 - sum( abs(imag(gg)-ggvx0).^2 ) / sum( abs(imag(gg)-mean(imag(gg))).^2);
%                 end
                [~, MInd1]=max(Rvt1);
                if isempty(MInd1)
                    MInd1=1;
                end
                if ~isempty(MInd1) &&  floor(Rvt1(MInd1)*1000)>floor(Rvt0(MInd)*1000)
                    %                 vx0=min(min(vxmx0(MInd(1)),100e-6/abs(vz0)),30e-3);
                    vx0=min(vxmx1(MInd1),15e-3*mfR0);
                    R1=max(Rvt1);
                else
                    vx0=max(min(vxI,15e-3*mfR0),3e-3);
                    R1=RI;
                end
                if R1<=0
                    R1=0.5;
                end
                mfI0=mfI;
                Fmin_cstrn(:,1)=[ms0-0.05 ms0+0.5];   % Ms constrain
                Fmin_cstrn(:,2)=[max(mfR0-0.0,0) min(mfR0+0.05,1)];   % MfR constrain
                Fmin_cstrn(:,3)=[max(mfI0-0.05,0) min(mfI0+0.05,1)];   % MfI constrain
                Fmin_cstrn(:,4)=(vx0+[-10*min(max(0.1*(0.9/R1)^4,0),0.2),2*min(max(0.1*(0.9/R1)^4,0),0.2)]*vx0)*tau0(end)/(Sigma2(1));  % Vx constrain
                Fmin_cstrn(:,5)=[0 0]*1e-3*tau0(end)/(Sigma2(2));  % Vy constrain
%                 Fmin_cstrn(:,5)=(vx0+[-15*min(max(0.1*(0.9/R1)^4,0),0.5),2*min(max(0.1*(0.9/R1)^4,0),0.5)]*vx0)*tau0(end)/(Sigma2(2));  % Vy constrain
                Fmin_cstrn(:,6)=vz0+[-1 1]*1e-3;  % Vz constrain
            else
                mfI0=max(abs(imag(gg0)));
                mfR0=MfR(iz,ix);
%                 vz0=sign(vz0)*min(abs(vz0),2e-3);
                vy0=0;
                R1=Rvt0(MInd);
                if R1<=0
                    R1=0.3;
                end
                Fmin_cstrn(:,1)=[ms0-0.1 ms0+0.05];   % Ms constrain
                Fmin_cstrn(:,2)=[max(mfR0-0.00, 0) min(mfR0+0.05,1)];   % MfR constrain
                Fmin_cstrn(:,3)=[mfI0-0.05 mfI0+0.5];   % MfI constrain
                Fmin_cstrn(:,4)=(vx0+[-10*min(max(0.5*(0.9/R1)^4,0),1),min(5*max(0.5*(0.9/R1)^4,0),0.5*abs(vx0*1e3))]*1e-3)*tau0(end)/(Sigma2(1));  % Vx constrain
                Fmin_cstrn(:,5)=[0 0]*1e-3*tau0(end)/(Sigma2(2));  % Vy constrain
%                 Fmin_cstrn(:,5)=(vx0+[-10*min(max(0.5*(0.9/R1)^4,0),1),2*min(max(0.5*(0.9/R1)^4,0),1)]*1e-3)*tau0(end)/(Sigma2(2));  % Vy constrain
                Fmin_cstrn(:,6)=sign(vz0)*(abs(vz0)+[-1*abs(vz0) 0.5*abs(vz0)]);  % Vz constrain
            end
            %% fit real(g1)
%             fitE = @(c) sum( abs(c(1) + c(2).*exp( -(c(3)*tn).^2-(c(4)*tn).^2-(c(5)*t).^2/(Sigma2(3))^2).*cos(2*k0*c(5)*t) - real(gg) ).^2 );
%             fitC0 = [ms0,mfR0, vx0*tau0(end)/(Sigma2(1)), vy0*tau0(end)/(Sigma2(2)), vz0];
%             [fitC, fval] = fmincon(fitE, fitC0, [],[],[],[], ...
%                 [Fmin_cstrn(1,1) Fmin_cstrn(1,2) Fmin_cstrn(1,4) Fmin_cstrn(1,5) Fmin_cstrn(1,6)], ...
%                 [Fmin_cstrn(2,1) Fmin_cstrn(2,2) Fmin_cstrn(2,4) Fmin_cstrn(2,5) Fmin_cstrn(2,6)], ...
%                 [], optimset('Display','off','TolFun',1e-6,'TolX',1e-6));%
%             %% fit result
%             ms=fitC(1);
%             mfR=fitC(2);
%             vx = abs(fitC(3))*Sigma2(1)/tau0(end);
%             vy = abs(fitC(4))*Sigma2(2)/tau0(end);
%             vz = sign(vz0)*abs(fitC(5));
%             ggFit=ms + mfR*exp(-vx^2*t.^2/(Sigma2(1))^2-vy^2*t.^2/(Sigma2(2))^2-vz^2*t.^2/(Sigma2(3))^2).*cos(2*k0*vz*t)...
%                 +1i*mfI*exp(-vx^2*t.^2/(Sigma2(1))^2-vy^2*t.^2/(Sigma2(2))^2-vz^2*t.^2/(Sigma2(3))^2).*sin(2*k0*vz*t);
%             R(1,iz,ix) = 1 - sum( abs((gg)-ggFit).^2 ) / sum( abs((gg)-mean((gg))).^2);
%             Mf(:,iz,ix)=[mfI mfR];
%             GGf(iz,ix,1:L_ggFit)=ggFit;
%             Vx(iz,ix)=vx*1e3; Vy(iz,ix)=vy*1e3; Vz(iz,ix)=vz*1e3;
            %% fit complex (g1)
            fitE = @(c) sum( abs(c(1) + c(2).*exp( -(c(4)*tn).^2-(c(5)*tn).^2-(c(6)*t).^2/(Sigma2(3))^2).*cos(2*k0*c(6)*t)+...
                1i*c(3).*exp( -(c(4)*tn).^2-(c(5)*tn).^2-(c(6)*t).^2/(Sigma2(3))^2).*sin(2*k0*c(6)*t)- gg ).^2 );
            fitC0 = [ms0,mfR0,mfI0, vx0*tau0(end)/(Sigma2(1)), vy0*tau0(end)/(Sigma2(2)), vz0];
            [fitC, fval] = fmincon(fitE, fitC0, [],[],[],[], ...
                [Fmin_cstrn(1,1) Fmin_cstrn(1,2) Fmin_cstrn(1,3) Fmin_cstrn(1,4) Fmin_cstrn(1,5) Fmin_cstrn(1,6)], ...
                [Fmin_cstrn(2,1) Fmin_cstrn(2,2) Fmin_cstrn(2,3) Fmin_cstrn(2,4) Fmin_cstrn(2,5) Fmin_cstrn(2,6)], ...
                [], optimset('Display','off','TolFun',1e-6,'TolX',1e-6));%
            %% fit result
            ms=fitC(1);
            mfR=fitC(2);
            mfI=fitC(3);
            vx = abs(fitC(4))*Sigma2(1)/tau0(end);
            vy = abs(fitC(5))*Sigma2(2)/tau0(end);
            vz = sign(vz0)*abs(fitC(6));
            ggFit=ms + mfR*exp(-vx^2*t.^2/(Sigma2(1))^2-vy^2*t.^2/(Sigma2(2))^2-vz^2*t.^2/(Sigma2(3))^2).*cos(2*k0*vz*t)...
                +1i*mfI*exp(-vx^2*t.^2/(Sigma2(1))^2-vy^2*t.^2/(Sigma2(2))^2-vz^2*t.^2/(Sigma2(3))^2).*sin(2*k0*vz*t);
            R(iz,ix) = 1 - sum( abs((gg)-ggFit).^2 ) / sum( abs((gg)-mean((gg))).^2);
            Mf(:,iz,ix)=[mfI mfR];
            GGf(iz,ix,1:L_ggFit)=ggFit;
            Vx(iz,ix)=vx*1e3; Vy(iz,ix)=vy*1e3; Vz(iz,ix)=vz*1e3;
        else
            Vx(iz,ix)=0; Vy(iz,ix)=0; Vz(iz,ix)=0;
            Mf(:,iz,ix)=[0 0];
            R(iz,ix) = 1;
            GGf(iz,ix,:)=gg0;
        end
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
%%%% Jonghwan's method
% GGa = unwrap(angle(GG-repmat(Ms,[1 1 nTau])),[],3);
% T = repmat(reshape(tau,[1 1 nTau]),[nz nx 1]);
% Vz0 = sum(GGa.*T,3) ./ sum(T.^2,3) /q;
