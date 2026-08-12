%% determine Vx0 and PVz0
% return vz0 in m/s
function [Vx0,Vz1,PVz0,MfI0,R0,GGini]=DtmVx0nPVz0(gg, Vz0, Ms0, MfR0, k0,A_FWHM,rFrame)
%% I. DAQ parameter
nTau=length(gg);
tau=[1:nTau]/rFrame; % time lag, s
Sigma=A_FWHM*0.7/(2*sqrt(2*log(2))); % intensity-based sigma
Sigma2=2*Sigma;
%% II. Vx and PVz mesh
mVz=Vz0+[-1:0.5:3]*1e-3;
if Vz0*1e3>3 || -1*Vz0*1e3>5 
%     mVx=[0:0.5:min(150/(abs(Vz0)*1e3),25)]*1e-3;
    mVx=[0:0.5:min(200/(abs(Vz0)*1e3),25)]*1e-3;
    mPVz=max(min(abs(Vz0)/20e-3,6e-3/abs(Vz0)),0.2)+[-0.1:0.025:0.4];
%     mPVz=min(max(min(abs(Vz0)/20e-3,6e-3/abs(Vz0)),0.3),0.6)+[-0.1:0.05:0.2];
else
    mVx=[0:0.5:abs(Vz0)*1e3*2]*1e-3;
%     mVx=[0:0.5:10]*1e-3;
    mPVz=[0.2:0.1:0.8];
end

% mPVz=0.5; % test, constrain PVz=0.5 which makse Vx=0 at the axial flow vessels
Vy0=0;
clear Rvt0;
%% II.1 calculate the fitting accuracy for each pair
for imVz=1:length(mVz)
    for iPvz=1:length(mPVz)
        for iVx=1:length(mVx)
            mRGG=Ms0 + min(MfR0,1)*exp(-mVx(iVx)^2*tau.^2/(Sigma2(1))^2-Vy0^2*tau.^2/(Sigma2(2))^2-mVz(imVz)^2*tau.^2/(Sigma2(3))^2).*exp(-(k0*tau*mPVz(iPvz)*mVz(imVz)).^2).*cos(2*k0*mVz(imVz)*tau);
            imR(iPvz,iVx)=1 - sum( abs(real(gg)-mRGG.').^2 ) / sum( abs(real(gg)-mean(real(gg))).^2);
        end
    end
    [iMIndPvz, iMIndVx]=find(imR==max(imR(:)));
    if isempty(iMIndPvz)
        iMIndPvz=1;
        iMIndVx=1;
    end
    imaxR(imVz)=mean(max(imR(:)));
    imaxIndPvz(imVz)=iMIndPvz(1);
    imaxIndVx(imVz)=iMIndVx(1);
end
%% II.2 maximum mR indices
[MR0,MindVz0]=max(imaxR);
MindVz=MindVz0(1);
RR0=mean(MR0); % The best approximated fitting accuracy of the real GG
Vz1=mVz(MindVz);
Vx0=min(mVx(imaxIndVx(MindVz)),(150e-3/abs(Vz1))*1e-3*MfR0);
PVz0=mPVz(imaxIndPvz(MindVz));
%% III. determine MfI0 from the imag GG
mMfI=0:0.05:1;
for iMfI=1:length(mMfI)
    mIGG=mMfI(iMfI)*exp(-Vx0^2*tau.^2/(Sigma2(1))^2-Vy0^2*tau.^2/(Sigma2(2))^2-Vz0^2*tau.^2/(Sigma2(3))^2).*exp(-(k0*tau*PVz0*Vz1).^2).*sin(2*k0*Vz1*tau);
    mRI(iMfI)=1 - sum( abs(imag(gg)-mIGG.').^2 ) / sum( abs(imag(gg)-mean(imag(gg))).^2);
end
[mRI0,MindMfI]=max(mRI);
if isempty(MindMfI)
    MindMfI=5;
end
MindMfI0=MindMfI(1);
RI0=mean(mRI0); % The best approximated fitting accuracy of the real GG
MfI0=max(mMfI(MindMfI0),max(abs(imag(gg))));
%% IV. inital R0
GGini=Ms0+MfR0*exp(-Vx0^2*tau.^2/(Sigma2(1))^2-Vy0^2*tau.^2/(Sigma2(2))^2-Vz0^2*tau.^2/(Sigma2(3))^2).*exp(-(k0*tau*PVz0*Vz1).^2).*cos(2*k0*Vz1*tau)+...
    1i*MfI0*exp(-Vx0^2*tau.^2/(Sigma2(1))^2-Vy0^2*tau.^2/(Sigma2(2))^2-Vz0^2*tau.^2/(Sigma2(3))^2).*exp(-(k0*tau*PVz0*Vz1).^2).*sin(2*k0*Vz1*tau);
R0=1 - sum( abs((gg)-GGini.').^2 ) / sum( abs((gg)-mean((gg))).^2);

% %% determine Vx0 and PVz0
% % return vz0 in m/s
% function [Vx0,Vz1,PVz0,MfI0,R0]=DtmVx0nPVz0(gg, Vz0, Ms0, MfR0, k0,A_FWHM,rFrame)
% %% I. DAQ parameter
% nTau=length(gg);
% tau=[1:nTau]/rFrame; % time lag, s
% Sigma=A_FWHM*0.7/(2*sqrt(2*log(2))); % intensity-based sigma
% Sigma2=2*Sigma;
% %% II. Vx and PVz mesh
% mVz=Vz0+[-2:0.5:2]*1e-3;
% if Vz0*1e3>3 || -1*Vz0*1e3>5 
%     mVx=[0:0.5:min(150/(abs(Vz0)*1e3),25)]*1e-3;
% %     mVx=[0:0.5:min(200/(abs(Vz0)*1e3),25)]*1e-3;
% %     mPVz=min(max(min(abs(Vz0)/20e-3,6e-3/abs(Vz0)),0.4),0.6)+[-0.1:0.05:0.2];
%     mPVz=max(min(abs(Vz0)/20e-3,6e-3/abs(Vz0)),0.2)+[-0.1:0.025:0.4];
% else
%     mVx=[0:0.5:abs(Vz0)*1e3*2]*1e-3;
% %     mPVz=[0.5:0.1:1];
%     mPVz=[0.2:0.1:0.8];
% end
% % mPVz=0.5; % test, constrain PVz=0.5 which makse Vx=0 at the axial flow vessels
% Vy0=0;
% clear Rvt0;
% %% II.1 calculate the fitting accuracy for the mesh
% for imVz=1:length(mVz)
%     for iPvz=1:length(mPVz)
%         for iVx=1:length(mVx)
%             mRGG=Ms0 + min(MfR0,1)*exp(-mVx(iVx)^2*tau.^2/(Sigma2(1))^2-Vy0^2*tau.^2/(Sigma2(2))^2-mVz(imVz)^2*tau.^2/(Sigma2(3))^2).*exp(-(k0*tau*mPVz(iPvz)*mVz(imVz)).^2).*cos(2*k0*mVz(imVz)*tau);
%             imR(iPvz,iVx)=1 - sum( abs(real(gg)-mRGG.').^2 ) / sum( abs(real(gg)-mean(real(gg))).^2);
%         end
%     end
%     [iMIndPvz, iMIndVx]=find(imR==max(imR(:)));
%     if isempty(iMIndPvz)
%         iMIndPvz=1;
%         iMIndVx=1;
%     end
%     imaxR(imVz)=mean(max(imR(:)));
%     imaxIndPvz(imVz)=iMIndPvz(1);
%     imaxIndVx(imVz)=iMIndVx(1);
% end
% %% II.2 maximum mR indices
% [MR0,MindVz0]=max(imaxR);
% MindVz=MindVz0(1);
% RR0=mean(MR0); % The best approximated fitting accuracy of the real GG
% Vz1=mVz(MindVz);
% Vx0=min(mVx(imaxIndVx(MindVz)),(150e-3/abs(Vz1))*1e-3*MfR0);
% PVz0=mPVz(imaxIndPvz(MindVz));
% %% III. determine MfI0 from the imag GG
% mMfI=max(abs(imag(gg))):0.05:1;
% for iMfI=1:length(mMfI)
%     mIGG=mMfI(iMfI)*exp(-Vx0^2*tau.^2/(Sigma2(1))^2-Vy0^2*tau.^2/(Sigma2(2))^2-Vz0^2*tau.^2/(Sigma2(3))^2).*exp(-(k0*tau*PVz0*Vz1).^2).*sin(2*k0*Vz1*tau);
%     mRI(iMfI)=1 - sum( abs(imag(gg)-mIGG.').^2 ) / sum( abs(imag(gg)-mean(imag(gg))).^2);
% end
% [mRI0,MindMfI]=max(mRI);
% if isempty(MindMfI)
%     MindMfI=5;
% end
% MindMfI0=MindMfI(1);
% RI0=mean(mRI0); % The best approximated fitting accuracy of the real GG
% MfI0=max(mMfI(MindMfI0),max(abs(imag(gg))));
% %% IV. inital R0
% GGini=Ms0+MfR0*exp(-Vx0^2*tau.^2/(Sigma2(1))^2-Vy0^2*tau.^2/(Sigma2(2))^2-Vz0^2*tau.^2/(Sigma2(3))^2).*exp(-(k0*tau*PVz0*Vz1).^2).*cos(2*k0*Vz1*tau)+...
%     1i*MfI0*exp(-Vx0^2*tau.^2/(Sigma2(1))^2-Vy0^2*tau.^2/(Sigma2(2))^2-Vz0^2*tau.^2/(Sigma2(3))^2).*exp(-(k0*tau*PVz0*Vz1).^2).*sin(2*k0*Vz1*tau);
% R0=1 - sum( abs((gg)-GGini.').^2 ) / sum( abs((gg)-mean((gg))).^2);
