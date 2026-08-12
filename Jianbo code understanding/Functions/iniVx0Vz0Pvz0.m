%% determine the initial guess of Vx0,Vz0,Pvz0
% Input:
    % GG, 2D array, (nVox,nTau), nVox=nz*nx
    % Vz0: 1D array, [nVox,1], initial value, m/s
    % Ms0: 1D array, [nVox,1], initial value for Ms
    % MfR0: 1D array, [nVox,1], initial value for Mf
    % PRSSinfo.C: sound speed, m/s
    % PRSSinfo.FWHM: (X, Y, Z) spatial resolution, Full Width at Half Maximum of point spread function, m
    % PRSSinfo.rFrame: sIQ frame rate, Hz
    % PRSSinfo.f0: Transducer center frequency, Hz
    % PRSSinfo.MpVz: maximu pVz
% Output: 
    % Vx: initial guess for Vx, in m/s
    % Vz: initial guess for Vz, in m/s
    % PVz: initial guess for PVz, in m/s
    % MfI: initial guess for Mf_imag
    % R: fitting accuracy with initial guesses
% Jianbo Tang, 20190812
function [Vx,Vz,PVz,MfI,R]=iniVx0Vz0Pvz0(GG, Vz0, Ms0, MfR0, PRSSinfo)
%% I. DAQ parameter
[nVox,nTau]=size(GG);
tau=[1:nTau]/PRSSinfo.rFrame; % time lag, s
Sigma=PRSSinfo.FWHM*0.7/(2*sqrt(2*log(2))); % intensity-based sigma
Sigma2=2*Sigma;
lambda0=PRSSinfo.C/PRSSinfo.f0;        % wavlength
k0 = 2*pi/lambda0;   % wave number
%% approx decay time
GGtd=sign(max(min(abs(GG),[],2),0.2)-abs(GG));
[~,Ind]=max(GGtd,[],2);
Tdc=Ind/PRSSinfo.rFrame;
%% II. Vx, Vz and PVz mesh
NmVx=10;
mVx0=[0:NmVx-1]/(NmVx-1);
if PRSSinfo.MpVz==0
%     [mVx,mVx00,~]=ndgrid((30*3e-3./Tdc*1e-3),mVx0,tau);
    [mVx,mVx00,~]=ndgrid((40e-3),mVx0,tau);
else
    sVmsk=(Vz0>-5e-3).*(Vz0<3e-3);
    VxMax1=abs(Vz0)*2;                      % when -5<Vz0<3 mm/s
    VxMax2=min(180./abs(Vz0*1e3),35)*1e-3;  % when Vz0>3 or <-5 mm/s
    VxMax=(VxMax1.*sVmsk)+VxMax2.*(1-sVmsk);
    [mVx,mVx00,~]=ndgrid(VxMax,mVx0,tau);
end
mVx=mVx.*mVx00;
StepVx=mVx(:,2,1)-mVx(:,1,1);

NmPvz=floor(PRSSinfo.MpVz/0.1+1);
mPvz=linspace(0,PRSSinfo.MpVz,NmPvz); %
% StepPvz=mPvz(:,2)-mPvz(:,1);

[mMs, ~]=ndgrid(Ms0,tau);
[mMf, Tau]=ndgrid(MfR0,tau);
%% real part of g1 
NmVz=5;
mVz0=[NmVz-3:0.8:NmVz+1]/(NmVz);
[mVz,mVz00]=ndgrid(Vz0,mVz0);
mVz1=mVz.*mVz00;
StepVz=mVz1(:,2)-mVz1(:,1);
for iVz=1:NmVz
    iVz0=mVz1(:,iVz);
    [mVz, ~]=ndgrid(iVz0,tau);
    for iVx=1:NmVx
        for iPvz=1:NmPvz
            GGr(:,:,iVx,iPvz)=mMs + mMf.*exp(-permute(mVx(:,iVx,:),[1,3,2]).^2.*Tau.^2/(Sigma2(1))^2-mVz.^2.*Tau.^2/(Sigma2(3))^2)...
                .*exp(-(k0*Tau.*mPvz(iPvz).*mVz).^2).*cos(2*k0*mVz.*Tau);
        end
    end
    RR=1 - (sum( abs(repmat(real(GG),[1,1,NmVx,NmPvz])-GGr).^2,2))...
        ./ repmat(sum( abs(real(GG)-mean(real(GG),2)).^2,2),[1,1,NmVx,NmPvz]);
    RR=permute(RR,[1,3,4,2]);
    [mR0(:,iVz),indR(:,iVz)]=max(RR(:,:),[],2);
end
[R,IndMRVz]=max(mR0,[],2);
IR=mR0-R;
IR(IR==0)=1; IR(IR~=1)=0;
InD=max(IR.*indR,[],2);
[MIvx, MIPvz]=ind2sub([NmVx,NmPvz],InD);
Vx=(mVx(:,1,1)+(MIvx-1).*StepVx); % m/s
PVz(:,1)=mPvz(MIPvz);
Vz=mVz1(:,1)+(IndMRVz-1).*StepVz; % m/s
%% III. determine MfI0 from the imag GG
NmMfI=5;
mMfI0=[NmMfI-0.5:0.75:NmMfI+3]/(NmMfI);
[mMfI,mMfI00]=ndgrid(max(abs(imag(GG)),[],2),mMfI0);
mMfI1=mMfI.*mMfI00;
StepMfI=mMfI1(:,2)-mMfI1(:,1);
for imMfI=1:NmMfI
    GGi(:,:,imMfI)=mMfI1(:,imMfI).*exp(-Vx.^2.*tau.^2/(Sigma2(1))^2-Vz.^2.*tau.^2/(Sigma2(3))^2).*exp(-(k0*tau.*PVz.*Vz).^2).*sin(2*k0*Vz.*tau);
end
RRi=1 - (sum( abs(repmat(imag(GG),[1,1,NmMfI])-GGi).^2,2))...
    ./ repmat(sum( abs(imag(GG)-mean(imag(GG),2)).^2,2),[1,1,NmMfI]);
[RI0,MindMfI]=max(RRi,[],3);
MindMfI0=MindMfI(:,1);
MfI=min(mMfI1(:,1)+(MindMfI0-1).*StepMfI,1); % m/s

