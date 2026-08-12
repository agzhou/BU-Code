%% vUS fitting
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
% output:
    % Ms: static component fraction, [nz,nx]
    % Mf: dynamic component fraction, [nz,nx,2], 2: [real,imag]
    % Vx: x-direction velocity component, [nz,nx], mm/s
    % Vz: axial-direction velocity component, [nz,nx], mm/s
    % V=sqrt(Vx.^2+Vz.^2), [nz,nx], mm/s
    % R: fitting accuracy, [nz,nx]
    % GGf: GG fitting results, [nz,nx, nTau]
% Jianbo Tang, 20190812
function [Vz,Vx,pVz,Ms,Mf,R,GGf]=GG2vUS_GPU(GG, Vz0, Ms0, MfR0, PRSSinfo)
%% I. DAQ parameter
[nVox,nTau]=size(GG);
tau=[1:nTau]/PRSSinfo.rFrame; % time lag, s
Sigma=PRSSinfo.FWHM*0.7/(2*sqrt(2*log(2))); % intensity-based sigma
Sigma2=2*Sigma;
lambda0=PRSSinfo.C/PRSSinfo.f0;        % wavlength
k0 = 2*pi/lambda0;   % wave number
%% approx decay time
% GGtd=sign(max(min(abs(GG),[],2),0.2)-abs(GG));
% [~,Ind]=max(GGtd,[],2);
% Tdc=Ind/PRSSinfo.rFrame;
%% II. Vx, Vz and PVz mesh
NmVx=10;
mVx0=linspace(0,1,NmVx); 
if PRSSinfo.MpVz==0
%     [mVx,mVx00,~]=ndgrid((30*3e-3./Tdc*1e-3),mVx0,tau);
    [mVx,mVx00,~]=ndgrid((40e-3),mVx0,tau);
    NmPVz=1;
    mPVz=0; StepPvz=0;
else
    sVmsk=(Vz0>-5e-3).*(Vz0<5e-3);
    VxMax1=abs(Vz0)*3;                      % when -5<Vz0<3 mm/s
    VxMax2=min(180./abs(Vz0*1e3),30)*1e-3;  % when Vz0>3 or <-5 mm/s
    VxMax=((VxMax1.*sVmsk)+VxMax2.*(1-sVmsk)).*MfR0;
    VxMax=movmean(VxMax,5,1);
    [mVx,mVx00,~]=ndgrid(VxMax,mVx0,tau);
    
    NmPVz=5;
    mPVz=linspace(0.3,PRSSinfo.MpVz,NmPVz);
% %     mPVzMax1=max(min(abs(Vz0)/20e-3,6e-3./abs(Vz0)),0.2);
% %     mPVzMax2=0.8;
% %     mPVzMax=((mPVzMax2.*sVmsk)+mPVzMax1.*(1-sVmsk));
%     mPVzMax=ones(nVox,1)*0.8;
%     [mPVz,mPVz00,~]=ndgrid(mPVzMax,mPVz0,tau);
%     mPVz=mPVz.*mPVz00;
%     StepPvz=mPVz(:,2,1)-mPVz(:,1,1);
end
mVx=mVx.*mVx00;
StepVx=mVx(:,2,1)-mVx(:,1,1);
clear mVx00 mPVz00 mPVzMax

[mMs, ~]=ndgrid(Ms0,tau);
[mMf, Tau]=ndgrid(MfR0,tau);
%% real part of g1 
NmVz=9;
mVz0=linspace(0.5, 1.2, NmVz);                      
[mVz,mVz00]=ndgrid(Vz0,mVz0);
mVz1=mVz.*mVz00;
clear mVz00;
StepVz=mVz1(:,2)-mVz1(:,1);
RR=zeros(nVox,NmVz,NmVx,NmPVz,'gpuArray');
for iVz=1:NmVz
    iVz0=mVz1(:,iVz);
    [mVz, ~]=ndgrid(iVz0,tau);
    for iVx=1:NmVx
        for iPvz=1:NmPVz
            GGr(:,:,iPvz)=mMs + mMf.*exp(-permute(mVx(:,iVx,:),[1,3,2]).^2.*Tau.^2/(Sigma2(1))^2-mVz.^2.*Tau.^2/(Sigma2(3))^2)...
                .*exp(-(k0*Tau.*mPVz(iPvz).*mVz).^2).*cos(2*k0*mVz.*Tau);
        end
        RR(:,iVz,iVx,:)=1 - permute((sum( abs(repmat(real(GG),[1,1,NmPVz])-GGr).^2,2)...
            ./ repmat(sum( abs(real(GG)-mean(real(GG),2)).^2,2),[1,1,NmPVz])),[1,2,4,3]);
    end
end
clear GGr;
[R0,InD]=max(RR(:,:),[],2);
[MIvz, MIvx, MIPvz]=ind2sub([NmVz, NmVx,NmPVz],InD);
Vz0=mVz1(:,1)+(MIvz-1).*StepVz; % m/s
Vx0=(mVx(:,1,1)+(MIvx-1).*StepVx).*(R0>0)+((80e-3./abs(Vz0))*1e-3.*MfR0).*(R0<=0); % m/s
% Vx0(abs(Vx0)>35e-3)=10e-3;
% Vx0=min((mVx(:,1,1)+(MIvx-1).*StepVx),(200e-3./abs(Vz0))*1e-3.*MfR0);
% pVz0(:,1)=(mPVz(:,1,1)+(MIPvz-1).*StepPvz);
pVz0(:,1)=mPVz(MIPvz);
%% III. determine MfI0 from the imag GG
NmMfI=10;
mMfI0=linspace(0.6, 1.4, NmMfI);
[mMfI,mMfI00]=ndgrid(max(abs(imag(GG)),[],2),mMfI0);
mMfI1=mMfI.*mMfI00;
StepMfI=mMfI1(:,2)-mMfI1(:,1);
for imMfI=1:NmMfI
    GGi(:,:,imMfI)=mMfI1(:,imMfI).*exp(-Vx0.^2.*tau.^2/(Sigma2(1))^2-Vz0.^2.*tau.^2/(Sigma2(3))^2).*exp(-(k0*tau.*pVz0.*Vz0).^2).*sin(2*k0*Vz0.*tau);
end
RRi=1 - (sum( abs(repmat(imag(GG),[1,1,NmMfI])-GGi).^2,2))...
    ./ repmat(sum( abs(imag(GG)-mean(imag(GG),2)).^2,2),[1,1,NmMfI]);
[RI0,MindMfI]=max(RRi,[],3);
MindMfI0=MindMfI(:,1);
MfI0=min(mMfI1(:,1)+(MindMfI0-1).*StepMfI,1); % m/s
clear mVx mVx00 mPvz mMs mMf mMfI mMfI00 GGi RR
%% IV. g1 fit, finer grid
NmVx=5;
mVx0=linspace(0.8, 1.5, NmVx);
[mVx,mVx00,~]=ndgrid(Vx0,mVx0,tau);
mVx=mVx.*mVx00;
StepVx=mVx(:,2,1)-mVx(:,1,1);

NmPVz=5;
mpVz0=linspace(0.8,1.2,NmPVz);
[mpVz,mpVz00,~]=ndgrid(pVz0,mpVz0,tau);
mpVz=mpVz.*mpVz00;
StepPvz=mpVz(:,2,1)-mpVz(:,1,1);

[mMs, ~]=ndgrid(Ms0,tau);
[mMfR, ~]=ndgrid(MfR0,tau);
[mMfI, ~]=ndgrid(MfI0,tau);
[mVz, Tau]=ndgrid(Vz0,tau);
%% complex part of g1 
RR=zeros(nVox,NmVx,NmPVz,'gpuArray');
for iVx=1:NmVx
    for iPvz=1:NmPVz
        GGc(:,:,iPvz)=mMs + exp(-permute(mVx(:,iVx,:),[1,3,2]).^2.*Tau.^2/(Sigma2(1))^2-mVz.^2.*Tau.^2/(Sigma2(3))^2)...
            .*exp(-(k0*Tau.*permute(mpVz(:,iPvz,:),[1,3,2]).*mVz).^2).*(mMfR.*cos(2*k0*mVz.*Tau)+mMfI.*1i.*sin(2*k0*mVz.*Tau));
    end
    RR(:,iVx,:)=1 - permute((sum( abs(repmat((GG),[1,1,NmPVz])-GGc).^2,2)...
    ./ repmat(sum( abs((GG)-mean((GG),2)).^2,2),[1,1,NmPVz])),[1,2,3]);
end
clear GG GGc
[R,indR]=max(RR(:,:),[],2);
[MIvx, MIPvz]=ind2sub([NmVx,NmPVz],indR);
Vx=reshape((mVx(:,1,1)+(MIvx-1).*StepVx)*1e3, [PRSSinfo.Dim(1), PRSSinfo.Dim(2)]); % mm/s
pVz=reshape((mpVz(:,1,1)+(MIPvz-1).*StepPvz), [PRSSinfo.Dim(1), PRSSinfo.Dim(2)]); %
Vz=reshape(Vz0*1e3, [PRSSinfo.Dim(1), PRSSinfo.Dim(2)]); % mm/s
Ms=reshape(Ms0, [PRSSinfo.Dim(1), PRSSinfo.Dim(2)]);
R=reshape(R, [PRSSinfo.Dim(1), PRSSinfo.Dim(2)]);
Mf(:,:,1)=reshape(MfR0, [PRSSinfo.Dim(1), PRSSinfo.Dim(2)]);
Mf(:,:,2)=reshape(MfI0, [PRSSinfo.Dim(1), PRSSinfo.Dim(2)]);

[mVx1, ~]=ndgrid((mVx(:,1,1)+(MIvx-1).*StepVx),tau);
[mpVz1, ~]=ndgrid((mpVz(:,1,1)+(MIPvz-1).*StepPvz),tau);
GGf=mMs + exp(-mVx1.^2.*Tau.^2/(Sigma2(1))^2-mVz.^2.*Tau.^2/(Sigma2(3))^2)...
            .*exp(-(k0*Tau.*mpVz1.*mVz).^2).*(mMfR.*cos(2*k0*mVz.*Tau)+mMfI.*1i.*sin(2*k0*mVz.*Tau));
GGf=reshape(GGf,[PRSSinfo.Dim(1), PRSSinfo.Dim(2),PRSSinfo.Dim(3)]);


