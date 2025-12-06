%% determine Vx0 and PVz0 if forced to be 0
% Input:
    % GG: [nVox,nTau], nVox=nz*nx*ny
    % Vz0: 1D array, [nVox,1], initial value, m/s
    % Ms0: 1D array, [nVox,1], initial value for Ms
    % Mf0: 1D array, [nVox,1], initial value for Mf
    % PRSinfo: processing information
    % PRSinfo.FWHM: (T, Z), m
    % PRSinfo.fAline: DAQ Aline rate, Hz
    % PRSinfo.Lam: [light source center, wavelength bandwidth], m
    % PRSinfo.Dim: GG dimension: [nz,nx,ny,nTau]
% Output: 
    % Vt: 1D array, [nVox,1], initial guess for Vt, in m/s
    % Vz: 1D array, [nVox,1], initial guess for Vz, in m/s
    % D: 1D array, [nVox,1], initial guess for D, m^2/s
    % R: 1D array, [nVox,1], fitting accuracy with initial guesses
% Jianbo Tang, 20190731
function [Vt,Vz,D,R]=iniDLSOCT_GPU(GG, Vz0, Ms0, Mf0, PRSinfo)
%% I. DAQ parameter
Sigma=PRSinfo.FWHM*0.7/(2*sqrt(2*log(2))); % intensity-based sigma, full width at the 1/e maximum value
Sigma2=2*Sigma;
k0 = 2*pi/PRSinfo.Lam(1);   % wave number, /m
n=1.35; % refractive index
q=2*n*k0;
[nVox,nTau]=size(GG); % nVox=nz*nx*ny
nz=PRSinfo.Dim(1); nx=PRSinfo.Dim(2); ny=PRSinfo.Dim(3); 
if nz*nx*ny>100*400*2
    if rem(ny,2)==0
        nyPchk=2;
    else
        nyPchk=1;
    end
    nyChk=ny/nyPchk;
else
    nyPchk=ny;
    nyChk=ny/nyPchk;
end
nVoxPchk=nz*nx*nyPchk;
%% GPU array
GG=gpuArray(GG);
CR=(Mf0>0.1); %(abs(GG(:,1))>abs(GG(:,2))).*(abs(GG(:,1))>0.5); % threshold criteria
Vz0=gpuArray(Vz0);  % m/s
Ms0=gpuArray(Ms0);
Mf0=gpuArray(Mf0);
%% II. fit for a meshgrid of Vt and D
mVt=gpuArray(zeros(nVoxPchk,11,nTau)); mVz=gpuArray(zeros(nVoxPchk,5,nTau));
tau=[1:nTau]/PRSinfo.fAline; % time lag, s
for iChk=1:nyChk
    iVoxStart=(iChk-1)*nVoxPchk+1;
    iVoxEnd=iChk*nVoxPchk;
    StepD=5;  % um^2/s
    NmVt=11;
    mVt0=[0:NmVt-1]/(NmVt-1);
    [mVt,mVt00,~]=ndgrid(min(60./abs(Vz0(iVoxStart:iVoxEnd)*1e3),15)*1e-3.*Mf0(iVoxStart:iVoxEnd).^1.5,mVt0,tau);
    mVt=mVt.*mVt00;
    StepVt=mVt(:,2,1)-mVt(:,1,1);
    mD=[0:StepD:90]*1e-12; % m^2/s, diffusion coefficient
    NmD=length(mD); 
    [mVz, ~]=ndgrid(Vz0(iVoxStart:iVoxEnd),tau);
    [mMs, ~]=ndgrid(Ms0(iVoxStart:iVoxEnd),tau);
    [mMf, Tau]=ndgrid(Mf0(iVoxStart:iVoxEnd),tau);
    GGi=gpuArray(single(zeros(nVoxPchk,nTau,NmVt,NmD)));
%     %% complex g1
%     for iVt=1:NmVt
%         for iD=1:NmD
%             GGi(:,:,iVt,iD)=mMs + mMf.*exp(-permute(mVt(:,iVt,:),[1,3,2]).^2.*Tau.^2/(Sigma2(1))^2-mVz.^2.*Tau.^2/(Sigma2(2))^2)...
%                 .*exp(-q^2*Tau.*mD(iD)).*exp(1i*q*mVz.*Tau);
%         end
%     end
%     RR=1 - (sum( abs(repmat((GG(iVoxStart:iVoxEnd,:)),[1,1,NmVt,NmD])-GGi).^2,2))...
%         ./ repmat(sum( abs((GG(iVoxStart:iVoxEnd,:))-mean((GG(iVoxStart:iVoxEnd,:)),2)).^2,2),[1,1,NmVt,NmD]);
    %% real part of g1 
    for iVt=1:NmVt
        for iD=1:NmD
            GGi(:,:,iVt,iD)=mMs + mMf.*exp(-permute(mVt(:,iVt,:),[1,3,2]).^2.*Tau.^2/(Sigma2(1))^2-mVz.^2.*Tau.^2/(Sigma2(2))^2)...
                .*exp(-q^2*Tau.*mD(iD)).*cos(q*mVz.*Tau);
        end
    end
    RR=1 - (sum( abs(repmat(real(GG(iVoxStart:iVoxEnd,:)),[1,1,NmVt,NmD])-GGi).^2,2))...
        ./ repmat(sum( abs(real(GG(iVoxStart:iVoxEnd,:))-mean(real(GG(iVoxStart:iVoxEnd,:)),2)).^2,2),[1,1,NmVt,NmD]); 
    
    RR=permute(RR,[1,3,4,2]);
    [mR0,RI]=max(RR(:,:),[],2);
    [MIvt, MId]=ind2sub([NmVt,NmD],RI);
    Vt(iVoxStart:iVoxEnd,1)=gather(mVt(:,1,1)+(MIvt-1).*StepVt); % m/s
    D(iVoxStart:iVoxEnd,1)=gather((MId-1)*StepD.*CR(iVoxStart:iVoxEnd))*1e-12; % m^2/s
    R(iVoxStart:iVoxEnd,1)=gather(mR0);
end
Vz=gather(Vz0);

%%%%%%%%% II. fit for a meshgrid of Vt and D
% StepVt=1; % mm/s
% StepD=5;  % um^2/s
% mVt=[0:StepVt:16]*1e-3; % m/s
% mD=[0:StepD:80]*1e-12; % m^2/s, diffusion coefficient
% NmVt=length(mVt); NmD=length(mD);
% tau=[1:nTau]/PRSinfo.fAline; % time lag, s
% for iChk=1:nyChk
%     [mVz, Tau]=ndgrid(Vz0((iChk-1)*nVoxPchk+1:iChk*nVoxPchk),tau);
%     [mMs, Tau]=ndgrid(Ms0((iChk-1)*nVoxPchk+1:iChk*nVoxPchk),tau);
%     [mMf, Tau]=ndgrid(Mf0((iChk-1)*nVoxPchk+1:iChk*nVoxPchk),tau);
%     GGi=gpuArray(single(zeros(nVoxPchk,nTau,length(mVt),length(mD))));
% %     %% complex g1
% %     for iVt=1:length(mVt)
% %         for iD=1:length(mD)
% %             GGi(:,:,iVt,iD)=mMs + mMf.*exp(-mVt(iVt).^2.*Tau.^2/(Sigma2(1))^2-mVz.^2.*Tau.^2/(Sigma2(2))^2).*exp(-q^2*Tau.*mD(iD)).*exp(1i*q*mVz.*Tau);
% %         end
% %     end
% %     RR=1 - (sum( abs(repmat((GG((iChk-1)*nVoxPchk+1:iChk*nVoxPchk,:)),[1,1,length(mVt),length(mD)])-GGi).^2,2))...
% %         ./ repmat(sum( abs((GG((iChk-1)*nVoxPchk+1:iChk*nVoxPchk,:))-mean((GG((iChk-1)*nVoxPchk+1:iChk*nVoxPchk,:)),2)).^2,2),[1,1,length(mVt),length(mD)]);
%     %% real part of g1 
%     for iVt=1:length(mVt)
%         for iD=1:length(mD)
%             GGi(:,:,iVt,iD)=mMs + mMf.*exp(-mVt(iVt).^2.*Tau.^2/(Sigma2(1))^2-mVz.^2.*Tau.^2/(Sigma2(2))^2).*exp(-q^2*Tau.*mD(iD)).*cos(q*mVz.*Tau);
%         end
%     end
%     RR=1 - (sum( abs(repmat(real(GG((iChk-1)*nVoxPchk+1:iChk*nVoxPchk,:)),[1,1,length(mVt),length(mD)])-GGi).^2,2))...
%         ./ repmat(sum( abs(real(GG((iChk-1)*nVoxPchk+1:iChk*nVoxPchk,:))-mean(real(GG((iChk-1)*nVoxPchk+1:iChk*nVoxPchk,:)),2)).^2,2),[1,1,length(mVt),length(mD)]);
%     
%     RR=permute(RR,[1,3,4,2]);
%     [mR0,RI]=max(RR(:,:),[],2);
%     [MIvt, MId]=ind2sub([NmVt,NmD],RI);
%     Vt((iChk-1)*nVoxPchk+1:iChk*nVoxPchk,1)=gather((min((MIvt-1)*StepVt*1e-3,...
%         min(60./abs(Vz0((iChk-1)*nVoxPchk+1:iChk*nVoxPchk)),15).*Mf0((iChk-1)*nVoxPchk+1:iChk*nVoxPchk).^1.5*1e-3))...
%         .*CR((iChk-1)*nVoxPchk+1:iChk*nVoxPchk)); % m/s
%     D((iChk-1)*nVoxPchk+1:iChk*nVoxPchk,1)=gather((MId-1)*StepD.*CR((iChk-1)*nVoxPchk+1:iChk*nVoxPchk))*1e-12; % m^2/s
%     R((iChk-1)*nVoxPchk+1:iChk*nVoxPchk,1)=gather(mR0);
% end
% Vz0=gather(Vz0);

%%%%%%%%%%% II. fit for a meshgrid of Vz, Vt and D
% mVt=gpuArray(zeros(nVoxPchk,11,nTau)); mVz=gpuArray(zeros(nVoxPchk,5,nTau));
% tau=[1:nTau]/PRSinfo.fAline; % time lag, s
% for iChk=1:nyChk
%     iVoxStart=(iChk-1)*nVoxPchk+1;
%     iVoxEnd=iChk*nVoxPchk;
%     %% Vt initial grid
%     NmVt=11;
%     mVt0=[0:NmVt-1]/(NmVt-1);
%     [mVt,mVt00,~]=ndgrid(min(60./abs(Vz0(iVoxStart:iVoxEnd)*1e3),15)*1e-3.*Mf0(iVoxStart:iVoxEnd).^1.5,mVt0,tau);
%     mVt=mVt.*mVt00;
%     StepVt=mVt(:,2,1)-mVt(:,1,1);
%     %% Vz initial grid
%     NmVz=3;
%     mVz0=linspace(0.8,1.2, NmVz);
%     [mVz,mVz00,~]=ndgrid(Vz0(iVoxStart:iVoxEnd),mVz0,tau);
%     mVz=mVz.*mVz00;
%     StepVz=mVz(:,2,1)-mVz(:,1,1);
%     %% D initial array
%     StepD=5;  % um^2/s
%     mD=[0:StepD:80]*1e-12; % m^2/s, diffusion coefficient
%     NmD=length(mD);
%     
%     GGi=gpuArray(single(zeros(nVoxPchk,nTau,NmVz,NmVt,NmD)));
%     [mMs, ~]=ndgrid(Ms0(iVoxStart:iVoxEnd),tau);
%     [mMf, Tau]=ndgrid(Mf0(iVoxStart:iVoxEnd),tau);
%     %% complex g1
%     for iVz=1:NmVz
%         for iVt=1:NmVt
%             for iD=1:NmD
%                 GGi(:,:,iVz,iVt,iD)=mMs + mMf.*exp(-permute(mVt(:,iVt,:),[1,3,2]).^2.*Tau.^2/(Sigma2(1))^2-permute(mVz(:,iVz,:),[1,3,2]).^2.*Tau.^2/(Sigma2(2))^2)...
%                     .*exp(-q^2*Tau.*mD(iD)).*exp(1i*q*permute(mVz(:,iVz,:),[1,3,2]).*Tau);
%             end
%         end
%     end
%     RR=1 - (sum( abs(repmat((GG(iVoxStart:iVoxEnd,:)),[1,1,NmVz, NmVt,NmD])-GGi).^2,2))...
%         ./ repmat(sum( abs((GG(iVoxStart:iVoxEnd,:))-mean((GG(iVoxStart:iVoxEnd,:)),2)).^2,2),[1,1,NmVz, NmVt,NmD]);
% %     %% real part of g1 
% %     for iVz=1:NmVz
% %         for iVt=1:NmVt
% %             for iD=1:NmD
% %                 GGi(:,:,iVz,iVt,iD)=mMs + mMf.*exp(-permute(mVt(:,iVt,:),[1,3,2]).^2.*Tau.^2/(Sigma2(1))^2-permute(mVz(:,iVz,:),[1,3,2]).^2.*Tau.^2/(Sigma2(2))^2)...
% %                     .*exp(-q^2*Tau.*mD(iD)).*cos(q*permute(mVz(:,iVz,:),[1,3,2]).*Tau);
% %             end
% %         end
% %     end
% %     RR=1 - (sum( abs(repmat(real(GG(iVoxStart:iVoxEnd,:)),[1,1,NmVz, NmVt,NmD])-GGi).^2,2))...
% %         ./ repmat(sum( abs(real(GG(iVoxStart:iVoxEnd,:))-mean(real(GG(iVoxStart:iVoxEnd,:)),2)).^2,2),[1,1,NmVz, NmVt,NmD]);
% 
%     RR=permute(RR,[1,3,4,5,2]);
%     [mR0,RI]=max(RR(:,:),[],2);
%     [MIvz, MIvt, MId]=ind2sub([NmVz,NmVt,NmD],RI);
%     Vz(iVoxStart:iVoxEnd,1)=gather(mVz(:,1,1)+(MIvz-1).*StepVz); % m/s
%     Vt(iVoxStart:iVoxEnd,1)=gather(mVt(:,1,1)+(MIvt-1).*StepVt); % m/s
%     D(iVoxStart:iVoxEnd,1)=gather((MId-1)*StepD.*CR(iVoxStart:iVoxEnd))*1e-12; % m^2/s
%     R(iVoxStart:iVoxEnd,1)=gather(mR0);
% end