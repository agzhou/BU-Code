function [Ms, Mf, Vx, Vy, Vz, A, Rfit, GGf] = US_FitAcr(tau, GG, Res)

    % constant
    f0=18e6; % ultrasound transducer center frequency, Hz
    C=1540; 
    k0 = 2*pi*f0/C;
    q=2*k0;
    % variables
    [nz nx ntau] = size(GG);
    na = 2;
    Ms = zeros(nz,nx,na);  Mf = Ms;  D = Ms;  V = Ms;  A = Ms;  Rfit = zeros(nz,nx,na);  GGf = ones(nz,nx,ntau,na)*(1+1i);
    t = reshape(tau(1:end),[ntau 1]);		
    tn = t / tau(end);
		
    % 1st initial guess
%     Ms0 = min(max( abs( US_FindCOR(GG(:,:,floor(end/2):end)) ) ,0),1);
    Ms0=0.05;
    GGa = unwrap(angle(GG-repmat(Ms0,[1 1 ntau])),[],3);			
    T = repmat(reshape(tau,[1 1 ntau]),[nz nx 1]);
    Vz0 = sum(GGa.*T,3) ./ sum(T.^2,3) /q;
    clear GGa T;

    % 2nd initial guess
    X = [tau(1)^2 tau(1) 1; tau(2)^2 tau(2) 1; tau(3)^2 tau(3) 1];
    C = US_Multiply(inv(X),shiftdim(abs(GG(:,:,2:4)),2));
%     Me0 = 1 - min(max( shiftdim(C(3,:,:),1) ,0),1);
    Me0 =1-abs(GG(:,:,1));
    clear X C;
    Mf0 = max(1-Ms0-Me0,0);

%     nN = 45;
%     Ms1 = zeros(1,nN);
%     Mf1 = zeros(1,nN);
%     iN = 0;
%     for ms=1:9
%         for mf=1:9
%             if ms + mf <= 10
%                 iN = iN+1;
%                 Ms1(iN) = ms/10;
%                 Mf1(iN) = mf/10;
%             end
%         end
%     end
%     [DD VT VZ RR] = US_FitAcr_dv(tau, GG, Ms1, Mf1, q,h,hx);
%     [r IR] = max(RR,[],5);

    % fit at each voxel
    for iz=1:nz
        for ix=1:nx
                
                % variables
                g = squeeze(GG(iz,ix,1:end));
                msi = zeros(3,1);  mfi = msi;  di = msi;  vxi = msi; vyi = msi;  vzi = msi;
                msf0 = 1-Me0(iz,ix);
                vi0 = .5*Res(3)/tau(end);
                % three initial guess
                msi(1) = Ms0(iz,ix);  mfi(1) = Mf0(iz,ix);   vxi(1) = Vz0(iz,ix); vyi(1) = Vz0(iz,ix); vzi(1) = Vz0(iz,ix);
%                 msi(2) = Ms1(IR(iz,ix));  mfi(2) = Mf1(IR(iz,ix));   vxi(2) = VT(iz,ix,1,IR(iz,ix)); vyi(2) = VT(iz,ix,1,IR(iz,ix)); vzi(2) = VZ(iz,ix,IR(iz,ix));
                msi(2) = msf0/2;  mfi(2) = msf0/2;    vxi(2) = vi0/sqrt(2); vyi(2) = vi0/sqrt(3); vzi(2) = vi0/sqrt(3);

                % fminsearch
                for ja=1:1
%                     fitE2 = @(c) sum( abs( min(abs(c(1)),1) + min(abs(c(2)),1-min(abs(c(1)),1)).*exp( -(c(3)*tn).^2-(c(4)*tn).^2-(sum(unwrap(angle(g-min(abs(c(1)),1))).*t)/sum(t.^2)/(2*k0))^2*t.^2/(exp(-1)*Res(3))^2).*exp(1i*2*k0* sum(unwrap(angle(g-min(abs(c(1)),1))).*t)/sum(t.^2)/(2*k0) *t) - g ).^2 );
%                     fitC0 = [msi(ja) mfi(ja) vxi(ja)*tau(end)/(exp(-1)*Res(1)) vyi(ja)*tau(end)/(exp(-1)*Res(2))];
%                     fitC = fminsearch(fitE2, fitC0, optimset('Display','off'));
%                     ms = min(abs(fitC(1)),1);  mf = min(abs(fitC(2)),1-ms);  vx = abs(fitC(3))*exp(-1)*Res(1)/tau(end);vy = abs(fitC(4))*exp(-1)*Res(2)/tau(end);
%                     vz = sum(unwrap(angle(g-ms)).*t)/sum(t.^2)/(2*k0);
%                     v = sqrt(vx^2+vy^2+vz^2);
                    %%
                    fitE2 = @(c) sum( abs(Ms0 + Mf0.*exp( -(c(1)*tn).^2-(sum(unwrap(angle(g-Ms0)).*t)/sum(t.^2)/(2*k0))^2*t.^2/(exp(-1)*Res(3))^2).*exp(1i*2*k0* sum(unwrap(angle(g-Ms0)).*t)/sum(t.^2)/(2*k0) *t) - g ).^2 );
                    fitC0 = [vxi(ja)*tau(end)/(exp(-1)*Res(1))];
                    fitC = fminsearch(fitE2, fitC0, optimset('Display','off'));
                    ms = Ms0;  mf = Mf0;  vx = abs(fitC(1))*exp(-1)*Res(1)/tau(end);vy = 0; %abs(fitC(2))*exp(-1)*Res(2)/tau(end);
                    vz = sum(unwrap(angle(g-ms)).*t)/sum(t.^2)/(2*k0);
                    v = sqrt(vx^2+vy^2+vz^2);
                    %%

                    Ms(iz,ix,ja) = ms;  Mf(iz,ix,ja) = mf;  Vx(iz,ix,ja) = vx; Vy(iz,ix,ja) = vy; Vz(iz,ix,ja) = vz;
                    if (v > 0)  A(iz,ix,ja) = vz/v;  end;
                    
                    GGf(iz,ix,1:ntau,ja) = ms + mf*exp(-vx^2*t.^2/(exp(-1)*Res(1))^2-vy^2*t.^2/(exp(-1)*Res(2))^2-vz^2*t.^2/(exp(-1)*Res(3))^2).*exp(1i*2*k0*vz*t);  
                    gf(1)=1; gf(2:ntau+1)=squeeze(GGf(iz,ix,:,ja));
                    g0(1)=1; g0(2:ntau+1)= g;
                    Rfit(iz,ix,ja) = 1 - sum( abs(g0-gf).^2 ) / sum( abs(g0-mean(g0)).^2);	
%                     R(iz,ix,ja) = 1 - sum( abs(g0(2:ntau+1)-gf(2:ntau+1)).^2 ) / sum( abs(g0(2:ntau+1)-mean(g0(2:ntau+1))).^2);	
%                     R(iz,ix,ja) = 1 - sum( abs(GG(iz,ix,:)-GGf(iz,ix,:,ja)).^2 ,4) / sum( abs(GG(iz,ix,:)-mean(GG(iz,ix,:),4)).^2 ,4);									

                end

        end

    end

