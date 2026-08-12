%% determine Vz0
% return vz0 in m/s
function [Vz0]=DtmVz0(gg, rFrame, lambda, nItp)
% nItp: resampling 
% rFrame: sampling frequency, Hz
% lambda: waelength, m
if nargin<4
    nItp=10;
end
%% I. resampling
nTau=length(gg);
sTau=linspace(1,nTau,nTau);
rsTau=linspace(1,nTau,nTau*nItp);
ggR=movmean(interp1(sTau,real(gg),rsTau,'linear'),50)';
ggI=movmean(interp1(sTau,imag(gg),rsTau,'linear'),50)';
rFrameItp=rFrame*nItp;
%% II. determine Vz0
if max(abs(imag(gg(1:floor(10e-3*rFrame)))))>=0.1
    %% 2. determine Vz0 value
    [r, lags]=xcorr(ggR);
    % figure,plot(lags, r)
    L_r=length(r);
    rDiffSign=[sign(diff(r(ceil(end/2):end)));-1*sign(diff(r(end-1:end)))];
    indSignChange=2+find(rDiffSign(3:end)>0)-1;
    FirstVally_r_index=indSignChange(1);
    gg0RIvs=-1*ggR;
    [peakValues, valleyInd] = findpeaks(gg0RIvs(11:end));
    if ~isempty(valleyInd) && numel(peakValues)<5
        HalfCyc=valleyInd(1)+10;
    else
        HalfCyc=FirstVally_r_index;
    end
    %% to correct Vz0 for some special cases
    if HalfCyc==floor(L_r/2) %|| max(abs(ggR(floor(HalfCyc*1.2):end)))>ggR(HalfCyc)*0.75
        vz0=0.002;
    else
        Tvz=HalfCyc/rFrameItp*2; % one cycle period
        vz0=lambda/(2*Tvz); % absolute value of velocity, m/s
    end
    %% 3. determine flow direction, negative: upwards flow; positive: downwards flow
    SignImag=sign(ggI);
    SignImagChange=find(SignImag~=SignImag(1));
    if ~isempty(SignImagChange) && SignImagChange(1)>5
        VzSign=SignImag(1);
    else
        VzSign=sign(mean((ggI(1:floor(HalfCyc/nItp/2)))));
    end
    if abs(VzSign)~=1
        VzSign=1;
    end
    %% 4. Vz0 coded with flow direction
    Vz0=vz0*VzSign;
else
    VzSign=sign(mean(ggI(1:5)));
    Vz0=0.0005*VzSign;  % m/s
    Tvz=0;
end



% %% determine Vz0
% % return vz0 in m/s
% function [Vz0]=DtmVz0(gg, rFrame, lambda, nItp)
% % nItp: resampling 
% % rFrame: sampling frequency, Hz
% % lambda: waelength, m
% if nargin<4
%     nItp=10;
% end
% %% I. resampling
% nTau=length(gg);
% sTau=linspace(1,nTau,nTau);
% rsTau=linspace(1,nTau,nTau*nItp);
% % ggR=movmean(interp1(sTau,real(gg),rsTau,'linear'),30)';
% % ggI=movmean(interp1(sTau,imag(gg),rsTau,'linear'),30)';
% ggR=interp1(sTau,movmean(real(gg),5),rsTau,'linear')';
% ggI=interp1(sTau,movmean(imag(gg),5),rsTau,'linear')';
% rFrameItp=rFrame*nItp;
% 
% %% II. determine Vz0
% if max(abs(imag(gg(1:floor(10e-3*rFrame)))))>=0.1
%     %% 2. determine Vz0 value
%     [r, lags]=xcorr(ggR);
%     % figure,plot(lags, r)
%     L_r=length(r);
%     rDiffSign=[sign(diff(r(ceil(end/2):end)));-1*sign(diff(r(end-1:end)))];
%     indSignChange=2+find(rDiffSign(3:end)>0)-1;
%     FirstVally_r_index=indSignChange(1);
%     gg0RIvs=-1*ggR;
%     [peakValues, valleyInd] = findpeaks(gg0RIvs(11:end));
%     if ~isempty(valleyInd) && numel(peakValues)<5
%         HalfCyc=valleyInd(1)+10;
%     else
%         HalfCyc=FirstVally_r_index;
%     end
%     %% to correct Vz0 for some special cases
%     if HalfCyc==floor(L_r/2) %|| max(abs(ggR(min(floor(HalfCyc*2),nTau*nItp):end)))>abs(ggR(HalfCyc))*0.75
%         vz0=0.002;
%     else
%         Tvz=HalfCyc/rFrameItp*2; % one cycle period
%         vz0=lambda/(2*Tvz); % absolute value of velocity, m/s
%     end
%     %% 3. determine flow direction, negative: upwards flow; positive: downwards flow
%     SignImag=sign(ggI);
%     SignImagChange=find(SignImag~=SignImag(1));
%     if ~isempty(SignImagChange) && SignImagChange(1)>5
%         VzSign=SignImag(1);
%     else
%         VzSign=sign(mean((ggI(1:floor(HalfCyc/nItp/2)))));
%     end
%     if abs(VzSign)~=1
%         VzSign=1;
%     end
%     %% 4. Vz0 coded with flow direction
%     Vz0=vz0*VzSign;
% else
%     VzSign=sign(mean(ggI(1:5)));
%     Vz0=0.001*VzSign;  % m/s
%     Tvz=0;
% end


