%% US g1 fit for in vivo data, fit negative and postive frequency signal separately
% input: 
    % sIQ: bulk motion removed data, [nz,nx,nt]
    % PRSSinfo: data processing parameters, including 
        % PRSSinfo.FWHM: (X, Y, Z) spatial resolution, Full Width at Half Maximum of point spread function, m
        % PRSSinfo.rFrame: sIQ frame rate, Hz
        % PRSSinfo.f0: Transducer center frequency, Hz
        % PRSSinfo.C: Sound speed in the sample, m/s
        % PRSSinfo.g1nT: g1 calculation sample number
        % PRSSinfo.g1nTau: maximum number of time lag
        % PRSSinfo.SVDrank: SVD rank [low high]
        % PRSSinfo.HPfC:  High pass filtering cutoff frequency, Hz
        % PRSSinfo.NEQ: do noise equalization? 0: no noise equalization; 1: apply noise equalization
        % PRSSinfo.rfnScale: spatial refind scale
        % PRSSinfo.MpVz: maximu pVz
        % PRSSinfo.useMsk: 1: use ULM data as spatial mask; 0: no spatial mask
        % PRSSinfo.ulmMsk: ULM-based spatial constrain mask
            % [nz,nx,3], 1: up flow (positive frequency); 2 down flow (negative
            % frequency); 3: all flow 
            % ulmMsk=1 otherwise
 % output:
    % Ms: static component fraction, [nz,nx]
    % Mf: dynamic component fraction, [nz,nx,2], 2: [real,imag]
    % Vx: x-direction velocity component, [nz,nx], mm/s
    % Vz: axial-direction velocity component, [nz,nx], mm/s
    % V=sqrt(Vx.^2+Vz.^2), [nz,nx], mm/s
    % pVz: Vz distribution (sigma-Vz), [nz,nx]
    % R: fitting accuracy, [nz,nx]
    % CR: freqCR.*pnCR.*MskCR.*ggCR
    % GGf: gg fitting results, [nz,nx, nTau]
 % Jianbo Tang, 20190820
function Vcz=sIQ2Vcz_NP(sIQ, PRSSinfo)
%% O. constant
lambda0=PRSSinfo.C/PRSSinfo.f0;        % wavlength
k0 = 2*pi/lambda0;   % wave number
PRSSinfo.FWHM=[PRSSinfo.FWHM(1) 1000e-6 PRSSinfo.FWHM(2)]; % just to put the FWHM_y, any number other than 0
Sigma=PRSSinfo.FWHM*0.7/(2*sqrt(2*log(2))); % intensity-based sigma
Sigma2=2*Sigma;
nItpVz0=10;          % for Vz0 determination
dt = 1/PRSSinfo.rFrame;      % frame interval, s
tau = [1:PRSSinfo.g1nTau]*dt; % time lag, s
tn = tau / tau(end);
%% I. determine spectrum power ratio and signal (|f|<1000Hz) to noise ratio
sIQ=gpuArray(sIQ);
[nz0,nx0,nt]=size(sIQ);
% I.0 all frequency signal and SNR
PRSSinfo.g1StartT=1;
GG0 = sIQ2GG(sIQ, PRSSinfo); % g1 of whole frequency signal
fCoor=linspace(-PRSSinfo.rFrame/2,PRSSinfo.rFrame/2,nt)';
fCoorSig=zeros(size(fCoor));
fCoorSig(abs(fCoor)<1100)=1; % signal frequency range
fCoorSig=circshift(fCoorSig,nt/2);
fIQ=(fft(sIQ,nt,3)); % no fft shift
% I.1 frequency-based SNR
zPix=max(floor(nz0*0.1),1):1:floor(nz0*1);
fSNR0=squeeze(sum(abs(fIQ.*repmat(permute(fCoorSig,[3 2 1]),[nz0 nx0 1])),3))./squeeze(sum(abs(fIQ),3)); % SNR of oringla data
zfSNR=PRSSinfo.useMsk*(mean(fSNR0,2)-1.3*std(fSNR0,[],2))+(1-PRSSinfo.useMsk)*(mean(fSNR0,2)-0.9*(1+([1:nz0]./(5*nz0)).^2)'.*std(fSNR0,[],2));
fC=polyfit(zPix,zfSNR(max(floor(nz0*0.1),1):floor(nz0*1))',1);
fSNRthd0=repmat(polyval(fC,[1:nz0])',[1, nx0])*0.95;     
% I.2 all frequency spatial mask
if PRSSinfo.useMsk==1
    MskCR=(abs(PRSSinfo.ulmMsk(:,:,3))>0.5);
else
    MskCR=1;
end
% I.3 binary mask - all frequency
aCR=(((fSNR0>fSNRthd0)+abs(GG0(:,:,1))>0.4)>0).*(abs(GG0(:,:,1))>(PRSSinfo.useMsk*0.1+(1-PRSSinfo.useMsk)*0.15)); % all frequency signal-based thresholding
clear GG0 sIQ
%% II. positive&negative frequncy signal vUS processing
% PRSSinfo.MpVz=0.8; % maximu pVz
DispPrss={'Positive frequency signal', 'Negative frequency signal'};
for iNP=1:2
    iFIQ=zeros(size(fIQ),'gpuArray');
    iFIQ(:,:,(floor(nt/2))*(iNP-1)+1:floor(nt/2)*iNP)=fIQ(:,:,(floor(nt/2))*(iNP-1)+1:floor(nt/2)*iNP);
    %% III. acceptable signal criteria, CR
    % III. 1 freqCR
    iFIQ_S=iFIQ.*repmat(permute(fCoorSig,[3 2 1]),[nz0 nx0 1]);
    fSNR0=squeeze(sum(abs(iFIQ_S),3))./squeeze(sum(abs(iFIQ),3)); % SNR of frequency signal
    zfSNR=PRSSinfo.useMsk*(mean(fSNR0,2)-1.1*std(fSNR0,[],2))+(1-PRSSinfo.useMsk)*(mean(fSNR0,2)-0.9*(1+([1:nz0]./(5*nz0)).^2)'.*std(fSNR0,[],2));
    fC=polyfit(zPix,zfSNR(max(floor(nz0*0.1),1):floor(nz0*1))',1);
    fSNRthd0=repmat(polyval(fC,[1:nz0])',[1, nx0])*(0.9+(iNP-1)*0.1);
    fCR=(fSNR0>fSNRthd0);
    fCR0=(fSNR0>0.55);
    % III. 2 pnCR
    ipnRatio=sum(abs(iFIQ(:,:,(floor(nt/2))*(iNP-1)+1:floor(nt/2)*iNP))-repmat(median(abs(fIQ(:,:,:)),3),[1 1 floor(nt/2)]),3); % spectrum power of positive frequency
    ipnRatio=ipnRatio./(sum(abs(fIQ(:,:,:))-repmat(median(abs(fIQ(:,:,:)),3),[1 1 nt]),3));
    pnCR=(ipnRatio>((0.3-PRSSinfo.useMsk*0.1)+(-0.1+0.1*(iNP-1))));
    % III. 3 iGG and ggCR
    iIQ=(ifft(iFIQ,nt,3));
    iGG=sIQ2GG(iIQ, PRSSinfo); % g1 of p or n frequency signal
%     disp(['Color Doppler Processing - ', datestr(datetime('now'))]);
    iVcz=gather(ColorDoppler(iIQ,PRSSinfo)); % color Doppler
    ggCR=(abs(iGG(:,:,1))>(PRSSinfo.useMsk*0.1+(1-PRSSinfo.useMsk)*0.2));
    % III. 4 CR for positive or negative frequency signal
    iCR=gather(((pnCR+fCR0)>0).*fCR.*ggCR.*aCR);  % acceptable signal criteria
    Vcz0(:,:,iNP)=iVcz.*iCR;
end
% III. 3 MaskCR
if PRSSinfo.rfnScale>1
    for iNP=1:2
        if PRSSinfo.useMsk==1
            mskCR=PRSSinfo.ulmMsk(:,:,iNP);
        else
            mskCR=1;
        end

        Vcz(:,:,iNP)=imresize(Vcz0(:,:,iNP),[nz0,nx0]*PRSSinfo.rfnScale,'nearest').*mskCR; % spatial interpolation
        Vcz(isnan(abs(Vcz)))=0;
    end

else
    Vcz=Vcz0;
end
Vcz=-1*Vcz;

            
