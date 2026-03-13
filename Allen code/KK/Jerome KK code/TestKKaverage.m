%%
dipstart

%% Read raw data 

FileNameBase='D:\OBUS_GE9L\';

%VMat='CIRS3';
VMat='abdomen4';

numTXangles=15;

[RawData,TXangle,TXdelay,TXapod]=RXunscrambleN(append(FileNameBase,VMat),numTXangles);  %raw data

Pitch = 0.00023;    % m
SamplingFrequency = 4*(250/48)*10^6;   % Hz
c = 1540;   % m/s
NA=0.2;     % numerical aperture
s=2*Pitch*SamplingFrequency/c;  %aspect ratio
xSize=192;

zBegin=100;  
zSize=1000;

RawDataComplex=DataHilbert(RawData);

%% Conventional complex DAS

[Recon,ReconSum]=BeamformAnglesComplex(RawData,zBegin,zSize,TXangle,NA,s);  %beamformed data
aReconSum=abs(mean(Recon,3)).^2;

dipshow3D(10,(aReconSum).^0.125)

%% KK-DAS with coherent or incoherent averaging
iNum=8;  %number of instances

rangeRXangles =12;    %half total range in degrees  (TX range = +/-12 degrees)
%numRXangles = 32;      %not needed here

dtheta = (2*rangeRXangles*pi/180)/(xSize-1);
startAngle = rangeRXangles*pi/180;

angles=-startAngle:dtheta:startAngle;  %calculate 192 distinct angles

randindex=randperm(192);  
randangles=angles(randindex);  %randomly permute angles

AnglesPerInstance=xSize/iNum;
RXA=zeros(AnglesPerInstance,iNum);

% Two methods to bin angles

% Method 1:
% for n=1:iN
%     RXA(:,n)=randangles((n-1)*AnglesPerInstance+1:n*AnglesPerInstance);
% end

% Method 2
for n=1:AnglesPerInstance
    randindex((n-1)*iNum+1:n*iNum)=(n-1)*iNum+randperm(iNum);
end
for n=1:iNum
    RXA(:,n)=angles(randindex(n:iNum:xSize));
end

% Cycle through instances
iKK=zeros(zSize,xSize,iNum); %intensity
eKK=zeros(zSize,xSize,iNum)+i*zeros(zSize,xSize,iNum);  % field
for ii=1:iNum   %cycle through instances
ii
RXangle=squeeze(RXA(:,ii));

RawDataKKComplex=DataCompressKK(RawDataComplex,RXangle,s);
[eConstructComplex,eConstructSumComplex,eDelay]=...
    BeamformKK(RawDataKKComplex,zBegin,zSize,xSize,TXangle,RXangle,s);
KKreconComplex=sum(eConstructComplex,[3,4]);

eKK(:,:,ii)=KKreconComplex;
iKK(:,:,ii)=abs(KKreconComplex).^2;
end

iKKcoh=abs(mean(eKK,3)).^2;  %coherent summation over instances
iKKinc=mean(abs(eKK).^2,3);  %incoherent summation over instances
vKK=var(eKK,[],3);  % field variance
KKcoh=(iKKcoh-vKK);  % this doesn't work well
KKinc=(iKKinc-vKK);  % this seems to correct haze in iKKinc and do better than iKKcoh

dipshow3D(10,iKKcoh.^0.125)  
dipshow3D(20,iKKinc.^0.125)
dipshow3D(30,vKK.^0.125)
dipshow3D(40,KKcoh.^0.125)
dipshow3D(50,KKinc.^0.125)

