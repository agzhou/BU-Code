%%
dipstart

%% Read raw data 

FileNameBase='D:\OBUS_GE9L\';

%VMat='Sim15Cos2_Point';
%VMat='Sim15_Point';
%VMat='CIRS1';
VMat='abdomen4';

numTXangles=15;

[RawData,TXangle,TXdelay,TXapod]=RXunscrambleN(append(FileNameBase,VMat),numTXangles);  %raw data

Pitch = 0.00023;    % m
SamplingFrequency = 4*(250/48)*10^6;   % Hz
c = 1540;   % m/s
NA=0.2;     % numerical aperture
s=2*Pitch*SamplingFrequency/c;  %aspect ratio
xSize=192;

zBegin=1;  
zSize=1000;

%zBegin=600;  %good for resolution target
%zSize=250;   %good for resolution target

%zBegin=20;  %good for abdomen
%zSize=450;   %good for abdomen

RawDataComplex=DataHilbert(RawData);


%% Conventional complex DAS
tic
[Recon,ReconSum]=BeamformAnglesComplex(RawData,zBegin,zSize,TXangle,NA,s);  %beamformed data

aReconSum=abs(ReconSum);
toc

dipshow3D(10,(aReconSum).^0.25)

%% KK-DAS Complex
iN=16;  %number of images for min processing
kkmin=10^9;   %init kmin processing
% Step 1: Define RXangle for KK-DAS data compression
rangeRXangles =12;    %half total range in degrees  (TX range = +/-12 degrees)
numRXangles = 17;      %number of angle steps  (TX num = 15)


vKKraw=zeros(zSize,xSize,iN);
vKKmin=zeros(zSize,xSize,iN);
% aa=1.4*log(360/pi/rangeRXangles/xSize)/log(2/numRXangles);
% aaa=aa:-(aa-1/aa)/iN:1/aa;

for ii=1:iN
ii
%equidistant RX 
RXangle = zeros(numRXangles,1);
if (numRXangles > 1)
    dtheta = (2*rangeRXangles*pi/180)/(numRXangles-1);
    startAngle = -rangeRXangles*pi/180;
else
    dtheta = 0;
    startAngle=0;
end 

for n=1:numRXangles
    RXangle(n)=startAngle+(n-1)*dtheta;
end

%add linear randomness to RX angles 
if ii>1
    rr=4*dtheta*(rand(numRXangles,1)-0.5);
    RXangle=RXangle+rr;
end

% Step 2: Compress data
RawDataKKComplex=DataCompressKK(RawDataComplex,RXangle,s);
% Step 3: Beamform
[eConstructComplex,eConstructSumComplex,eDelay]=BeamformKK(RawDataKKComplex,zBegin,zSize,xSize,TXangle,RXangle,s);
KKreconComplex=sum(eConstructComplex,[3,4]);  %sum over all k's

dipshow3D(25,abs(KKreconComplex).^0.25);

kkmin=min(abs(KKreconComplex),kkmin);
dipshow3D(23,kkmin.^0.25);

vKKraw(:,:,ii)=abs(KKreconComplex).^0.25;  %frame for video
vKKmin(:,:,ii)=kkmin.^0.25;                %frame for video


%calculate histogram of kdiff
dd=0;
for xx=1:numTXangles
    for yy=1:numRXangles
        dd(xx,yy)=TXangle(xx)-RXangle(yy);
    end
end
figure(110)
histogram(dd(:),numTXangles*numRXangles);

end

dipshow3D(34,vKKraw)  %raw images used for min processing
dipshow3D(35,vKKmin)  %progressive result of min processing

%stack2avi(vKKmin,'testvideo',4); 


