%%
% dipstart

%% Read raw data 

% FileNameBase='D:\OBUS_GE9L\';

currentDir = matlab.desktop.editor.getActiveFilename; 
currentDir = regexp(currentDir, filesep, 'split');
dataFilePath = fullfile(currentDir{1:find(contains(currentDir,"Ultrasound"),1)},"Datasets\");

%VMat='Sim15Cos2_Point';
%VMat='Sim15_Point';
% VMat='CIRS1';
VMat='abdomen4';

numTXangles=15;

% Jerome Version
% [RawData,TXangle,TXdelay,TXapod]=RXunscrambleN(append(FileNameBase,VMat),numTXangles);  %raw data

% Google Drive version
[RawData,TXangle,TXdelay,TXapod]=RXunscrambleN(append(dataFilePath,'Jerome Data\',VMat),numTXangles);  %raw data

Pitch = 0.00023;    % m
SamplingFrequency = 4*(250/48)*10^6;   % Hz
c = 1540;   % m/s
NA=0.6;     % numerical aperture
s=2*Pitch*SamplingFrequency/c;  %aspect ratio
xSize=192;

zBegin=1;
zSize=1500;

%% Conventional DAS 

[Recon,ReconSum]=BeamformAngles(RawData,zBegin,zSize,TXangle,NA,s);  %beamformed data
[cRecon,aRecon,pRecon]=DataHilbert(Recon);  %analytic representation of beamformed data (complex)
cReconSum=sum(cRecon,3);
aReconSum=abs(cReconSum);

%% display routines with dipimage (scroll with 'n' and 'p')
% dipshow3D(10,RawData)
% dipshow3D(11,(aRecon).^0.25)
% dipshow3D(12,(aReconSum).^0.25)

figure
plotGammaScaleImage(aReconSum,0.25)
% axis image


%% KK-DAS 

% Step 1: Define RXangle for KK-DAS data compression
rangeRXangles = 10;    %half total range in degrees  (TX range = +/-12 degrees)
numRXangles = 10;      %number of angle steps  (TX num = 15)
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

% Step 2: Compress data
RawDataKK=DataCompressKK(RawData,RXangle,s);

% Step 3: Beamform
[eConstruct,eConstructSum]=BeamformKK(RawDataKK,zBegin,zSize,xSize,TXangle,RXangle,s);
KKrecon=DataHilbert(eConstructSum);

%% Dipimage display functions
% dipshow3D(20,RawData)
% dipshow3D(21,RawDataKK)
% dipshow3D(22,eConstructSum)
% dipshow3D(23,eConstruct)
% dipshow3D(25,abs(KKrecon).^0.25);

figure
plotGammaScaleImage(KKrecon,0.25)
% axis image

%% Below not needed
%% Partial KK-DAS (Scan RX angle only (TX fixed))

TXindex=8;
RawDataKK=DataCompressKK(RawData(:,:,TXindex),RXangle,s);
[eConstruct,eConstructSum]=BeamformKK(RawDataKK,zBegin,zSize,xSize,TXangle(TXindex),RXangle,s);
KKrecon=DataHilbert(eConstructSum);

dipshow3D(30,RawData(:,:,TXindex))
dipshow3D(31,RawDataKK)
dipshow3D(32,eConstructSum)
dipshow3D(33,eConstruct)
dipshow3D(35,abs(KKrecon).^0.25);


%% Partial KK-DAS (Scan TX angle only (RX fixed))

RXangle0=0*pi/180;
RawDataKK=DataCompressKK(RawData,RXangle0,s);
[eConstruct,eConstructSum]=BeamformKK(RawDataKK,zBegin,zSize,xSize,TXangle,RXangle0,s);
KKrecon=DataHilbert(eConstructSum);

dipshow3D(40,RawData)
dipshow3D(41,RawDataKK)
dipshow3D(42,eConstructSum)
dipshow3D(43,eConstruct)
dipshow3D(45,abs(KKrecon).^0.25);


%% Partial KK-DAS (Both TX and RX fixed)

TXindex=1;
RXindex=1;
RawDataKK=DataCompressKK(RawData(:,:,TXindex),RXangle(RXindex),s);
[eConstruct,eConstructSum]=BeamformKK(RawDataKK,zBegin,zSize,xSize,TXangle(TXindex),RXangle(RXindex),s);
KKrecon=DataHilbert(eConstructSum);

dipshow3D(50,RawData(:,:,TXindex))
dipshow3D(51,RawDataKK)
dipshow3D(52,eConstructSum)
dipshow3D(53,eConstruct)
dipshow3D(55,abs(KKrecon).^0.25);


