%%
dipstart

%% Read raw data 

% FileNameBase='D:\OBUS_GE9L\';
% FileNameBase='U:\eng_research_biomicroscopy\Projects\Ultrasound\Datasets\Jerome Data\';
FileNameBase='U:\eng_research_biomicroscopy\Projects\Ultrasound\Datasets\Plexiglass Data\';

%VMat='CIRS3';
% VMat='abdomen4';
VMat='PlexiPhantom';

numTXangles=15;

[RawData,TXangle,TXdelay,TXapod]=RXunscrambleN(append(FileNameBase,VMat),numTXangles);  %raw data

Pitch = 0.00023;    % m
SamplingFrequency = 4*(250/48)*10^6;   % Hz
c = 1540;   % m/s
NA=0.2;     % numerical aperture
s=2*Pitch*SamplingFrequency/c;  %aspect ratio
xSize=192;

zBegin=1000; 
zSize=1000; 

RawDataComplex=DataHilbert(RawData);

%% Conventional DAS
% [Recon,ReconSum]=BeamformAnglesComplex(RawData,zBegin,zSize,TXangle,NA,s);  %beamformed data
% aReconSum=abs(ReconSum);
% 
% dipshow3D(11,(aReconSum).^0.25)


%% KK-DAS 
%%

%% Define RXangles and compress KK-DAS data
rangeRXangles =12;    %half total range in degrees  (TX range = +/-12 degrees)
numRXangles =19 ;      %number of angle steps  (TX num = 15)
 
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

RawDataKKComplex=DataCompressKK(RawDataComplex,RXangle,s); 

%% Determine zones for SoS estimation
numZoneZ=5;
widthZoneZ=1; %gaps between zones if width<1; zones overlap is width>1
numZoneX=1;
widthZoneX=1;

% Define zones in Z
zoneSizeZ=floor(zSize/numZoneZ*widthZoneZ)-1;
if widthZoneZ<=1
    zoneSeparationZ=floor(zSize/numZoneZ);
    zoneCenterStartZ=floor(zoneSeparationZ/2)+zBegin;
else
    zoneSeparationZ=floor((zSize-zoneSizeZ)/(numZoneZ-1));
    zoneCenterStartZ=ceil(zoneSizeZ/2)+zBegin;
end
for n=1:numZoneZ
    zoneCenterZ(n)=zoneCenterStartZ+(n-1)*zoneSeparationZ;
    zoneBeginZ(n)=zoneCenterZ(n)-round(zoneSizeZ/2);
end

% Define zones in X (this part isn't useful/needed because of my beamformer)
zoneSizeX=floor(xSize/numZoneX*widthZoneX)-1;  
if widthZoneX<=1
    zoneSeparationX=floor(xSize/numZoneX);
    zoneCenterStartX=floor(zoneSeparationX/2)+1;
else
    zoneSeparationX=floor((xSize-zoneSizeX)/(numZoneX-1));
    zoneCenterStartX=ceil(zoneSizeX/2)+1;
end
for n=1:numZoneX
    zoneCenterX(n)=zoneCenterStartX+(n-1)*zoneSeparationX;
    zoneBeginX(n)=zoneCenterX(n)-round(zoneSizeX/2);
end
    
if numZoneX<3, zoneCenterX=[46,96,146]; end %Needed becasue of my beamformer

[meshZoneX,meshZoneZ]=meshgrid(zoneCenterX,zoneCenterZ);

%% Cycle through zones and finbd SoS adjust in each zone
%Continues to step while metric increases.
%When metric decreases, does parabolic fit using last three steps  

cStep=0.01;  %search step size (i.e. relative change in SoS)
numIterate=10;  %max number of search iterations
a=2;   %power law defining metric (2=intensity; 4 = variance
cEst=zeros(numZoneZ,numZoneX);

c0=0.98; % starting adjust for SoS search. Normally c0=1 (i.e. no adjust)
for ZoneZ=1:numZoneZ
    ZoneZ
    if ZoneZ>1, c0=cEst(ZoneZ-1,2); end %Use SoS from previous zone to start estimate for next zone
   
eConstructComplex=BeamformKK(RawDataKKComplex,zoneBeginZ(ZoneZ),zoneSizeZ,xSize,TXangle,RXangle,s/c0);
KKrecon=sum(eConstructComplex,[3,4]);
%metric0=squeeze(sum(abs(KKrecon).^a,'all')); %initial point 1
metric0=sum(abs(KKrecon).^4,'all')/sum(abs(KKrecon).^2,'all')^2;

c1=c0+cStep;
eConstructComplex=BeamformKK(RawDataKKComplex,zoneBeginZ(ZoneZ),zoneSizeZ,xSize,TXangle,RXangle,s/c1);
KKrecon=sum(eConstructComplex,[3,4]);
%metric1=squeeze(sum(abs(KKrecon).^a,'all')); %intial point 2
metric1=sum(abs(KKrecon).^4,'all')/sum(abs(KKrecon).^2,'all')^2;

% Choose direction for hill climbing
if metric1>metric0
    cVec=[c1,c0,c0];
    metricVec=[metric1,metric0,metric0];
else
    cStep=-cStep;
    c2=c0+cStep;
    eConstructComplex=BeamformKK(RawDataKKComplex,zoneBeginZ(ZoneZ),zoneSizeZ,xSize,TXangle,RXangle,s/c2);
    KKrecon=sum(eConstructComplex,[3,4]);
    %metric2=squeeze(sum(abs(KKrecon).^a,'all'));
    metric2=sum(abs(KKrecon).^4,'all')/sum(abs(KKrecon).^2,'all')^2;
    cVec=[c2,c0,c1];
    metricVec=[metric2,metric0,metric1];
end

%Perfom hill climbing until metric starts decreasing. When it does, do parabolic fit
%using last 3 steps
for nIterate=1:numIterate
    if metricVec(1)>metricVec(2)
        cVec=circshift(cVec,1);
        metricVec=circshift(metricVec,1);
        cVec(1)=cVec(2)+cStep;
        eConstructComplex=BeamformKK(RawDataKKComplex,zoneBeginZ(ZoneZ),zoneSizeZ,xSize,TXangle,RXangle,s/cVec(1));
        KKrecon=sum(eConstructComplex,[3,4]);
        %metricVec(1)=squeeze(sum(abs(KKrecon).^a,'all'));
        metricVec(1)=sum(abs(KKrecon).^4,'all')/sum(abs(KKrecon).^2,'all')^2;      
    else
        break;
    end
        
end
parabolaMat=[cVec.^2;cVec;[1,1,1]]'; 
parabolaVec=pinv(parabolaMat)*metricVec'; %vector defining parabolic fit
cEst(ZoneZ,2)=-parabolaVec(2)/2/parabolaVec(1); %vertex of parabola, correpsoing to location of max

cEst(ZoneZ,1)=cEst(ZoneZ,2); %needed for my beamformer
cEst(ZoneZ,3)=cEst(ZoneZ,2); %needed for my beamformer
 
end
cEst

%% Interpolate cRatio (fill in gaps between zones)
[meshInterpX,meshInterpZ]=meshgrid(1:xSize,zBegin:zBegin+zSize-1);
cEstInterp=interp2(meshZoneX,meshZoneZ,cEst,meshInterpX,meshInterpZ,'spline');
cAdjust=cEstInterp;
dipshow3D(100,cAdjust)

%% Check result without and with SoS adjust
eConstructComplex=BeamformKK(RawDataKKComplex,zBegin,zSize,xSize,TXangle,RXangle,s);
KKrecon=squeeze(sum(eConstructComplex,[3,4]));

eConstructComplex=BeamformKKadjust(RawDataKKComplex,zBegin,zSize,xSize,TXangle,RXangle,s,cAdjust);
KKreconAdjust=sum(eConstructComplex,[3,4]);

dipshow3D(200,abs(KKrecon).^0.25)
dipshow3D(201,abs(KKreconAdjust).^0.25)

%% Below not needed
%%
zBegin1=1600;
zSize1=200;
kkMax=15;
DeltaErr=0.01;
KKreconI=0;
for kk=1:kkMax
    err(kk)=1+(kk-1-(kkMax-1)/2)*DeltaErr;
eConstructComplex=BeamformKK(RawDataKKComplex,zBegin1,zSize1,xSize,TXangle,RXangle,s/err(kk));
KKrecon=sum(eConstructComplex,[3,4]);
KKreconI(kk)=sum(abs(KKrecon).^2,'all');
KKreconV(kk)=sum(abs(KKrecon).^4,'all');
end

KKreconC=KKreconV./KKreconI.^2;

[mm,ee]=max(KKreconC,[],'all','linear');
err(ee)

figure(18)
plot(KKreconI)
figure(19)
plot(KKreconV)
figure(20)
plot(KKreconC)

eConstructComplex=BeamformKK(RawDataKKComplex,zBegin,zSize,xSize,TXangle,RXangle,s);
KKrecon=sum(eConstructComplex,[3,4]);
dipshow3D(21,abs(KKrecon).^0.25)

eConstructComplex=BeamformKK(RawDataKKComplex,zBegin,zSize,xSize,TXangle,RXangle,s/err(ee));
KKrecon=sum(eConstructComplex,[3,4]);
dipshow3D(22,abs(KKrecon).^0.25)


%%
vTest=zeros(1,zSize);
vTest(zoneCenterZ)=cEst(:,1);
%vTest=circshift(vTest,-zBegin);
vertexInterp=interp1(zoneCenterZ,cEst(:,1),0.5:(zSize-0.5),'linear');
%cRatio=interp1(zoneCenter,vertex,(zBegin+0.5):(zBegin+zSize-0.5),'pchip');

figure(5)
hold off
plot(vTest)
hold on
plot(vertexInterp)
