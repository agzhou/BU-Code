function [eConstruct,eConstructSum,eDelay] = BeamformKKadjust(RawDataKK, tBegin, tSpan, xSize, TXangle, RXangle, s,cAdjust)

tSize = size(RawDataKK,1);
TXSize= size(RawDataKK,2);
RXSize= size(RawDataKK,3);
zSize=tSpan;
  
eConstruct=zeros(zSize,xSize,TXSize,RXSize);
eDelay=zeros(zSize,xSize,TXSize,RXSize); %added
eConstructSum=zeros(zSize,xSize);
zRXcos=zeros(tSpan,RXSize);
zTXcos=zeros(tSpan,TXSize);
xRXsin=zeros(xSize,RXSize);
xTXsin=zeros(xSize,TXSize);
xMin=ones(tSpan,RXSize);
xMax=ones(tSpan,RXSize)*xSize;

cc=cos(TXangle)*cos(RXangle');


for nR=1:RXSize
    zRXcos(:,nR)=(tBegin:tBegin+tSpan-1)*cos(RXangle(nR))/2;
    sinRX=abs(sin(RXangle(nR)));
    if RXangle(nR)>0
        xRXsin(:,nR)=s*sinRX*((1:xSize)-1)/2;
        xMin(:,nR)=1+floor((tBegin:tBegin+tSpan-1)*sinRX/s);
    else
        xRXsin(:,nR)=s*sinRX*((xSize:-1:1)-1)/2;
        xMax(:,nR)=xSize-floor((tBegin:tBegin+tSpan-1)*sinRX/s);
    end
end


for nT=1:TXSize
    zTXcos(:,nT)=(tBegin:tBegin+tSpan-1)*cos(TXangle(nT))/2;
    if TXangle(nT)>0
        xTXsin(:,nT)=s*abs(sin(TXangle(nT)))*((1:xSize)-1)/2;
    else
        xTXsin(:,nT)=s*abs(sin(TXangle(nT)))*((xSize:-1:1)-1)/2;
    end
end
    

for nT=1:TXSize
    for nR=1:RXSize
        %for nt=tBegin:tBegin+tSpan-1
        for nt=1:tSpan
            for nx=xMin(nt,nR):xMax(nt,nR)
                delay=zRXcos(nt,nR)+zTXcos(nt,nT)+(xRXsin(nx,nR)+xTXsin(nx,nT))./cAdjust(nt,nx);
                %eConstruct(nt,nx,nT,nR)=RawDataKK(round(delay),nT,nR);
                eConstruct(nt,nx,nT,nR)=RawDataKK(round(delay),nT,nR)*cc(nT,nR);
                %eDelay(nt,nx,nT,nR)=exp(i*2*pi*delay/4); %added
            end
        end
    end
end

eConstructSum=sum(eConstruct,[3,4]);





