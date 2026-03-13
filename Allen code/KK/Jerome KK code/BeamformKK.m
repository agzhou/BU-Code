function [eConstruct,eConstructSum,delay] = BeamformKK(RawDataKK, tBegin, tSpan, xSize, TXangle, RXangle, s)

tSize = size(RawDataKK,1);
TXSize= size(RawDataKK,2);
RXSize= size(RawDataKK,3);
zSize=tSpan;
  
eConstruct=zeros(zSize,xSize,TXSize,RXSize);
eConstructSum=zeros(zSize,xSize);
zRXcos=zeros(tSpan,RXSize);
zTXcos=zeros(tSpan,TXSize);
xRXsin=zeros(xSize,RXSize);
xTXsin=zeros(xSize,TXSize);
xMin=ones(tSpan,RXSize);
xMax=ones(tSpan,RXSize)*xSize;


for nR=1:RXSize
    zRXcos(:,nR)=(tBegin:tBegin+tSpan-1)*cos(RXangle(nR))/2;
    sinRX=abs(sin(RXangle(nR)));
    if RXangle(nR)>0
        xRXsin(:,nR)=s*sinRX*((1:xSize)-0)/2;
        xMin(:,nR)=1+floor((tBegin:tBegin+tSpan-1)*sinRX/s);
    else
        xRXsin(:,nR)=s*sinRX*((xSize:-1:1)-0)/2;
        xMax(:,nR)=xSize-floor((tBegin:tBegin+tSpan-1)*sinRX/s);
    end
end


for nT=1:TXSize
    zTXcos(:,nT)=(tBegin:tBegin+tSpan-1)*cos(TXangle(nT))/2;
    if TXangle(nT)>0
        xTXsin(:,nT)=s*abs(sin(TXangle(nT)))*((1:xSize)-0)/2;
    else
        xTXsin(:,nT)=s*abs(sin(TXangle(nT)))*((xSize:-1:1)-0)/2;
    end
end
    

for nT=1:TXSize
    for nR=1:RXSize
        for nt=1:tSpan
            for nx=xMin(nt,nR):xMax(nt,nR)
                delay=zRXcos(nt,nR)+xRXsin(nx,nR)+zTXcos(nt,nT)+xTXsin(nx,nT);
                eConstruct(nt,nx,nT,nR)=RawDataKK(round(delay),nT,nR);
            end
        end
    end
end

eConstructSum=sum(eConstruct,[3,4]);

end




