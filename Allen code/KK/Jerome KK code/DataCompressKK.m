function RawDataKK=DataCompressKK(Data,RXangle,s)
% Compresses RawData(t,x,TXangle) -> RawDataKK(t,TXangle,RXangle)  
% s is "aspect ratio"
% From TestKKadaptive.m: s=2*Pitch*SamplingFrequency/c;  %aspect ratio

tSize=size(Data,1);
xSize=size(Data,2);
TXSize=size(Data,3);
RXSize=numel(RXangle);

RawDataKK=zeros(tSize,TXSize,RXSize);
DataTemp=zeros(tSize,xSize);

for nR=1:RXSize
    
    slope=s*sin(RXangle(nR))/2;

for nT=1:TXSize
       
       for nx=1:xSize
           if slope>0
               nShift=((nx-1)*abs(slope));
           else
               nShift=((xSize-nx)*abs(slope));
           end
           
           DataTemp(:,nx)=circshift(Data(:,nx,nT),round(nShift));
       end 
       
       RawDataKK(:,nT,nR)=sum(DataTemp,2);
        
end

end

