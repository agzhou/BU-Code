% 1D move window average
function MWAdata=MoveAvg(data,MW_length,DIM)
if nargin<2
    MW_length=10;
    DIM=1;
end
if nargin<3
    DIM=1;
end

[D1,D2,D3]=size(data);
if DIM==1
    for in=1:D1
        MWAdata(in,:,:)=squeeze(mean(data(in:min(in+MW_length,D1),:,:),1));
    end
elseif DIM==2
    for in=1:D2
        MWAdata(:,in,:)=squeeze(mean(data(:,in:min(in+MW_length,D2),:),2));
    end
elseif DIM==3
    for in=1:D3
        MWAdata(:,:,in)=squeeze(mean(data(:,:,in:min(in+MW_length,D3)),3));
    end
end