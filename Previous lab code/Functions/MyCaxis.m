function cAxis=MyCaxis(imgData,LowThre, HighThre)
if nargin <2
    LowThre=0.08;
    HighThre=0.1;
elseif nargin<3
    HighThre=0.1;
end
cAxis=[min(imgData(:))*(1+sign(min(imgData(:)))*LowThre), max(imgData(:))*(1-sign(max(imgData(:)))*HighThre)];