% 1D move window average
function MWAdata=movavg1D(data,MW_length)
if nargin<2
    MW_length=10;
end

Length_data=length(data);
for in=1:Length_data
    MWAdata(in)=mean(data(in:min(in+MW_length,Length_data)));
end