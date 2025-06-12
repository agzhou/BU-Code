function corrmap = CorrMap(data, fitdata)
% data: 3 dimention, [nx, ny, nt];
% fitdata: nt*1 vector
% corrmap: nx*ny 
% calculate Pearson Correlation Coefficient between spatialtemperal data and a temperal data;
% last modified: 06/02/2021, Bingxue Liu
if numel(size(data))==3
fitdataN = (fitdata-mean(fitdata))./sqrt(sum((fitdata-mean(fitdata)).^2));
dataN = (data-mean(data,3))./sqrt(sum((data-mean(data,3)).^2, 3));
corrmap = sum(dataN.*permute(fitdataN,[1 3 2]),3);
else if numel(size(data))==2
        fitdataN = (fitdata-mean(fitdata))./sqrt(sum((fitdata-mean(fitdata)).^2)); 
        dataN = (data-mean(data))./sqrt(sum((data-mean(data)).^2));
        corrmap = sum(dataN.*fitdataN);
    end
end
