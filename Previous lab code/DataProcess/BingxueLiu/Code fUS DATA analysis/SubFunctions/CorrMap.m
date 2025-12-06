function corrmap = CorrMap(data, fitdata)
% data: 3 dimention, [nx, ny, nt];
% fitdata: nt*1 vector
% corrmap: nx*ny 
% calculate Pearson Correlation Coefficient between spatialtemperal data and a temperal data;
% last modified: 06/02/2021, Bingxue Liu
fitdataN = (fitdata-mean(fitdata))./sqrt(mean((fitdata-mean(fitdata)).^2));
dataN = (data-mean(data,3))./sqrt(mean((data-mean(data,3)).^2, 3));
corrmap = mean(dataN.*permute(fitdataN,[1 3 2]),3);