%% function for processing IQ data to V, vUS
% cluter rejection is based on singular value decomposition (SVD)
function SCC_IQ2GG(datapath, filename)
% IQ: IQ data
% PRMT: data processing parameter
% vUS: Obtained vUS results
load([datapath,'vUS-PRSinfo.mat'])
ROI=PRSinfo.ROI;
%% SVD process 1 (direct SVD use MATLAB)
disp('Loading data ...');
load ([datapath, filename]);
disp('Data loaded!');
%%
disp('SVD Processing ...');
cIQ=(IQ(min(ROI(:,1)):max(ROI(:,1)),min(ROI(:,2)):max(ROI(:,2)),:));
% [sIQ, Noise]=SVDfilter(cIQ,PRSinfo.SignalRank); % sIQ: signal IQ
[sIQ, sIQHP, sIQHHP, eqNoise]=IQ2sIQ(IQ,PRSSinfo);
        [PDI0]=sIQ2PDI(sIQ);  % PDI processing
        [PDIHP0]=sIQ2PDI(sIQHP); % PDI processing
        [PDIHHP0]=sIQ2PDI(sIQHHP);  % PDI processing
disp('sIQ to GG ...');
GG=IQ2g1(sIQHP,PRSinfo.g1startT,PRSinfo.g1nT,PRSinfo.g1nTau);
disp('vUS Processing ...');

%% 
SavePath=['/',strjoin(pathInfo(1:end-2),'/'),'/RESULT-',pathInfo{end-1},'-GG/'];
% if ~exist(SavePath)
    mkdir(SavePath);
% end

save([SavePath,'GG',filename(3:end)],'-v7.3','GG','PDI','PDIHP','PDIHHP','eqNoise','PRSSinfo','P');
disp('vUS data saved!')
