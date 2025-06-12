function vtpXX = calTime2Peak(tcvein, trial)
% Bingxue Liu, 02/04/2022    
window = [trial.nRest+1, trial.nRest+trial.nStim+2]; 
ninterp = 10;
interpx = interp1(1:trial.nlength,tcvein, 1:1/ninterp:trial.nlength, "spline");
tcvein_bar = rescale(interpx,'InputMin',mean(interpx(1:trial.nRest*ninterp)),'InputMax',max(interpx));
[v0,vtpXX0] = max(tcvein_bar(1:window(2)*ninterp));
[~,vtpXX1] = min(abs(tcvein_bar(1:vtpXX0)-v0*0.9)); % time to peak;
[~,vtpXX2] = min(abs(tcvein_bar(1:vtpXX0)-v0*0.3)); % onset time;
vtpXX = [vtpXX1;vtpXX2]/ninterp-trial.nRest;