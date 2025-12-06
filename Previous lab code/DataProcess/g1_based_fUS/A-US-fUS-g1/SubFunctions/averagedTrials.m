function Results = averagedTrials(roiGGV, trial)

trialroiGGV = reshape(medfilt1(roiGGV(trial.nBase-trial.nRest+1: end-trial.nRest),3,'truncate'), [trial.nlength,trial.n]);
ratiotrialroiGGV = trialroiGGV./mean(trialroiGGV(1:trial.nRest,:))*100;
mratiotrialroiGGV = mean(ratiotrialroiGGV,2);
stdratiotrialroiGGV = std(ratiotrialroiGGV,1,2);
semratiotrialroiGGV = std(ratiotrialroiGGV,1,2)./sqrt(trial.n);

Results.data = roiGGV;
Results.trial = trialroiGGV;
Results.ratio = ratiotrialroiGGV;
Results.m = mratiotrialroiGGV;
Results.std = stdratiotrialroiGGV;
Results.sem = semratiotrialroiGGV;

end