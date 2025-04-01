function trials = organizeTrials(Daq_Frequency,trial_length,num_trials,Data)
    for i = 0:(num_trials -1)
        startpoint = (Daq_Frequency*trial_length*i) + 1;%Start of each trials
        trials(i+1).ts = Data.ts(startpoint:startpoint +(Daq_Frequency*trial_length)-1);
        trials(i+1).frameReadOut = Data.frameReadOut(startpoint:startpoint +(Daq_Frequency*trial_length)-1);
        trials(i+1).stimulusTrigger = Data.stimulusTrigger(startpoint:startpoint +(Daq_Frequency*trial_length)-1);
        trials(i+1).ledTrigger = Data.ledTrigger(startpoint:startpoint +(Daq_Frequency*trial_length)-1);
        trials(i+1).baslerExposure = Data.baslerExposure(startpoint:startpoint +(Daq_Frequency*trial_length)-1);

        trials(i+1).ts = trials(i+1).ts + i.*trial_length;
        trials(i+1).speaker = Data.speaker(startpoint:startpoint +(Daq_Frequency*trial_length)-1);
        trials(i+1).lever = Data.lever(startpoint:startpoint +(Daq_Frequency*trial_length)-1);
    end
end