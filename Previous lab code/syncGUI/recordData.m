function recordData(t,D)
    global Data
    t;
    Data.ts = [Data.ts; t];
    Data.frameReadOut = [Data.frameReadOut; D(:,1)];
    Data.stimulusTrigger = [Data.stimulusTrigger; D(:,2)];
    Data.ledTrigger = [Data.ledTrigger; D(:,3)];
    Data.speaker = [Data.speaker; D(:,4)];
    Data.lever = [Data.lever; D(:,5)];
    Data.baslerExposure = [Data.baslerExposure; D(:,6)];
end
