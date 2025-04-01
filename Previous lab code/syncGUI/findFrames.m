function out = findFrames(trials)
startFrame = 0;
for numtrial = 1:length(trials)
    numtrial;
    % Get Data
    Data = trials(numtrial);
    %Threshhold Stimulus input
    stimOn = find(Data.stimulusTrigger>3);
    %Threshold Frame input
    prevVal = zeros(length(Data.frameReadOut),1);
    prevVal(2:end) = Data.frameReadOut(1:end-1);
    frameDiff = Data.frameReadOut - prevVal;
    frameDiff(1) = 0;
    frameOut = find(frameDiff>1.4);
    Iframe = find(diff(Data.ts(frameOut))<= .001);
    for i = 1:length(Iframe)
        frameOut = [frameOut(1:Iframe(i)-1);frameOut(Iframe(i)+1:end)];
        disp("removed Frame")
    end
    min(diff(Data.ts(frameOut)));
    max(diff(Data.ts(frameOut)));
    % Threshhold LED input for On and Off
    ledDiff = diff(Data.ledTrigger);
    ledDiff = [Data.ledTrigger(1);ledDiff];
%     prevled = zeros(length(Data.ledTrigger),1);
%     prevled(3:end) = Data.ledTrigger(1:end-2);
%     ledDiff2 = Data.ledTrigger - prevled;
    %ledDiff = (ledDiff+ledDiff2)./2;
    ledOn = [];
    ledOut = find(ledDiff>1.3);
    ledStop = find(ledDiff<-1.3);
    [minOut,Iout] = min(diff(Data.ts(ledOut)));
    [minStop,Istop] = min(diff(Data.ts(ledStop)));
    if (minOut < 0.01)
        ledOut = [ledOut(1:Iout);ledOut(Iout+2:end)];
    end
    if (minStop <0.01)
        ledStop = [ledStop(1:Istop);ledStop(Istop+2:end)];
    end
    if length(ledOut)>length(ledStop)
        ledStop = [trials(numtrial).ts(ledStop);trials(numtrial).ts(end)];
        ledOn = Data.ts(ledOut);
    elseif length(ledOut)<length(ledStop)
        ledOn = [Data.ts(1);trials(numtrial).ts(ledOut)];
        ledStop = trials(numtrial).ts(ledStop);
    else
        ledOn = Data.ts(ledOut);
        ledStop = trials(numtrial).ts(ledStop);
    end
    trials(numtrial).ledOn = ledOn;
    trials(numtrial).ledOff = ledStop;
    length(ledOut);
    %Getting time stamps
    stimTime = Data.ts(stimOn);
    frameTime = Data.ts(frameOut);
    ledTime = ledOn;
%     max(diff(ledTime))
%     min(diff(ledTime))
%     max(diff(ledStop))
%     min(diff(ledStop))
    max(ledStop-ledOn);
    min(ledStop-ledOn);
    onFrames = [];
    LEDs = [];
    trials(numtrial).topLed = [];
    trials(numtrial).topTime = [];
    % Determine place of each frame in LED cycle
    for j = 1:length(frameOut)
        led = 0;
        currTime = frameTime(j);
        topLed = find(currTime<=ledStop,1);
        topTime = ledStop(topLed);
        trials(numtrial).topLed = [trials(numtrial).topLed;topLed];
        trials(numtrial).topTime = [trials(numtrial).topTime;topTime];
        if isempty(topLed) 
            %Frame is after last falling edge
            if Data.ts(end) == ledStop(end)
                led = 1;
            else
                endDist = currTime - ledStop(end);
                if endDist <= 0.0333
                    led = 2;
                else
                    led = 3;
                end
            end
        elseif topLed > 1
            botTime = ledOn(topLed);
            %Frame falls during on trigger
            if currTime >= botTime
                led = 1;
            else
                timeDiff = botTime - currTime;
                %Width between end of prev led1 and start of curr led1
                BGdiv = (botTime - ledStop(topLed-1))./2;
                if timeDiff >= BGdiv
                    led = 2;
                else
                    led = 3;
                end
            end
        else
            if (ledOn(1) <= currTime)
                led = 1;
            else
                 endDist =  ledOn(1) - currTime;
                 if endDist> .0333
                     led = 2;
                 else
                     led =3;
                 end
            end
        end
%         if size(botFrame == 0)
%             trigDiff = ledTime(1) - currTime;
%         elseif (botFrame == length(ledTime))
%             trigDiff = currTime - botFrame;
%         else
%             trigDiff = ledTime(botFrame+1) - ledTime(botFrame);
%         end
       LEDs = [LEDs;led];
    end
    for i = 1:length(stimTime)
        currTime = stimTime(i);
        frame = find(frameTime == currTime);
            onFrames = [onFrames;frame];
    end
    trials(numtrial).startFrame = min(onFrames)+ startFrame;
    trials(numtrial).stopFrame = max(onFrames)+ startFrame;
    trials(numtrial).frameLED = horzcat(frameTime,LEDs);
    %LEDs(end-1:end);
    out = trials;
    startFrame = startFrame + length(frameOut);
end
end