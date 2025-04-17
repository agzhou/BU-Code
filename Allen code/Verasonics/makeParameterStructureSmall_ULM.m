
Psmall.endDepthMM = endDepthMM;
% Psmall.fps_target = fps_target;
Psmall.PRF = PRF;
Psmall.maxAngle = maxAngle;
Psmall.na = na;

if exist('numFramesPerBuffer', 'var')
    Psmall.numFramesPerBuffer = numFramesPerBuffer;
else
    Psmall.numFramesPerBuffer = numFramesPerSF;
end

Psmall.bufferIndex = 0;
Psmall.frameRate = frameRate;

if exist('bufferDutyCycle', 'var')
    Psmall.bufferDutyCycle = bufferDutyCycle;
end