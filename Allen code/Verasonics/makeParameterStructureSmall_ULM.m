
Psmall.endDepthMM = endDepthMM;
Psmall.fps_target = fps_target;
Psmall.PRF = PRF;
Psmall.maxAngle = maxAngle;
Psmall.na = na;
Psmall.numFramesPerBuffer = numFramesPerBuffer;

Psmall.bufferIndex = 0;
Psmall.frameRate = frameRate;

if exist('bufferDutyCycle', 'var')
    Psmall.bufferDutyCycle = bufferDutyCycle;
end