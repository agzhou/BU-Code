P.angles = angles;
P.endDepth = endDepth;
P.endDepthMM = endDepthMM;
P.Event = Event;
P.PRF = PRF;
P.frameRate = frameRate;
% P.flow_v_mm_s = flow_v_mm_s;
P.maxAcqLength_adjusted = maxAcqLength_adjusted;
P.maxAngle = maxAngle;
P.Media = Media;
P.na = na;
% P.nf = nf;
P.numFramesPerBuffer = numFramesPerBuffer;

P.numElements = numElements;
P.Receive = Receive;
P.Resource = Resource;
P.SeqControl = SeqControl;
P.startDepth = startDepth;
P.startDepthMM = startDepthMM;
P.TGC = TGC;
P.Trans = Trans;
P.TW = TW;
P.TX = TX;
P.wl = wl;
% P.L = L;

P.TPC = TPC;
P.samplesPerWave = Receive(1).samplesPerWave;

if exist('bufferDutyCycle', 'var')
    P.bufferDutyCycle = bufferDutyCycle;
end

% Functional params
P.apis = apis;
P.vts = vts;
P.numTrials = numTrials;
P.daqrate = P.daqrate;
P.Mcr_d = Mcr_d;
P.Mcr_fcp = Mcr_fcp;