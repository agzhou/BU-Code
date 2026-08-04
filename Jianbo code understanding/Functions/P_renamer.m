function [P_post] = P_renamer(P_acq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%P_renamer reformats ultrasound acquisition parameters from new acquisition
%code (1/2025) for compatibility with existing ULM postprocessing code
%P_acq - acquisition parameters
%P_post - acquisition parameters, old format (compatible with
%SCC_JOBGEN_US_GENERATOR GUI)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P.PDIFR = 0.5;
P.numSupFrames = 1;
P.CCFR = P_acq.frameRate;
P.numCCframes = P_acq.numFramesPerBuffer;
P.CCangle = P_acq.maxAngle;
ANGLES = P_acq.angles;
P.numAngles = length(ANGLES);
TW_Parameters = P_acq.TW.Parameters;
P.TWfrequency = TW_Parameters(1);
P.TWnHC = TW_Parameters(3);
P.SampleMode = P_acq.Receive(1).sampleMode;
P.nSperWave = P_acq.samplesPerWave;
P.t2NextPlaneDAQ = P_acq.SeqControl(1).argument;
P.tIntPDI = P_acq.SeqControl(4).argument;
P.t2NextSupPlaneDAQ = P_acq.SeqControl(4).argument;
P.numAcqs = (P_acq.numFramesPerBuffer)*length(ANGLES);
P.nCh = P_acq.Resource.Parameters.numRcvChannels;
P.vSound = P_acq.Resource.Parameters.speedOfSound;
P.LensDelay = P_acq.Trans.lensCorrection;
P.frequency = P_acq.Trans.frequency;
P.pitch = P_acq.Trans.spacingMm;
P.actZsamples = P_acq.Receive(1).endSample;
%P.MiniIntSupFrame
P.PeakDelay = P_acq.TW.peak;
P.startAngle = ANGLES(1);
for k = 1:5
P.TXDelay(k,:) = P_acq.TX(k).Delay;
end
P.dAngle = abs(ANGLES(2)-ANGLES(1));
P.nSmplPerWvlnth = P_acq.samplesPerWave;
P.startDepth = P_acq.startDepth;
P.startDepthMM = P_acq.startDepthMM;
P.endDepth = P_acq.endDepth;
P.endDepthMM = P_acq.endDepthMM;
P.maxDepth = P_acq.maxAcqLength_adjusted;
P.wavelength = P_acq.wl;
P.dxImg = 1000*0.5*P_acq.Trans.spacingMm; %um
P.xCoor = [0 : P.dxImg : (P_acq.Resource.Parameters.numRcvChannels)*(P_acq.Trans.spacingMm)]/1000;
P.dzImg = 1000*0.5*P_acq.wl; %um
P.zCoor = [(P.startDepth)*(P_acq.wl) : P.dzImg : (P.endDepth)*(P_acq.wl)];
%
P_post = P;
end