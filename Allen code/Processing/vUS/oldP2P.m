%% Description:
%   Convert Jianbo's old-format acquisition/DAQ parameters struct (the
%   "P" struct attached to old lab sIQ/IQ data, with fields like CCFR,
%   wavelength, vSound, etc.) into Allen's current P struct format (the
%   struct saved in params.mat alongside acquired IQ data, with fields
%   like frameRate, wl, etc.), for use with vUS_2D.m and other
%   Allen code\Processing scripts.
%
%   This is the inverse of P_renamer.m (Jianbo code understanding\Functions),
%   which converts Allen's P struct into Jianbo's old format. Only the
%   fields P_renamer.m derives directly from P_acq are reversed here;
%   fields P_renamer.m computes from other fields (dxImg, xCoor, dzImg,
%   zCoor, numAcqs, t2NextSupPlaneDAQ, nSmplPerWvlnth) are redundant with
%   fields already covered below and are not reconstructed separately.
%
% Inputs:
%   P_old: Jianbo's old-format acquisition struct, with (at least):
%       P_old.CCFR, numCCframes, CCangle, numAngles, startAngle, dAngle,
%       TWfrequency, TWnHC, SampleMode, nSperWave, t2NextPlaneDAQ,
%       tIntPDI, nCh, vSound, LensDelay, frequency, pitch, actZsamples,
%       PeakDelay, TXDelay, startDepth, startDepthMM, endDepth,
%       endDepthMM, maxDepth, wavelength
%
% Outputs:
%   P: Allen's acquisition parameters struct (as normally loaded from
%       params.mat), with fields frameRate, numFramesPerBuffer, maxAngle,
%       angles, na, wl, startDepth, startDepthMM, endDepth, endDepthMM,
%       maxAcqLength_adjusted, samplesPerWave, Resource, Trans, Receive,
%       TW, SeqControl, TX

function [P] = oldP2P(P_old)
    P.frameRate = P_old.CCFR;
    P.numFramesPerBuffer = P_old.numCCframes;
    P.maxAngle = P_old.CCangle;

    P.na = P_old.numAngles;
    % dAngle only stores abs(angles(2) - angles(1)) in P_old, so the sign
    % is not recoverable; angles are reconstructed assuming ascending order
    P.angles = P_old.startAngle + (0:P_old.numAngles - 1) * P_old.dAngle;

    % P_old.wavelength is NOT used here: its units are inconsistent across
    % Jianbo's pipeline versions (mm in the original acquisition scripts'
    % vSound/frequency*1e-3, but a raw meters passthrough in P_renamer.m's
    % reverse-compatibility path). Instead, wl is derived directly from
    % vSound and frequency (Trans.frequency, conventionally in MHz for
    % Verasonics), the same way PRSSinfo.f0 is built from a MHz input in
    % MAIN_vUS_invivo.m ("f0 = frequency_MHz * 1e6"). This avoids a
    % possible silent 1000x unit error.
    P.wl = P_old.vSound / (P_old.frequency * 1e6); % m

    P.startDepth = P_old.startDepth;
    P.startDepthMM = P_old.startDepthMM;
    P.endDepth = P_old.endDepth;
    P.endDepthMM = P_old.endDepthMM;
    P.maxAcqLength_adjusted = P_old.maxDepth;

    P.samplesPerWave = P_old.nSperWave;

    P.Resource.Parameters.numRcvChannels = P_old.nCh;
    P.Resource.Parameters.speedOfSound = P_old.vSound;

    P.Trans.lensCorrection = P_old.LensDelay;
    P.Trans.frequency = P_old.frequency;
    P.Trans.spacingMm = P_old.pitch;

    P.Receive(1).sampleMode = P_old.SampleMode;
    P.Receive(1).endSample = P_old.actZsamples;

    P.TW.peak = P_old.PeakDelay;
    P.TW.Parameters = [P_old.TWfrequency, NaN, P_old.TWnHC]; % Parameters(2) is not recoverable from P_old

    P.SeqControl(1).argument = P_old.t2NextPlaneDAQ;
    P.SeqControl(4).argument = P_old.tIntPDI;

    for k = 1:size(P_old.TXDelay, 1)
        P.TX(k).Delay = P_old.TXDelay(k, :);
    end
end
