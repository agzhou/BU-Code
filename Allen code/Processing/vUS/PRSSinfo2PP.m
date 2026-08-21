%% Description:
%   Convert Jianbo's PRSSinfo struct (the processing-parameters struct
%   passed into sIQ2vUS_NPDV.m) into the PP struct used by vUS_2D.m and
%   its helper functions (stackData, unstackData, spectralSNR,
%   g1BasedSNR, pnSpectralSNR, findVzPhaseDiff, InitvUS2DParamsWithMesh,
%   g1vUS2D_Jac).
%
% Inputs:
%   PRSSinfo: Jianbo's processing-parameters struct, with (at least):
%       PRSSinfo.rFrame: sIQ frame rate, Hz
%       PRSSinfo.f0: Transducer center frequency, Hz
%       PRSSinfo.C: Sound speed in the sample, m/s
%       PRSSinfo.g1nTau: maximum number of time lags
%   sIQ: [z voxels, x voxels, frames] complex data matrix that PRSSinfo
%       describes (used to get zp, xp, nf, matching what
%       sIQ2vUS_NPDV.m itself derives via size(sIQ))
%
% Outputs:
%   PP: Allen's 2D vUS processing-parameters struct, with fields zp, xp,
%       nf, nTau, xDim, zDim, fDim, dimensionality, faxis, freqMask,
%       frameRate, wl, k0 (see createStruct.m call in vUS_2D.m)

function [PP] = PRSSinfo2PP(PRSSinfo, sIQ)
    [zp, xp, nf] = size(sIQ);

    PP.zp = zp;
    PP.xp = xp;
    PP.nf = nf;
    PP.nTau = PRSSinfo.g1nTau;

    PP.xDim = 2; % Dimension of the data corresponding to x (lateral direction)
    PP.zDim = 1; % Dimension of the data corresponding to z (axial direction)
    PP.fDim = 3; % Dimension of the data corresponding to frequency (or time)
    PP.dimensionality = 2; % 2D data

    PP.faxis = linspace(-PRSSinfo.rFrame/2, PRSSinfo.rFrame/2, nf)';
    PP.freqMask = abs(PP.faxis) > 1100; % [Hz] -- matches the fixed 1100 Hz system-noise band masked out in both vUS_2D.m and sIQ2vUS_NPDV.m (fRangeSignal)

    PP.frameRate = PRSSinfo.rFrame;
    PP.wl = PRSSinfo.C / PRSSinfo.f0; % lambda0 = C / f0, as computed at the top of sIQ2vUS_NPDV.m
    PP.k0 = 2*pi / PP.wl;
end
