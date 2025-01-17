%% This function generates a referrence index matrix containing the contribution sampling points index of all channels for every (x, y) coordinates
% iP: emission and receive information for the ith plane wave
    % iP: transducer parameters
    % iP.angle=P.startAngle+(i-1)*P.dAngle: plane wave emission angle
    % iP.vSound=P.vSound: speed of sound
    % iP.nCh=P.nCh: number of probe channels
    % iP.pitch=P.pitch: pitch between transducer elements, mm
    % iP.Wavelength=P.vSound*1e3/(P.frequency*1e6), mm
    % iP.dSample=iP.Wavelength/4, distance between samplling points, mm, the default is 4 sampling points per
    % iP.nZsample= P.actZsamples: number of samples for each pulse-echo acquisition
    % wavelength
    % iP.delay=P.startDepth*iP.Wavelength, delay for the first acquie point
% NA=rAperture/zDistance: reconstruction numerical aperture
% P: DAQ information
    % P.startAngle: first angle of CC plane wave emision 
    % P.dAngle: angle increment of plane wave emision
    % P.numAngles: number of angles for Coherence compounding
    % P.CCangle: angle range for coherence compounding, in degree
    % P.startDepth: acquistion start depth, in wavelength
    % P.actZsamples: number of samples for each pulse-echo acquisition
    % P.nCh: number of probe channels
    % P.pitch: pitch between transducer elements, mm
    % P.frequency: center frequency, MHz
    % P.vSound: speed of sound, m/s
% (xCoor, zCoor): the coordinates of the image to be formed
%% output
% IndCtriMatrix: sampling points index matrix of all for beamformed image
% ApodChn: apodized channel/element set by numerical aperture
function [IndCtriMatrix,ApodChn]=RefIndMatrix(P,xCoor,zCoor, NA)
tic;
if nargin<4
    NA=1;
end
%% iP: transducer parameters
iP.vSound=P.vSound; % speed of sound
iP.nCh=P.nCh; % number of probe channels
iP.pitch=P.pitch; % pitch between transducer elements, mm
iP.Wavelength=P.vSound*1e3/(P.frequency*1e6); % mm
iP.dSample=iP.Wavelength/(P.nSmplPerWvlnth*P.nRFref); % distance between samplling points, mm, the default is 4 sampling points per wavelength
iP.nZsample= P.actZsamples*P.nRFref; % number of samples for each pulse-echo acquisition
iP.delay=P.startDepth*iP.Wavelength; % delay for the first acquire point
%% beamforming image coordinates
nx=length(xCoor);
nz=length(zCoor);
%% beamforming index matrix
IndCtriMatrix=zeros(nz,nx,P.nCh,P.numAngles);
ApodChn=zeros(nz,nx,P.nCh,P.numAngles);
for iAgl=1:P.numAngles
    iP.angle=P.startAngle+(iAgl-1)*P.dAngle;
    for ix=1:nx            % main loop, coordinates of the x points
        for iz=1:nz
            [IndCtriMatrix(iz,ix,:,iAgl),ApodChn(iz,ix,:,iAgl)]=CtribtPixel(iP,NA, xCoor(ix),zCoor(iz));      % calculate delays and apodization
        end%
    end
end
toc