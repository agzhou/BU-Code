%% a single plane wave beamforming (iP)
% RF: acquired RF data for a single plane wave detection
    % size(RF)=[iP.nZsample,iP.nCh]
% iP: emission and receive information for the ith plane wave
    % iP: transducer parameters
    % iP.angle=P.startAngle+(i-1)*P.dAngle: plane wave emission angle
    % iP.vSound=P.vSound: speed of sound
    % iP.nCh=P.nCh: number of probe channels
    % iP.pitch=P.pitch: pitch between transducer elements, mm
    % iP.Wavelength=P.vSound*1e3/(P.frequency*1e6), mm
    % iP.dSample=iP.Wavelength/4, distance between samplling points, mm, the default is 4 sampling points per
    % iP.nZsample= P.maxZsamples: number of samples for each pulse-echo acquisition
    % wavelength
    % iP.delay=P.startDepth*iP.Wavelength, delay for the first acquie point
% NA=rAperture/zDistance: reconstruction numerical aperture
% P: DAQ information
    % P.startAngle: first angle of CC plane wave emision 
    % P.dAngle: angle increment of plane wave emision
    % P.numAngles: number of angles for Coherence compounding
    % P.CCangle: angle range for coherence compounding, in degree
    % P.startDepth: acquistion start depth, in wavelength
    % P.maxZsamples: number of samples for each pulse-echo acquisition
    % P.nCh: number of probe channels
    % P.pitch: pitch between transducer elements, mm
    % P.frequency: center frequency, MHz
    % P.vSound: speed of sound, m/s
% (xCoor, zCoor): the coordinates of the image to be formed
function iBF=PWBeamform(RF,iP,xCoor,zCoor, NA)
if nargin<5
    NA=1;
end
nx=length(xCoor);
nz=length(zCoor);
iBF=zeros(nz,nx);      % init beamforming matrix                  
for ix=1:nx            % main loop, coordinates of the x points
    for iz=1:nz
        [indLTravel,Apod]=CtribtPixel(iP,NA, xCoor(ix),zCoor(iz));      % calculate delays and apodization             %
        for iCh=1:iP.nCh
            iBF(iz,ix)=iBF(iz,ix)+RF(indLTravel(iCh),iCh)*Apod(iCh);   % delay and sum the RF data 
        end
    end
end

iBF=hilbert(iBF);    %Hilbert transform, analitical signal