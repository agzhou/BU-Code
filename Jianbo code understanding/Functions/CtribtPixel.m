%% obtained contribution pixel index of all channels for position (x, z)
% P.nCh: number of probe channels
% P.pitch: pitch between transducer elements, mm
% P.frequency: center frequency, MHz
% P.vSound: speed of sound, m/s
% P.startAngle: first angle of CC plane wave emision 
% P.dAngle: angle increment of plane wave emision 
% P.numAngles: number of angles for Coherence compounding
% P.CCangle: angle range for coherence compounding, in degree
% P.startDepth: acquistion start depth, in wavelength
% P.maxZsamples: number of samples for each pulse-echo acquisition

function [indLtravel,Apod]=CtribtPixel(iP,NA,x,z)
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
% (x, z): current beamforming pixel
       
xElem = (0.5:iP.nCh)*iP.pitch;                                              % transducer element transverse (x) coordinates
if iP.angle<0
    LTravel= (cos(iP.angle)*z+(xElem(end)-x)*sin(abs(iP.angle))+sqrt(z^2+(xElem-x).*(xElem-x)));  % distance between the pixel(x, z) and each transducer element (xCoor, 0), mm
else
    LTravel= (cos(iP.angle)*z+x*sin(iP.angle)+sqrt(z^2+(xElem-x).*(xElem-x)));  % distance between the pixel(x, z) and each transducer element (xCoor, 0), mm
end
indLtravel= round((LTravel-iP.delay)/iP.dSample);                           % index of the contribution sampling point for every element
Apod= abs((x-xElem)/z) < NA ;                                               % apodisation aperture
indLtravel(find(indLtravel<=0 | indLtravel>=iP.nZsample))=1;                % limits (contribution sampling point index=1 if data does not exist)
