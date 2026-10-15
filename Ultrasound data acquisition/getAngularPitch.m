% Calculate the proper angular pitch to reject grating lobes outside of
% some lateral aperture

% Inputs in m, output in rad or deg

function [dAngle] = getAngularPitch(wl, lateralAperture, degOrRad)
    if nargin < 3
        degOrRad = 'rad';
    end

    dAngle = wl / lateralAperture;
    if isequal(degOrRad, 'deg')
        dAngle = rad2deg(dAngle);
    end
end