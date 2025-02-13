% Interpolate between two points [z, x] which represent the start and end
% of a track between paired bubbles in consecutive frames
% Need to round it in case the output is not a natural number
function [zcInterp, xcInterp] = ULM_interp2D(sourceCoords, targetCoords)
    
    zcInterpStart = sourceCoords(1);
    zcInterpEnd = targetCoords(1);
    % Note: the order (source or target) between the start and end has to be consistent between the z and x interpolations
    xcInterpStart = sourceCoords(2); 
    xcInterpEnd = targetCoords(2);

    numInterpPoints = max([abs(zcInterpEnd - zcInterpStart) + 1, abs(xcInterpEnd - xcInterpStart) + 1]);
    zcInterp = linspace(zcInterpStart, zcInterpEnd, numInterpPoints)'; % do this so you get the right number of points
    zcInterp = round(zcInterp);
    xcInterp = linspace(xcInterpStart, xcInterpEnd, numInterpPoints)'; % do this so you get the right number of points
    xcInterp = round(xcInterp);
end