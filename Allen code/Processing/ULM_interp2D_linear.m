% Interpolate between two points [z, x] which represent the start and end
% of a track between paired bubbles in consecutive frames
% Need to round it in case the output is not a natural number
% function [pts] = ULM_interp2D_linear(sourceCoords, targetCoords, vTemp)
function [pts] = ULM_interp2D_linear(sourceCoords, targetCoords, vTemp, ti)
pts = [];
    for nvp = 1:length(vTemp)
        zcInterpStart = sourceCoords(nvp, 1);
        zcInterpEnd = targetCoords(nvp, 1);
        % Note: the order (source or target) between the start and end has to be consistent between the z and x interpolations
        xcInterpStart = sourceCoords(nvp, 2); 
        xcInterpEnd = targetCoords(nvp, 2);
    
        numInterpPoints = max([abs(zcInterpEnd - zcInterpStart) + 1, abs(xcInterpEnd - xcInterpStart) + 1]);
        if numInterpPoints > 20
            warning(['bad', num2str(ti)])
            disp(zcInterpEnd - zcInterpStart)
%             disp(xcInterpEnd - xcInterpStart)
        end
        zcInterp = linspace(zcInterpStart, zcInterpEnd, numInterpPoints)'; % do this so you get the right number of points
        zcInterp = round(zcInterp);
        xcInterp = linspace(xcInterpStart, xcInterpEnd, numInterpPoints)'; % do this so you get the right number of points
        xcInterp = round(xcInterp);

        pts = [pts; [zcInterp, xcInterp, repmat(vTemp(nvp), length(zcInterp), 1)]];
    end
end