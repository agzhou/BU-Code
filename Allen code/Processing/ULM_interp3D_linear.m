% Interpolate between two points [x, y, z] which represent the start and end
% of a track between paired bubbles in consecutive frames.
% Need to round it in case the output is not a natural number.

function [pts] = ULM_interp3D_linear(sourceCoords, targetCoords, speedTemp)
    pts = [];
    for nvp = 1:size(speedTemp, 1)
        xcInterpStart = sourceCoords(nvp, 1);
        xcInterpEnd = targetCoords(nvp, 1);
        % Note: the order (source or target) between the start and end has to be consistent between the z and x interpolations
        ycInterpStart = sourceCoords(nvp, 2); 
        ycInterpEnd = targetCoords(nvp, 2);
        zcInterpStart = sourceCoords(nvp, 3); 
        zcInterpEnd = targetCoords(nvp, 3);
    
        numInterpPoints = max([abs(xcInterpEnd - xcInterpStart) + 1, abs(ycInterpEnd - ycInterpStart) + 1, abs(zcInterpEnd - zcInterpStart) + 1]);
%         if numInterpPoints > 20
%             warning(['bad', num2str(ti)])
%             disp(zcInterpEnd - zcInterpStart)
% %             disp(xcInterpEnd - xcInterpStart)
%         end
        xcInterp = linspace(xcInterpStart, xcInterpEnd, numInterpPoints)'; % do this so you get the right number of points
        xcInterp = round(xcInterp);

        ycInterp = linspace(ycInterpStart, ycInterpEnd, numInterpPoints)'; % do this so you get the right number of points
        ycInterp = round(ycInterp);

        zcInterp = linspace(zcInterpStart, zcInterpEnd, numInterpPoints)'; % do this so you get the right number of points
        zcInterp = round(zcInterp);
        

        pts = [pts; [xcInterp, ycInterp, zcInterp, repmat(speedTemp(nvp), numInterpPoints, 1)]];
    end
end