% Go from bVelocityM to bSum - bubble density map
function [bSum] = densityMap3D(bVelocityM, img_size, startFrame)
    bSum = zeros(img_size(1), img_size(2), img_size(3)); % Initialize the bubble density map variable
    
    for n = startFrame:length(bVelocityM) % Go through each track collection and plot
        tempBuf = bVelocityM{n};
        for tn = 1:size(tempBuf, 1) % Go through each track number tn
            trackTemp = squeeze(tempBuf(tn, :, :))';
    
            % Add 1 to a pixel's count if the current track intersects it
            bSum(trackTemp(1, 1), trackTemp(1, 2), trackTemp(1, 3)) = bSum(trackTemp(1, 1), trackTemp(1, 2), trackTemp(1, 3)) + 1;
            for iti = 1:size(trackTemp, 1) % inside track index
                bSum(trackTemp(iti, 4), trackTemp(iti, 5), trackTemp(iti, 6)) = bSum(trackTemp(iti, 4), trackTemp(iti, 5), trackTemp(iti, 6)) + 1;
            end
        end
    end
end