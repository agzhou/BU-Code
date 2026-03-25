% Turn the single-dimension angle list into a full list of x and y angles
function angles = listToAngles(anglesList)
    [anglesX, anglesY] = meshgrid(anglesList, anglesList);
    angles = [anglesX(:), anglesY(:)];
end