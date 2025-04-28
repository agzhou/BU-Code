function [cmap] = colormap_ULM

    numSteps = 10; % How many steps total for each color increment, then times 3 for rgb
    cmap = []; % Initialize the colormap

    % Black to blue to green to red. Black corresponds to small values.
    for c = 3:-1:1 % 1 - R, 2 - G, 3 - B
        temp = zeros(numSteps, 3);
        temp(:, c) = linspace(0, 1, numSteps);
        cmap = cat(1, cmap, temp);
    end
    




end