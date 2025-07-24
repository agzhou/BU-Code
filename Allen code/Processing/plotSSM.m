%% Description:
%       Calculate and plot the spatial similarity matrix as defined in
%       Baranger et al., 2018 (https://doi.org/10.1109/tmi.2018.2789499)

function [SSM] = plotSSM(U)
    SSM = corrcoef(abs(U));
    figure; imagesc(SSM); title('Spatial Similarity Matrix'); colorbar; clim([0, 1]); axis square


end