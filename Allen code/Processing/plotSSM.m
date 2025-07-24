function [] = plotSSM(U)
    SSM = corrcoef(abs(U));
    figure; imagesc(SSM); title('Spatial Similarity Matrix'); colorbar; clim([0, 1]); axis square


end