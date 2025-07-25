%% Description:
%       Calculate and plot the spatial similarity matrix as defined in
%       Baranger et al., 2018 (https://doi.org/10.1109/tmi.2018.2789499)

% Inputs:
%       U: spatial singular vectors in a matrix
%       (optional) showSSM: true or false, to plot the SSM or not

function [SSM] = plotSSM(U, varargin)
    if nargin > 1
        showSSM = varargin{1};
    end
    SSM = corrcoef(abs(U));
    if showSSM
        figure; imagesc(SSM); title('Spatial Similarity Matrix'); colorbar; clim([0, 1]); axis square
    end

end