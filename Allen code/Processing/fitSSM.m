%% Description: 
%       Find the optimal cutoffs for the Spatial Similarity Matrix (SSM).
%       See the supplemental materials for Baranger et al., 2018 
%       (https://doi.org/10.1109/tmi.2018.2789499)
%
%       This method calculates the normalized correlation between the SSM
%       and all the possibilities of juxtaposed squares as defined by the a
%       and b parameters. The maximum of the normalized correlation map
%       corresponds to the optimal a and b parameters.

% Inputs:
%       SSM (symmetric matrix, can use plotSSM.m to get it)
%       (optional) showSSM: true or false, to plot the SSM or not
%       (optional) globalBounds: [smallest a, largest a, largest b] to use for the parameter search
% Outputs:
%       XN: normalized correlation map
%       a_opt: optimal a parameter
%       b_opt: optimal b parameter


function [XN, a_opt, b_opt] = fitSSM(SSM, varargin)

    np = size(SSM, 1); % get the size/"# points" of the SSM, which is expected to be square

    % Initialize the correlation map between the SSM and the family of juxtaposed squares
    X = zeros(np);  % Unnormalized correlation map
    XN = zeros(np); % Normalized correlation map
    SSM_term = SSM - mean(SSM, 'all'); % Pre-calculate the term for speed
    X_factor = 1/np^2; % Normalization factor for X

    globalBounds = [1, np, np]; % Default global bounds

    if nargin > 1
        showSSM = varargin{1};
        if nargin > 2
            globalBounds = varargin{2};
        end
    end

    % Iteratively calculate the unnormalized and normalized correlation
    % maps between the SSM and the families of squares
    for a = globalBounds(1):globalBounds(2)
%         disp(a)
%         for b = 1:np
        for b = a:globalBounds(3) % by definition, b > a
            alpha = createAlpha(np, a, b);
            X(a, b) = X_factor .* sum(SSM_term .* (alpha - mean(alpha, 'all')), 'all');
            XN(a, b) = X(a, b) ./ sqrt( X_factor .* sum(SSM_term .* SSM_term, 'all') ) ...
                ./ sqrt( X_factor .* sum((alpha - mean(alpha, 'all')) .* (alpha - mean(alpha, 'all')), 'all') );
        end
    end

    % Find the max of the normalized correlation map and get the indices
    [~, XN_max_index] = max(XN, [], 'all');

    % Turn the indices into the optimal a and b parameters
    [a_opt, b_opt] = ind2sub(size(XN), XN_max_index);

    % Plot the SSM with the squares
    if showSSM
        plotSSMDivisions(SSM, a_opt, b_opt)
    end
end



%% Helper functions
% Define the "family of juxtaposed squares alpha"
function alpha = createAlpha(np, a, b)
    alpha = zeros(np);
    alpha(1:a, 1:a) = 1;
    alpha(a:b, a:b) = 1;
end

% Optional stuff to plot the optimized divisions
function plotSSMDivisions(SSM, a_opt, b_opt)
    np = size(SSM, 1); % get the size/"# points" of the SSM, which is expected to be square
    figure; imagesc(SSM); title('Spatial Similarity Matrix'); colorbar; clim([0, 1]); axis square
    
    line_color = 'black';

    % First square
    line(a_opt .* ones(np, 1), 1:np, 'Color', line_color)
    line(1:np, a_opt .* ones(np, 1), 'Color', line_color)
    % Second square
    line(b_opt .* ones(np, 1), 1:np, 'Color', line_color)
    line(1:np, b_opt .* ones(np, 1), 'Color', line_color)
end