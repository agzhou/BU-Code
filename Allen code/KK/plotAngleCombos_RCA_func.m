% Plot histograms of k = k_out - k_in for different receive angle schemes
% (in 2D case: k_x = k_x,out - k_x,in)
% And here, plotting delta_theta = theta_out - theta_in is equivalent for
% small angles.

function [varargout] = plotAngleCombos_RCA_func(anglesTX, anglesRX)
    maxAngle = max([anglesTX(:); anglesRX(:)]);
    binlimits = rad2deg([-maxAngle*2, maxAngle*2]); % Test histogram bin limits
    % binwidth = maxAngle * 2 / size(anglesTX, 1) / size(anglesRX, 1);
    binwidth = rad2deg(1 * pi/180);

    delta_theta = rad2deg(calcDeltaThetaRCA(anglesTX, anglesRX));

    if nargout > 0
        varargout{1} = delta_theta; % Return delta_theta as an optional output
    end
    
    % k_x = f*2*pi .* sind(delta_theta) ./ c;
    % figure; imagesc(delta_theta); colorbar; axis image; xlabel('TX angle'); ylabel('RX angle'); title('Delta theta [deg]')
    % figure; imagesc(k_x); colorbar; axis image; xlabel('TX angle'); ylabel('RX angle'); title('k_x [radians/m]')
    
    % Plot histograms
    % figure; histogram2(delta_theta(:, 1), delta_theta(:, 2), BinMethod="integers"); title('Delta theta counts'); xlabel('Delta theta [deg]'); ylabel('Counts')
    figure; 
    h = histogram2(delta_theta(:, 1), delta_theta(:, 2), ...
        'XBinLimits', binlimits, ...
        'YBinLimits', binlimits, ...
        'BinWidth', binwidth, ...
        'DisplayStyle', 'tile', ...
        'ShowEmptyBins', 'on');
    
    h.EdgeColor = [0.88, 0.88, 0.88];   % light gray grid
    
    title('Delta theta counts'); 
    xlabel('Delta theta (x) [deg]'); 
    xlim([-25, 25])
    ylabel('Delta theta (y) [deg]')
    ylim([-25, 25])
    zlabel('Counts')
    axis image
    colorbar
    
    % Make zero white
    c = parula(256);
    c = [1 1 1; c];   % prepend white
    colormap(c)
    
    view(90, 90)
    
    % figure; histogram(k_x); title('Delta k_x'); xlabel('Delta k_x [radians/m]'); ylabel('Counts')
    % % Plot histograms
    % figure; histogram(delta_theta, BinMethod="integers"); title("Delta theta counts, with j = " + num2str(j)); xlabel('Delta theta [deg]'); ylabel('Counts')
    % % figure; histogram(k_x); title("Delta k_x counts, with j = " + num2str(j)); xlabel('Delta k_x [radians/m]'); ylabel('Counts')
end

%% Helper functions
function delta_theta = calcDeltaThetaRCA(anglesTX, anglesRX)
    naTX = size(anglesTX, 1)/2; % # of transmit angles in one dimension
    naRX = size(anglesRX, 1)/2; % # of receive angles in one dimension
    delta_theta = [];
    % CR
    for tai = 1:naTX
        delta_theta = [delta_theta; anglesRX(1:naRX, :) - anglesTX(tai, :)];
    end
    % RC
    for tai = naTX + 1:naTX*2
        delta_theta = [delta_theta; anglesRX(naRX + 1:naRX*2, :) - anglesTX(tai, :)];
    end
end