function plotMIPs(data, gamcp) % expects 4D input (x, y, z, frames)
    % gamcp = gamma compression power
    
    figure; imagesc(squeeze(max(data, [], 1))' .^ gamcp); colormap hot; colorbar
    title('xz MIP')
    xlabel('y pixels')
    ylabel('z pixels')

    figure; imagesc(squeeze(max(data, [], 2))' .^ gamcp); colormap hot; colorbar
    title('yz MIP')
    xlabel('x pixels')
    ylabel('z pixels')

    figure; imagesc(squeeze(max(data, [], 3))' .^ gamcp); colormap hot; colorbar
    title('xy MIP')
    xlabel('x pixels')
    ylabel('y pixels')

end