function plotPoints(pts, SP)
    figure;
    scatter3(pts(:, 1), pts(:, 2), pts(:, 3), '.')
%     axis square
%     axis equal

    % xlim([min(pts(:, 1)), max(pts(:, 1))])
    % ylim([min(pts(:, 2)), max(pts(:, 2))])
    % zlim([min(pts(:, 3)), max(pts(:, 3))])

%     xlim([SP.xstart - SP.vesselDiam/2, SP.xstart + SP.vesselDiam/2])
%     ylim([SP.ystart - SP.vesselDiam/2, SP.ystart + SP.vesselDiam/2])
%     zlim([SP.zstart - SP.vesselLength/2, SP.zstart + SP.vesselLength/2])

    xlabel("x")
    ylabel("y")
    zlabel("z")

    % xlabel("x [m]")
    % ylabel("y [m]")
    % zlabel("z [m]")

%     fontsize(20, "points")
end