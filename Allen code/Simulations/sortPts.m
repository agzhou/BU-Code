function pts = sortPts(pts, dim)

    [~, ind_pts_sorted] = sort(pts(:, dim), 1);
    pts = pts(ind_pts_sorted, :);

end