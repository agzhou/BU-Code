function Media = sortMedia(Media, dim)


    [~, ind_Media_sorted] = sort(Media.MP(:, dim), 1);
    Media = Media.MP(ind_Media_sorted, :);
    %%
%     temp_Media(:, dim) = temp_Media(:, dim) + dist_per_frame_wl;

end