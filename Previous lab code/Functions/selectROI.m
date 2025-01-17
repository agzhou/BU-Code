function BW = selectROI(coefmap_GG0)

newSlect = 1;
if newSlect==1
    figure; imagesc(coefmap_GG0); axis image;colormap(jet);%caxis([0,1]);
    [loc_x,loc_y]=ginput(6); % [x z]
    BW=roipoly(coefmap_GG0,loc_x,loc_y);
end
Fig = figure; set(Fig,'Position',[600 600 500 350])
imagesc(BW.*coefmap_GG0); axis image;%caxis([-1,1]); 
colorbar;colormap('jet'); title('ROI');

end