function Fuse3Images(img_bk,img_olap1, img_olap2, cRangeBk,cRangeOlap1,cRangeOlap2, shhold, isInver)
%% %%%%%%%%%%%%% example %%%%%%%%%%%%%%%%%%%%%%%%
% xCoor=1:512;
% zCoor=1:320;
% figure,Fuse2Images(Vz(:,:,1),Vz(:,:,1),[-30 30],[-30 30],xCoor,zCoor,2.5)
%% 1) Colormap %%%%%%%%%%%%%%%%%%%%%%%%%%%%
[VzCmap, VzCmapDn,VzCmapUp]=Colormaps_fUS;
%% 2) create axes for pcolor and store handle
hAxes1 = axes;
% plot background figure
imagesc(img_bk); 
%create color bar and set range for color
caxis(hAxes1,cRangeBk);
%set colormap for pcolor axes
%colormap(hAxesP,Vzcmap);
colormap(hAxes1,gray);
axis on
axis equal
axis tight
axis off
colorbar('Ticks',[-10 10], 'TickLabels',{[],[]} );
% title('Fused Image')
%% 3) create axes for the overlap axes
hold on;
hAxes2 = axes;
% plot overlap figure
h1=imagesc(img_olap1);
alpha(h1,0.7*double((img_olap1)>shhold)); % -mean(mean(img_olap)), transparency settings
caxis(hAxes2,cRangeOlap1);
%set colormap for overlap axes
if isInver == 0
colormap(hAxes2,VzCmapDn);
else
colormap(hAxes2,flipud(VzCmapUp));
end
axis equal
axis tight
colorbar('Ticks',0:0.2:1, 'TickLabels',{[],[]});
%set visibility for axes to 'off' so it appears transparent
axis(hAxes2,'off')
axis off

%% 4) create axes for the overlap axes
hold on;
hAxes3 = axes;
% plot overlap figure
h1=imagesc(img_olap2);
alpha(h1,0.7*double((img_olap2)>shhold)); % -mean(mean(img_olap)), transparency settings
caxis(hAxes3,cRangeOlap2);
%set colormap for overlap axes
if isInver == 0
colormap(hAxes3,flipud(VzCmapUp));
else
colormap(hAxes3,VzCmapDn);
end
axis equal
axis tight
colorbar('Ticks',0:0.2:1, 'TickLabels',{[],[]});
%set visibility for axes to 'off' so it appears transparent
axis(hAxes3,'off')
axis off
%link the three overlaying axes so they match at all times to remain accurate
linkaxes([hAxes1,hAxes2, hAxes3]);