function Fuse2Images(img_bk,img_olap, cRangeBk,cRangeOlap, xCoor, yCoor, shhold)
[VzNegCmap, VzPotCmap]=Colormaps_fUS;
%% 1) create figure and store handle
hF = figure;
%% 2) create axes for pcolor and store handle
hAxes1 = axes;
%set colormap for pcolor axes
% colormap(hAxesP,Vzcmap);
colormap(hAxes1,VzPotCmap);
% plot background figure
imagesc(xCoor,yCoor,img_bk); 
%create color bar and set range for color
caxis(hAxes1,cRangeBk);
axis on
axis equal
axis tight
xlabel('x [mm]')
ylabel('z [mm]')
title('Fused Image')
%% 3) create axes for the overlap axes
hAxes2 = axes;
%set visibility for axes to 'off' so it appears transparent
axis(hAxes2,'off')
%set colormap for overlap axes
colormap(VzNegCmap);
% plot overlap figure
h1=imagesc(xCoor,yCoor,img_olap);
alpha(h1,0.6*double((img_olap)>shhold)); % -mean(mean(img_olap)), transparency settings
caxis(hAxes2,cRangeOlap);
axis equal
axis tight
axis off
%link the two overlaying axes so they match at all times to remain accurate
linkaxes([hAxes1,hAxes2]);


% %% 4) overlap with the third image
% hAxes3 = axes;
% %set visibility for axes to 'off' so it appears transparent
% axis(hAxes3,'off')
% %set colormap for overlap axes
% % colormap(hAxesCM,Dcmap);
% colormap(hAxes3,cool);
% % plot overlap figure
% h1=imagesc(img_olap2);
% alpha(h1,double(abs(img_olap2)>Dshhold2)); % -mean(mean(img_olap))
% 
% caxis(hAxes3,cRangeOlap2);
% axis off
% axis equal
% axis tight
% %link the two overlaying axes so they match at all times to remain accurate
% linkaxes([hAxes1,hAxes2, hAxes3]);



%% define colormap %%%
% cm1a = [32:-1:1]'*[0 1 0]/32; cm1a(:,3) = 1; cm1b = [32:-1:1]'*[0 0 1]/32;
% cm2 = hot(64);
% % colormap([cm1a; cm1b; cm2])
% %%%%% gray scale background %%%%
% img_bk_log=log(img_bk)-min(min(log(img_bk)));
% img_bk_gray=ind2rgb(round(img_bk_log*100),gray(256));
% figure,
% imagesc(img_bk_gray); hold on; h1=imagesc(img_olap);
% alpha(h1,double(abs(img_olap-mean(mean(img_olap)))>40))
% colormap([cm1a; cm1b; cm2])
% caxis(cRange1)
% axis equal;
% axis off
% 
% %%%%% color scale overlaped image%%%%
% img_bk_col=img_bk;
% figure,
% imagesc(img_bk_col); hold on; h1=imagesc(img_olap);
% alpha(h1,double(abs(img_olap-mean(mean(img_olap)))>30))
% 
% colormap([cm1a; cm1b; cm2])
% caxis(cRange1)
% axis equal;
% axis off
