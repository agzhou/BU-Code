COLOR=[1 0 0;
    0.08 0.17 0.55
    0.31 0.31 0.31];
addpath D:\g1_based_fUS\A-US-fUS-g1\SubFunctions
% load actmap_GG0, actmap_GGV0, actmap_PDI0, actmap_GGFV0
% load coefmap_GG0, 
% load dir

%% 1. masked Correlation map: CBFv based fUS
viscrop = [20, 215];
img_bk = log(mean(mtrialPDI(:,:,:),3)./eqNoise); 
img_bk = img_bk(:,viscrop(1):viscrop(2));
BB=ones(3,3);
BB(2,2)=2;
BB=BB/10;
actmap_covn = convn(actmap_PDI0(:,viscrop(1):viscrop(2)),BB,'same');
% substitute actmap_GG0, actmap_GGV0, actmap_PDI0, actmap_GGFV0
img_olap = actmap_covn;
% img_olap = actmap_GG0(:,viscrop(1):viscrop(2));

cRangeBk = [min(img_bk(:)), max(img_bk(:))*0.93];
cRangeOlap = [0.1, 0.9];
shhold = 0.15;
fig = figure;
set(fig, 'Position', [300 300 500 400])
Fuse2Images(img_bk,img_olap, cRangeBk,cRangeOlap, shhold)
colorbar('eastoutside', 'Ticks', cRangeOlap);%,'TickLabels',{'0.3','0.9'});

%% 2. roi based time course
ratiotrialroi = ROIGGV.ratio;
mratiotrialroi = ROIGGV.m;
semratiotrialroi = ROIGGV.sem;
minc = 90;
maxc = 120;

% tCoor = [1:trial.nlength]-trial.nRest-1;
% figure; plot(ratiotrialroi);
% outtrial = find(abs(mean(ratiotrialroi-mratiotrialroi,1))>mean(mean(ratiotrialroi-mratiotrialroi,1),2)+7)
% outtrial = find(abs(mean(ratiotrialroi-mratiotrialroi,1))== max(abs(mean(ratiotrialroi-mratiotrialroi,1))))
% 
% ratiotrialroi1 = ratiotrialroi(:,[1:outtrial-1,outtrial+1:end]);
% mratiotrialroi1 = mean(ratiotrialroi1,2);
% semratiotrialroi1 = std(ratiotrialroi1,1,2)./sqrt(trial.n-numel(outtrial));

ratiotrialroi1 = ratiotrialroi;
mratiotrialroi1 = mratiotrialroi;
semratiotrialroi1 = semratiotrialroi;

fig = figure;
set(fig, 'Position', [700 700 400 300])
p1 = plot(tCoor, ratiotrialroi1, 'Color', [0 0 0 0.3]); 
hold on; plot(tCoor, mratiotrialroi1,'k','LineWidth', 1.3);
set(gca, 'ylim', [minc, maxc]);
xlabel('Time [s]')
ylabel('Relative change %')
set(gca, 'FontSize', 10);
hold on;
yyaxis right; ax = gca; plot(tCoor, smooth(stimhrf,4)*median(mratiotrialroi)/10+100,'Color','r'); 
ax.YColor = 'w'; set(gca, 'ylim',[minc,maxc], 'YTick', [], 'YTickLabel',[]);
hold on;
plot(0:trial.nStim,ones(trial.nStim+1,1)*95,'-k','LineWidth',2);

fig = figure;
set(fig, 'Position', [700 700 400 300])
Color.Shade='k';
Color.ShadeAlpha=0.25;
Color.Line='k';
ShadedErrorbar(mratiotrialroi1',semratiotrialroi1',tCoor,Color); 
set(gca, 'ylim', [minc, maxc]);
xlabel('Time [s]')
ylabel('Relative change %')
set(gca, 'FontSize', 10);
hold on;
yyaxis right; ax = gca; plot(tCoor, smooth(stimhrf,4)*median(mratiotrialroi)/10+100,'Color','r'); 
ax.YColor = 'w'; set(gca, 'ylim',[minc,maxc], 'YTick', [], 'YTickLabel',[]);
hold on;
plot(0:trial.nStim,ones(trial.nStim+1,1)*95,'-k','LineWidth',2);

% %% roi contour
% fig = figure;
% set(fig, 'Position', [300 300 500 400]) 
% h1 = imagesc(sign(abs(diff(BWGGFV0(:,viscrop(1):viscrop(2)))))); colorbar;colormap('gray');axis image;
% alpha(h1, 1*(abs(diff(BWGGFV0(:,viscrop(1):viscrop(2))))>0)); axis off;


%% 3.1.statiscal results: CBFv based fUS 
veinAmap_GG0 = veinAmap_GG0.*MSKGG;
arteriesAmap_GG0 = arteriesAmap_GG0.*MSKGG;

mvPC = mean(veinAmap_GG0(veinAmap_GG0~=0));
maPC = mean(arteriesAmap_GG0(arteriesAmap_GG0~=0));

fig = figure;
set(fig,'Position',[600 600 400 300])
h1 = histogram(veinAmap_GG0(veinAmap_GG0~=0),'BinWidth',0.03);
hold on;
h2 = histogram(arteriesAmap_GG0(arteriesAmap_GG0~=0),'BinWidth',0.03);
hold on;
plot([mvPC,mvPC], [0,40],'-b', 'LineWidth',1.5, 'Color',[0 0.5 1]);
hold on;
plot([maPC,maPC], [0,40],'-r','LineWidth',1.5, 'Color',[1 0.2 0.5]);
xlabel('Correlation coefficient');
ylabel('Occurance (pixels)' )

h1.FaceColor = [0.08 0.17 0.55];%[0 0.5 1];%[0.3010 0.7450 0.9330];%[0, 0, 1];
%h1.EdgeColor = 'w';
h1.FaceAlpha = 0.3;

h2.FaceColor = [1,0,0];%[1 0.2 0.5];%[0.8500 0.3250 0.0980];%[0, 1, 0];
%h2.EdgeColor = 'w';
h2.FaceAlpha = 0.3;

% h1.Normalization = 'probability';
% h1.BinWidth = 0.04;
% h2.Normalization = 'probability';
% h2.BinWidth = 0.04;
%% 3.2.statiscal results: CBV based fUS 
veinAmap_GGV0 = veinAmap_GGV0.*MSKGGV;
arteriesAmap_GGV0 = arteriesAmap_GGV0.*MSKGGV;

mvPC = mean(veinAmap_GGV0(veinAmap_GGV0~=0));
maPC = mean(arteriesAmap_GGV0(arteriesAmap_GGV0~=0));

fig = figure;
set(fig,'Position',[600 600 400 300])
h1 = histogram(veinAmap_GGV0(veinAmap_GGV0~=0),'BinWidth',0.03);
hold on;
h2 = histogram(arteriesAmap_GGV0(arteriesAmap_GGV0~=0),'BinWidth',0.03);
hold on;
plot([mvPC,mvPC], [0,60],'-b', 'LineWidth',1.5, 'Color',[0 0.5 1]);
hold on;
plot([maPC,maPC], [0,60],'-r','LineWidth',1.5, 'Color',[1 0.2 0.5]);
xlabel('Correlation coefficient');
ylabel('Occurance (pixels)' )

h1.FaceColor = [0.08 0.17 0.55];%[0 0.5 1];%[0.3010 0.7450 0.9330];%[0, 1, 0];
%h1.EdgeColor = 'w';
h1.FaceAlpha = 0.3;

h2.FaceColor = [1 0 0];%[1 0.2 0.5];%[0.8500 0.3250 0.0980];%[0, 1, 0];
%h2.EdgeColor = 'w';
h2.FaceAlpha = 0.3;
%% 3.3.statiscal results: CBFv based fUS 
ampGG_v = ampGG_v.*MSKGG;
ampGG_a = ampGG_a.*MSKGG;

mvPC = mean(ampGG_v (ampGG_v ~=0));
maPC = mean(ampGG_a(ampGG_a~=0));

fig = figure;
set(fig,'Position',[600 600 400 300])
h1 = histogram(ampGG_v (ampGG_v ~=0),'BinWidth',0.03);
hold on;
h2 = histogram(ampGG_a(ampGG_a~=0),'BinWidth',0.03);
hold on;
plot([mvPC,mvPC], [0,40],'-b', 'LineWidth',1.5, 'Color',[0 0.5 1]);
hold on;
plot([maPC,maPC], [0,40],'-r','LineWidth',1.5, 'Color',[1 0.2 0.5]);
xlabel('Correlation coefficient');
ylabel('Occurance (pixels)' )

h1.FaceColor = [0.08 0.17 0.55];%[0 0.5 1];%[0.3010 0.7450 0.9330];%[0, 0, 1];
%h1.EdgeColor = 'w';
h1.FaceAlpha = 0.3;

h2.FaceColor = [1,0,0];%[1 0.2 0.5];%[0.8500 0.3250 0.0980];%[0, 1, 0];
%h2.EdgeColor = 'w';
h2.FaceAlpha = 0.3;

% h1.Normalization = 'probability';
% h1.BinWidth = 0.04;
% h2.Normalization = 'probability';
% h2.BinWidth = 0.04;
%% 3.4.statiscal results: CBFv based fUS 
ampGGV_v = ampGGV_v.*MSKGGV;
ampGGV_a = ampGGV_a.*MSKGGV;

mvPC = mean(ampGGV_v (ampGGV_v ~=0));
maPC = mean(ampGGV_a(ampGGV_a~=0));

fig = figure;
set(fig,'Position',[600 600 400 300])
h1 = histogram(ampGGV_v (ampGGV_v ~=0),'BinWidth',0.03);
hold on;
h2 = histogram(ampGGV_a(ampGGV_a~=0),'BinWidth',0.03);
hold on;
plot([mvPC,mvPC], [0,40],'-b', 'LineWidth',1.5, 'Color',[0 0.5 1]);
hold on;
plot([maPC,maPC], [0,40],'-r','LineWidth',1.5, 'Color',[1 0.2 0.5]);
xlabel('Correlation coefficient');
ylabel('Occurance (pixels)' )

h1.FaceColor = [0.08 0.17 0.55];%[0 0.5 1];%[0.3010 0.7450 0.9330];%[0, 0, 1];
%h1.EdgeColor = 'w';
h1.FaceAlpha = 0.3;

h2.FaceColor = [1,0,0];%[1 0.2 0.5];%[0.8500 0.3250 0.0980];%[0, 1, 0];
%h2.EdgeColor = 'w';
h2.FaceAlpha = 0.3;

% h1.Normalization = 'probability';
% h1.BinWidth = 0.04;
% h2.Normalization = 'probability';
% h2.BinWidth = 0.04;
%% 3. arteries/veins masked Correlation map: CBFv based fUS
viscrop = [20, 215];
img_bk = log(mean(mtrialPDI(:,:,:),3)./eqNoise); 
img_bk = img_bk(:,viscrop(1):viscrop(2));
BB=ones(3,3);
BB(2,2)=2;
BB=BB/10;
actmap_covn1 = convn(arteriesAmap_GGFV0(:,viscrop(1):viscrop(2)),BB,'same');
actmap_covn2 = convn(veinAmap_GGFV0(:,viscrop(1):viscrop(2)),BB,'same');
% substitute actmap_GG0, actmap_GGV0, actmap_PDI0, actmap_GGFV0

img_olap1 = actmap_covn1;
img_olap2 = actmap_covn2;
cRangeBk = [min(img_bk(:)), max(img_bk(:))*0.93];
cRangeOlap1 = [0.1, 0.9];
cRangeOlap2 = [0.1, 0.9];
shhold = 0.15;
fig = figure;
set(fig, 'Position', [300 300 500 400])
Fuse3Images(img_bk,img_olap1, img_olap2, cRangeBk,cRangeOlap1,cRangeOlap2, shhold)
colorbar('eastoutside', 'Ticks', cRangeOlap);%,'TickLabels',{'0.3','0.9'});

%% 4. arteries/veins masked Amplitude map: CBFv based fUS

viscrop = [20, 215];
amp_a = ampGGV_a(:,viscrop(1):viscrop(2));
amp_v = ampGGV_v(:,viscrop(1):viscrop(2));
amp_a = amp_a + (amp_a==0);
amp_v = amp_v + (amp_v==0);

img_bk = log(mean(mtrialPDI(:,:,:),3)./eqNoise); 
img_bk = img_bk(:,viscrop(1):viscrop(2));
BB=ones(3,3);
BB(2,2)=2;
BB=BB/10;
actmap_covn1 = convn(amp_a,BB,'same');
actmap_covn2 = convn(amp_v,BB,'same');
% substitute actmap_GG0, actmap_GGV0, actmap_PDI0, actmap_GGFV0

% actmap_covn1 = ampGG_a(:,viscrop(1):viscrop(2));
% actmap_covn2 = ampGG_v(:,viscrop(1):viscrop(2));
minTickLable = min(min(amp_a(amp_a~=0)),min(amp_v(amp_v~=0)));
maxTickLable = max(max(amp_a(amp_a~=0)),max(amp_v(amp_v~=0)));
% minTick = min(min(actmap_covn1(actmap_covn1~=0)),min(actmap_covn2(actmap_covn2~=0)));
% maxTick = max(max(actmap_covn1(actmap_covn1~=0)),max(actmap_covn2(actmap_covn2~=0)));
minTick = 1;
maxTick = 1.3;
img_olap1 = actmap_covn1;
img_olap2 = actmap_covn2;
cRangeBk = [min(img_bk(:)), max(img_bk(:))*0.93];
cRangeOlap1 = [minTick, maxTick];%[0.5, 1.8];
cRangeOlap2 = [minTick, maxTick];%[0.5, 1.8];
shhold = 1.025;
fig = figure;
set(fig, 'Position', [300 300 500 400])
Fuse3Images(img_bk,img_olap1, img_olap2, cRangeBk,cRangeOlap1,cRangeOlap2, shhold)
colorbar('eastoutside', 'Ticks', [minTick, maxTick]);%, 'TickLabels',{num2str(minTickLable),num2str(maxTickLable)});

%% 4.1 arteries/veins masked Amplitude map: CBFv based fUS: original data
viscrop = [20, 215];
amp_a = ampGGF_a(:,viscrop(1):viscrop(2));
amp_v = ampGGF_v(:,viscrop(1):viscrop(2));

img_bk = log(mean(mtrialPDI(:,:,:),3)./eqNoise); 
img_bk = img_bk(:,viscrop(1):viscrop(2));

minc = min(min(amp_a(amp_a~=0)),min(amp_v(amp_v~=0)));
maxc = max(max(amp_a(amp_a~=0)),max(amp_v(amp_v~=0)));
minTick = minc;
maxTick = min(1.8, maxc);
minTickLable = minTick;
maxTickLable = maxTick;
% minTick = 0.5;
% maxTick = 1.5;
img_olap1 = amp_a;
img_olap2 = amp_v;
cRangeBk = [min(img_bk(:)), max(img_bk(:))*0.93];
cRangeOlap1 = [minTick, maxTick];%[0.5, 1.8];
cRangeOlap2 = [minTick, maxTick];%[0.5, 1.8];
shhold = 0.5;
fig = figure;
set(fig, 'Position', [300 300 500 400])
Fuse3Images(img_bk,img_olap1, img_olap2, cRangeBk,cRangeOlap1,cRangeOlap2, shhold)
colorbar('eastoutside', 'Ticks', [minTick, maxTick]);%, 'TickLabels',{num2str(minTickLable),num2str(maxTickLable)});


%% 5. arteries/veins roi based time course
ratiotrialroi = cat(3, ROIGGV.artery.ratio, ROIGGV.vein.ratio);
mratiotrialroi = [ROIGGV.artery.m, ROIGGV.vein.m];
semratiotrialroi = [ROIGGV.artery.sem, ROIGGV.vein.sem];
minc = 90;
maxc = 120;

tCoor = [1:trial.nlength]-trial.nRest-1;
% figure; plot(ratiotrialroi);
% outtrial = find(abs(mean(ratiotrialroi-mratiotrialroi,1))>mean(mean(ratiotrialroi-mratiotrialroi,1),2)+7)
% outtrial = find(abs(mean(ratiotrialroi-mratiotrialroi,1))== max(abs(mean(ratiotrialroi-mratiotrialroi,1))))
% ratiotrialroi1 = ratiotrialroi(:,[1:outtrial-1,outtrial+1:end],:);
% ratiotrialroi1 = ratiotrialroi;
% mratiotrialroi1 = mean(ratiotrialroi1,2);
% semratiotrialroi1 = std(ratiotrialroi1,1,2)./sqrt(trial.n-numel(outtrial));

ratiotrialroi1 = ratiotrialroi;
mratiotrialroi1 = mratiotrialroi;
semratiotrialroi1 = semratiotrialroi;

fig = figure;
set(fig, 'Position', [700 700 400 300])
Color.Shade=COLOR(1,:);
Color.ShadeAlpha=0.25;
Color.Line=COLOR(1,:);
ShadedErrorbar(mratiotrialroi1(:,1)',semratiotrialroi1(:,1)',tCoor,Color);% artery
set(gca, 'ylim', [minc, maxc]);
hold on
Color.Shade=COLOR(2,:);
Color.ShadeAlpha=0.25;
Color.Line=COLOR(2,:);
ShadedErrorbar(mratiotrialroi1(:,2)',semratiotrialroi1(:,2)',tCoor,Color); 
hold on;
plot(0:trial.nStim,ones(1+trial.nStim,1)*95,'-k','LineWidth',2);% vein
set(gca, 'ylim', [minc, maxc]);
xlabel('Time [s]')
ylabel('Relative change %')
set(gca, 'FontSize', 10);

%% contrast comparison for g1 optimal
cmin = 0.5; cmax = 2;
Fig = figure;
set(Fig, 'Position', [500 500 1200 500])
subplot(231);
imagesc(snrGG_final);caxis([cmin,cmax]); title('G1(rCBFv)'); axis image;colorbar;colormap('parula');axis off;
subplot(232)
imagesc(snrGGV_final);caxis([cmin,cmax]);title('G1(rCBV)'); axis image;colorbar;axis off;
subplot(233)
imagesc(snrGGF_final);caxis([cmin,cmax]);title('G1(rCBF)'); axis image;colorbar;axis off;
subplot(234)
imagesc(snrGGFV_final);caxis([cmin,cmax]);title('G1(optimal)'); axis image;colorbar;axis off;
subplot(235) 
imagesc(snrPDI_final);caxis([cmin,cmax]);title('PDI'); axis image;colorbar;axis off;

Fig = figure;
set(Fig, 'Position', [500 100 1200 500])
subplot(231); 
histogram(snrGG_final,50,'BinLimits',[cmin,cmax]); title(['G1(rCBFv) RMS Contrast：', num2str(std(snrGG_final, 1, "all"))]); 
subplot(232); 
histogram(snrGGV_final,50,'BinLimits',[cmin,cmax]);title(['G1(rCBV) RMS Contrast：', num2str(std(snrGGV_final, 1, "all"))]);
subplot(233)
histogram(snrGGF_final,50,'BinLimits',[cmin,cmax]);title(['G1(rCBF) RMS Contrast：', num2str(std(snrGGF_final, 1, "all"))]);
subplot(234); 
histogram(snrGGFV_final,50,'BinLimits',[cmin,cmax]);title(['G1(optimal) RMS Contrast：', num2str(std(snrGGFV_final, 1, "all"))]);
subplot(235); 
histogram(snrPDI_final,50,'BinLimits',[cmin,cmax]);title(['PDI RMS Contrast：', num2str(std(snrPDI_final, 1, "all"))]);

%%
Fig = figure;
set(Fig, 'Position', [500 500 1200 500])
subplot(231);
imagesc(snrGGdB);caxis([0,40]); title('G1(rCBFv)'); axis image;colorbar;colormap('gray');axis off;
subplot(232)
imagesc(snrGGVdB);caxis([0,40]);title('G1(rCBV)'); axis image;colorbar;axis off;
subplot(233)
imagesc(snrGGFdB);caxis([0,40]);title('G1(rCBF)'); axis image;colorbar;axis off;
subplot(234)
imagesc(snrGGFVdB);caxis([0,40]);title('G1(optimal)'); axis image;colorbar;axis off;
subplot(235) 
imagesc(snrPDIdB);caxis([0,40]);title('PDI'); axis image;colorbar;axis off;

Fig = figure;
set(Fig, 'Position', [500 100 1200 500])
subplot(231); 
histogram(snrGGdB,50,'BinLimits',[0,40]); title(['G1(rCBFv) RMS Contrast：', num2str(std(snrGGdB, 1, "all"))]); 
subplot(232); 
histogram(snrGGVdB,50,'BinLimits',[0,40]);title(['G1(rCBV) RMS Contrast：', num2str(std(snrGGVdB, 1, "all"))]);
subplot(233)
histogram(snrGGFdB,50,'BinLimits',[0,40]);title(['G1(rCBF) RMS Contrast：', num2str(std(snrGGFdB, 1, "all"))]);
subplot(234); 
histogram(snrGGFVdB,50,'BinLimits',[0,40]);title(['G1(optimal) RMS Contrast：', num2str(std(snrGGFVdB, 1, "all"))]);
subplot(235); 
histogram(snrPDIdB,50,'BinLimits',[0,40]);title(['PDI RMS Contrast：', num2str(std(snrPDIdB, 1, "all"))]);


% Fig = figure;
% set(Fig, 'Position', [500 500 1200 500])
% subplot(231);
% imagesc(snrGGdB.*(snrGGdB>17));caxis([0,40]); title('G1(rCBFv)'); axis image;colorbar;colormap('gray');axis off;
% subplot(232)
% imagesc(snrGGVdB.*(snrGGVdB>17));caxis([0,40]);title('G1(rCBV)'); axis image;colorbar;axis off;
% subplot(233)
% imagesc(snrGGFdB.*(snrGGFdB>17));caxis([0,40]);title('G1(rCBF)'); axis image;colorbar;axis off;
% subplot(234)
% imagesc(snrGGFVdB.*(snrGGFVdB>17));caxis([0,40]);title('G1(optimal)'); axis image;colorbar;axis off;
% subplot(235) 
% imagesc(snrPDIdB.*(snrPDIdB<30));caxis([0,40]);title('PDI'); axis image;colorbar;axis off;

%% R: rCBF./rCBV image
viscrop = [20, 215];
img_bk = log(mean(mtrialPDI(:,:,:),3)./eqNoise); 
img_bk = img_bk(:,viscrop(1):viscrop(2));
%actR = actR + (actR==0);
BB=ones(3,3);
BB(2,2)=2;
BB=BB/10;
actmap_covn = convn(actR(:,viscrop(1):viscrop(2)),BB,'same');
% substitute actmap_GG0, actmap_GGV0, actmap_PDI0, actmap_GGFV0
img_olap = actmap_covn;
%img_olap = actR(:,viscrop(1):viscrop(2));

cRangeBk = [min(img_bk(:)), max(img_bk(:))*0.93];
cRangeOlap = [1.5, 4];
shhold = 1.5;
fig = figure;
set(fig, 'Position', [300 300 500 400])
Fuse2Images(img_bk,img_olap, cRangeBk,cRangeOlap, shhold)
colorbar('eastoutside', 'Ticks', cRangeOlap);%,'TickLabels',{'0.3','0.9'});

