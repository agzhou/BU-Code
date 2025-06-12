%function plot_avCorrelationTimeCourse(datatype)
lib = {'G1(rCBFv)', 'G1(rCBV)','G1(rCBF)','G1(optimal)','PDI'} ;
for i = 1:numel(lib)
    datatype = lib{i};
switch datatype
    case 'G1(rCBFv)'
        arteriesAmap = arteriesAmap_GG0;
        veinAmap = veinAmap_GG0;
        ROI = ROIGG;
    case 'G1(rCBV)'
        arteriesAmap = arteriesAmap_GGV0;
        veinAmap = veinAmap_GGV0;
        ROI = ROIGGV;
    case 'G1(rCBF)'
        arteriesAmap = arteriesAmap_GGF0;
        veinAmap = veinAmap_GGF0;
        ROI = ROIGGF;
    case 'G1(optimal)'
        arteriesAmap = arteriesAmap_GGFV0;
        veinAmap = veinAmap_GGFV0;
        ROI = ROIGGFV;
    case 'PDI'
        arteriesAmap = arteriesAmap_PDI0;
        veinAmap = veinAmap_PDI0;
        ROI = ROIPDI;
end

%% 3. arteries/veins masked Correlation map: CBFv based fUS
viscrop = [20, 215];
img_bk = log(mean(mtrialPDI(:,:,:),3)./eqNoise); 
img_bk = img_bk(:,viscrop(1):viscrop(2));
BB=ones(3,3);
BB(2,2)=2;
BB=BB/10;
actmap_covn1 = convn(arteriesAmap(:,viscrop(1):viscrop(2)),BB,'same');
actmap_covn2 = convn(veinAmap(:,viscrop(1):viscrop(2)),BB,'same');
% substitute actmap_GG0, actmap_GGV0, actmap_PDI0, actmap_GGFV0

img_olap1 = actmap_covn1;
img_olap2 = actmap_covn2;
cRangeBk = [min(img_bk(:)), max(img_bk(:))*0.93];
cRangeOlap1 = [0.1, 0.9];
cRangeOlap2 = [0.1, 0.9];
shhold = 0.15;
fig = figure;
set(fig, 'Position', [300 300 500 400])
Fuse3Images(img_bk,img_olap1, img_olap2, cRangeBk,cRangeOlap1,cRangeOlap2, shhold, 0)
colorbar('eastoutside', 'Ticks', cRangeOlap2);%,'TickLabels',{'0.3','0.9'});
title(datatype);

saveas(fig, [FilePath, 'avCorr_', datatype,'.png']); 

%% 5. arteries/veins roi based time course
ratiotrialroi = cat(3, ROI.artery.ratio, ROI.vein.ratio);
mratiotrialroi = [ROI.artery.m, ROI.vein.m];
semratiotrialroi = [ROI.artery.sem, ROI.vein.sem];
minc = 90;
maxc = 130;

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
title(datatype);

%%
saveas(fig, [FilePath, 'avTimeC_', datatype,'.png']); 

end