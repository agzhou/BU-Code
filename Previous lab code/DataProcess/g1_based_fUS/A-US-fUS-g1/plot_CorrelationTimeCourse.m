%function plot_CorrelationTimeCourse(datatype)
lib = {'G1(rCBFv)', 'G1(rCBV)','G1(rCBF)','G1(optimal)','PDI'} ;
for i = 1:numel(lib)
    datatype = lib{i};
switch datatype
    case 'G1(rCBFv)'
        actmap = actmap_GG0;
        ROI = ROIGG;
    case 'G1(rCBV)'
        actmap = actmap_GGV0;
        ROI = ROIGGV;
    case 'G1(rCBF)'
        actmap = actmap_GGF0;
        ROI = ROIGGF;
    case 'G1(optimal)'
        actmap = actmap_GGFV0;
        ROI = ROIGGFV;
    case 'PDI'
        actmap = actmap_PDI0;
        ROI = ROIPDI;
end
 actmap_GG0 = max(actmap_GG,[],3);actmap_GGV0 = max(actmap_GGV,[],3);actmap_GGF0 = max(actmap_GGF,[],3);
actmap_GGFV0 = max(actmap_GGFV,[],3);actmap_PDI0 = max(actmap_PDI,[],3);
%%
COLOR=[1 0 0;
    0.08 0.17 0.55
    0.31 0.31 0.31];
addpath D:\g1_based_fUS\A-US-fUS-g1\SubFunctions
% load actmap_GG0, actmap_GGV0, actmap_PDI0, actmap_GGFV0
% load coefmap_GG0, 
% load dir

%% 1. masked Correlation map: CBFv based fUS
viscrop = [20, 215];
% img_bk = log(mean(mtrialPDI(:,:,:),3)./eqNoise); 
img_bk = (mean(mtrialPDI(:,:,:),3)./eqNoise).^0.35;
img_bk = img_bk(:,viscrop(1):viscrop(2));
BB=ones(3,3);
BB(2,2)=2;
BB=BB/10;
actmap_covn = convn(actmap(:,viscrop(1):viscrop(2)),BB,'same');
% substitute actmap_GG0, actmap_GGV0, actmap_PDI0, actmap_GGFV0
img_olap = actmap_covn;
% img_olap = actmap_GG0(:,viscrop(1):viscrop(2));

cRangeBk = [min(img_bk(:)), max(img_bk(:))*0.6];
cRangeOlap = [0.1, 0.9];
shhold = 0.15;
fig = figure;
set(fig, 'Position', [300 300 500 400])
Fuse2Images(img_bk,img_olap, cRangeBk,cRangeOlap, shhold)
colorbar('eastoutside', 'Ticks', cRangeOlap);%,'TickLabels',{'0.3','0.9'});
title(datatype);

saveas(fig, [FilePath, 'Corr_', datatype,'.png']); 

%% 2. roi based time course
ratiotrialroi = ROI.ratio;
mratiotrialroi = ROI.m;
semratiotrialroi = ROI.sem;
minc = 90;
maxc = 130;

tCoor = [1:trial.nlength]-trial.nRest-1;
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
plot(0:trial.nStim,ones(1+trial.nStim,1)*95,'-k','LineWidth',2);

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
plot(0:trial.nStim,ones(1+trial.nStim,1)*95,'-k','LineWidth',2);
title(datatype);

saveas(fig, [FilePath, 'TimeC_', datatype,'.png']); 

end
