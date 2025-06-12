% ROIplot
addpath D:\g1_based_fUS\manuscript_g1fUS
numdataset = 8;
for i = 1:8

    %[ROIGG, ROIGGV, ROIGGF, ROIGGFV, ROIPDI, Ratio, Peak, T2Peak, CNR, anpratio] = g1fUS_analysis_main_updated;
    [ROIGG, ROIGGV, ROIGGF, ROIGGFV, ROIPDI, ROIVV, ROIVVcz, Ratio, Peak, T2Peak, CNR, anpratio, actmap, ImgBG] = g1fUS_analysis_main_updated;

    CNRs.GG(:,:,i) = CNR.GG(1:140,1:228); 
    CNRs.GGV(:,:,i) = CNR.GGV(1:140,1:228);
    CNRs.GGFV(:,:,i) = CNR.GGFV(1:140,1:228);
    CNRs.PDI(:,:,i) = CNR.PDI(1:140,1:228);
    CNRs.GGF(:,:,i) = CNR.GGF(1:140,1:228);
    CNRs.VV(:,:,i) = CNR.VV(1:140,1:228); 
    CNRs.VVcz(:,:,i) = CNR.VVcz(1:140,1:228);

    anpratios.GG(:,i) = anpratio.GG(1)./(anpratio.GG(1)+anpratio.GG(2))*100;
    anpratios.GGV(:,i) = anpratio.GGV(1)./(anpratio.GGV(1)+anpratio.GGV(2))*100;
    anpratios.GGFV(:,i) = anpratio.GGFV(1)./(anpratio.GGFV(1)+anpratio.GGFV(2))*100;
    anpratios.GGF(:,i) = anpratio.GGF(1)./(anpratio.GGF(1)+anpratio.GGF(2))*100;
    anpratios.PDI(:,i) = anpratio.PDI(1)./(anpratio.PDI(1)+anpratio.PDI(2))*100;
    anpratios.VV(:,i) = anpratio.VV(1)./(anpratio.VV(1)+anpratio.VV(2))*100;
    anpratios.VVcz(:,i) = anpratio.VVcz(1)./(anpratio.VVcz(1)+anpratio.VVcz(2))*100;

    actmaps.GG(:,:,i) = actmap.GG(1:140,1:228); 
    actmaps.GGV(:,:,i) = actmap.GGV(1:140,1:228);
    actmaps.GGFV(:,:,i) = actmap.GGFV(1:140,1:228);
    actmaps.PDI(:,:,i) = actmap.PDI(1:140,1:228);
    actmaps.GGF(:,:,i) = actmap.GGF(1:140,1:228);
    actmaps.VV(:,:,i) = actmap.VV(1:140,1:228);
    actmaps.VVcz(:,:,i) = actmap.VVcz(1:140,1:228);

    ImgBGs(:,:,i) = ImgBG(1:140,1:228);

    for k = 1:3
    Peaks{k}.GG.vein(:,i) = Peak{k}.GG.vein;
    Peaks{k}.GG.artery(:,i) = Peak{k}.GG.artery;
    Peaks{k}.GGV.vein(:,i) = Peak{k}.GGV.vein;
    Peaks{k}.GGV.artery(:,i) = Peak{k}.GGV.artery;
    Peaks{k}.GGF.vein(:,i) = Peak{k}.GGF.vein;
    Peaks{k}.GGF.artery(:,i) = Peak{k}.GGF.artery;

    T2Peaks{k}.GG.vein(:,i) = T2Peak{k}.GG.vein;
    T2Peaks{k}.GG.artery(:,i) = T2Peak{k}.GG.artery;
    T2Peaks{k}.GGV.vein(:,i) = T2Peak{k}.GGV.vein;
    T2Peaks{k}.GGV.artery(:,i) = T2Peak{k}.GGV.artery;
    T2Peaks{k}.GGF.vein(:,i) = T2Peak{k}.GGF.vein;
    T2Peaks{k}.GGF.artery(:,i) = T2Peak{k}.GGF.artery;

    Ratios{k}.all(:,i) = Ratio{k}.all;
    Ratios{k}.vein(:,i) = Ratio{k}.vein;
    Ratios{k}.artery(:,i) = Ratio{k}.artery;

    TC{k}.GG(:,i) = ROIGG{k}.m;
    TC{k}.GGV(:,i) = ROIGGV{k}.m;
    TC{k}.GGF(:,i) = ROIGGF{k}.m;
    TC{k}.GGFV(:,i) = ROIGGFV{k}.m;
    TC{k}.PDI(:,i) = ROIPDI{k}.m;
    TC{k}.VV(:,i) = ROIVV{k}.m;
    TC{k}.VVcz(:,i) = ROIVVcz{k}.m;

    TC{k}.vein.GG(:,i) = ROIGG{k}.vein.m;
    TC{k}.artery.GG(:,i) = ROIGG{k}.artery.m;
    TC{k}.vein.GGV(:,i) = ROIGGV{k}.vein.m;
    TC{k}.artery.GGV(:,i) = ROIGGV{k}.artery.m;
    TC{k}.vein.GGF(:,i) = ROIGGF{k}.vein.m;
    TC{k}.artery.GGF(:,i) = ROIGGF{k}.artery.m;
    end
end
save(['V:\G1based fUS data\vUS\', 'ROI_alldatasets.mat'], 'TC','Peaks','T2Peaks','anpratios');

ROIlib = {'ROI(rCBFv)','ROI(rCBV)','ROI(rCBF)'};
COLOR=[0.08 0.17 0.55;
    1 0 0;
    0.31 0.31 0.31];
trial.nBase = 25;
trial.n = 10;
trial.nStim = 5;
trial.nRecover = 25;
trial.nRest = 5; % included in nRecover
trial.nlength = trial.nStim + trial.nRecover; 

%% bar plot: v/a #of pixels rCBFv rCBV rCBF in respective ROIs

fig = figure;
set(fig, 'Position', [600 600 350 300])
datam = ([mean(1-anpratios.GG,2),mean(anpratios.GG,2);mean(1-anpratios.GGV,2),mean(anpratios.GGV,2);mean(1-anpratios.GGF,2),mean(anpratios.GGF,2)])*100;
datastd = [std(1-anpratios.GG,1,2),std(anpratios.GG,1,2);std(1-anpratios.GGV,1,2),std(anpratios.GGV,1,2);std(1-anpratios.GGF,1,2),std(anpratios.GGF,1,2)]*100;
b = bar(datam,'grouped');
set(gca,'XTickLabel', {'rCBFv', 'rCBV', 'rCBF'});
hold on
[ngroups,nbars] = size(datam);
%     errorbar
for j = 1:nbars
    x(j,:) = b(j).XEndPoints;
end
errorbar(x', datam, datastd, 'k','LineStyle','none');
hold off
legend({'Vein','Artery'})
ylabel('Fraction of pixels activated %')
set(gca, 'YLim', [0, 100])
%%
fig = figure;
set(fig, 'Position', [600 600 350 300])
datam = ([mean(1-anpratios.GGFV,2),mean(anpratios.GGFV,2);mean(1-anpratios.PDI,2),mean(anpratios.PDI,2);mean(1-anpratios.GGF,2),mean(anpratios.GGF,2)])*100;
datastd = [std(1-anpratios.GGFV,1,2),std(anpratios.GGFV,1,2);std(1-anpratios.PDI,1,2),std(anpratios.PDI,1,2);std(1-anpratios.GGF,1,2),std(anpratios.GGF,1,2)]*100;
b = bar(datam,'grouped');
set(gca,'XTickLabel', {'G1(optimal)', 'PDI', 'rCBF'});
hold on
[ngroups,nbars] = size(datam);
%     errorbar
for j = 1:nbars
    x(j,:) = b(j).XEndPoints;
end
errorbar(x', datam, datastd, 'k','LineStyle','none');
hold off
legend({'Vein','Artery'})
ylabel('Fraction of pixels activated %')
set(gca, 'YLim', [0, 100])

%% bar plot: v/a peak rCBFv rCBV rCBF in different ROIs
for k = 1: 3
    fig = figure;
    set(fig, 'Position', [600 600 350 300])
    datam = ([mean(Peaks{k}.GG.vein,2),mean(Peaks{k}.GG.artery,2);mean(Peaks{k}.GGV.vein,2),mean(Peaks{k}.GGV.artery,2);mean(Peaks{k}.GGF.vein,2),mean(Peaks{k}.GGF.artery,2)]-1)*100;
    datastd = [std(Peaks{k}.GG.vein,1,2),std(Peaks{k}.GG.artery,1,2);std(Peaks{k}.GGV.vein,1,2),std(Peaks{k}.GGV.artery,1,2);std(Peaks{k}.GGF.vein,1,2),std(Peaks{k}.GGF.artery,1,2)]*100;
    b = bar(datam,'grouped');
    set(gca,'XTickLabel', {'rCBFv', 'rCBV', 'rCBF'});
    hold on
    [ngroups,nbars] = size(datam);
%     errorbar
    for j = 1:nbars
        x(j,:) = b(j).XEndPoints;
    end
    errorbar(x', datam, datastd, 'k','LineStyle','none');
    hold off
    legend({'Vein','Artery'})
    title(ROIlib{k})
    ylabel('Peak of relative change %')
    set(gca, 'YLim', [0, 35])
end
%% bar plot: v/a time2peak rCBFv rCBV rCBF in different ROIs
for k = 1: 3
    fig = figure;
    set(fig, 'Position', [300 300 350 300])
    datam = [mean(T2Peaks{k}.GG.vein,2),mean(T2Peaks{k}.GG.artery,2);mean(T2Peaks{k}.GGV.vein,2),mean(T2Peaks{k}.GGV.artery,2);mean(T2Peaks{k}.GGF.vein,2),mean(T2Peaks{k}.GGF.artery,2)];
    datastd = [std(T2Peaks{k}.GG.vein,1,2),std(T2Peaks{k}.GG.artery,1,2);std(T2Peaks{k}.GGV.vein,1,2),std(T2Peaks{k}.GGV.artery,1,2);std(T2Peaks{k}.GGF.vein,1,2),std(T2Peaks{k}.GGF.artery,1,2)];
    b = bar(datam,'grouped');
    set(gca,'XTickLabel', {'rCBFv', 'rCBV','rCBF'});
    hold on
    [ngroups,nbars] = size(datam);
%     errorbar
    for j = 1:nbars
        x(j,:) = b(j).XEndPoints;
    end
    errorbar(x', datam, datastd, 'k','LineStyle','none');
    hold off
    legend({'Vein','Artery'})
    title(ROIlib{k})
    ylabel('Time to peak [s]')
    set(gca, 'YLim', [0, 5])
end

%% bar plot: peak ratio rCBF/rCBV in different ROIs; v/a peak ratio rCBF/rCBV in different ROIs

    fig = figure;
    set(fig, 'Position', [400 400 450 350])
    datam = [mean(Ratios{1}.all,2);mean(Ratios{2}.all,2);mean(Ratios{3}.all,2)];
    datastd = [std(Ratios{1}.all,1,2);std(Ratios{2}.all,1,2);std(Ratios{3}.all,1,2)];
    bar(datam)
    set(gca,'XTickLabel', {'ROI(rCBFv)', 'ROI(rCBV)', 'ROI(rCBF)'});
    hold on
    er= errorbar(1:3,datam,datastd, 'k','LineStyle','none');
    ylabel('Ratio of rCBF peak to rCBV peak')
 
    fig = figure;
    set(fig, 'Position', [400 400 450 350])
    datam = [mean(Ratios{1}.vein,2), mean(Ratios{1}.artery,2);mean(Ratios{2}.vein,2), mean(Ratios{2}.artery,2);mean(Ratios{3}.vein,2), mean(Ratios{3}.artery,2)];  
    datastd = [std(Ratios{1}.vein,1,2), std(Ratios{1}.artery,1,2);std(Ratios{2}.vein,1,2), std(Ratios{2}.artery,1,2);std(Ratios{3}.vein,1,2), std(Ratios{3}.artery,1,2)];
    b = bar(datam,'grouped');
    set(gca,'XTickLabel', {'ROI(rCBFv)', 'ROI(rCBV)', 'ROI(rCBF)'});
    hold on
    [ngroups,nbars] = size(datam);
%     errorbar
    for j = 1:nbars
        x(j,:) = b(j).XEndPoints;
    end
    errorbar(x', datam, datastd, 'k','LineStyle','none');
    hold off
    legend({'Vein','Artery'})
    ylabel('Ratio of rCBF peak to rCBV peak')

%% averaged Time Course: rCBFv rCBV rCBF in different ROIs
tCoor = [1:trial.nlength]-trial.nRest-1;
COLOR0 = [0.9290 0.6940 0.1250;
    0.4660 0.6740 0.1880;
    0.6350 0.0780 0.1840];
for k = 1:3
fig = figure;
set(fig, 'Position', [700 700 400 300]) 
datam = [mean(TC{k}.GG,2),mean(TC{k}.GGV,2),mean(TC{k}.GGF,2)];
datastd = [std(TC{k}.GG,1,2),std(TC{k}.GGV,1,2),std(TC{k}.GGF,1,2)]/sqrt(numdataset);
for j = 1:3
    Color.Shade = COLOR0(j,:);
    Color.ShadeAlpha=0.25;
    Color.Line=COLOR0(j,:);
    ShadedErrorbar(datam(:,j)', datastd(:,j)', tCoor, Color);
end
hold on;
plot(0:trial.nStim,ones(1+trial.nStim,1)*95,'-k','LineWidth',2);% vein
legend({'','rCBFv','','rCBV','','rCBF'})
xlabel('Time [s]')
ylabel('Relative change %')
set(gca, 'FontSize', 10);
title(ROIlib{k});
end

%% figure 2.c 
lib = {'G1(rCBFv)', 'G1(rCBV)','G1(rCBF)','G1(optimal)','PDI'} ;
for i = 1:numel(lib)
    datatype = lib{i};
switch datatype
    case 'G1(rCBFv)'
%         actmap = actmap_GG0;
        ROI = TC{1}.GG;
    case 'G1(rCBV)'
%         actmap = actmap_GGV0;
        ROI = TC{2}.GGV;
    case 'G1(rCBF)'
%         actmap = actmap_GGF0;
        ROI = TC{3}.GGF;
    case 'G1(optimal)'
%         actmap = actmap_GGFV0;
        ROI = TC{3}.GGFV;
    case 'PDI'
%         actmap = actmap_PDI0;
        ROI = TC{3}.PDI;
end

ratiotrialroi = ROI;
mratiotrialroi = mean(ROI,2);
semratiotrialroi = std(ROI,1,2)/sqrt(numdataset);
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
% yyaxis right; ax = gca; plot(tCoor, smooth(stimhrf,4)*median(mratiotrialroi)/10+100,'Color','r'); 
% ax.YColor = 'w'; set(gca, 'ylim',[minc,maxc], 'YTick', [], 'YTickLabel',[]);
% hold on;
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
% yyaxis right; ax = gca; plot(tCoor, smooth(stimhrf,4)*median(mratiotrialroi)/10+100,'Color','r'); 
% ax.YColor = 'w'; set(gca, 'ylim',[minc,maxc], 'YTick', [], 'YTickLabel',[]);
% hold on;
plot(0:trial.nStim,ones(1+trial.nStim,1)*95,'-k','LineWidth',2);
title(datatype);
end

%% averaged Time Course: v/a rCBFv in different ROIs
tCoor = [1:trial.nlength]-trial.nRest-1;
% COLOR = ['r','b','k'];
for k = 1:3
fig = figure;
set(fig, 'Position', [500 500 400 300])
datam = [mean(TC{k}.vein.GG,2),mean(TC{k}.artery.GG,2)];
datastd = [std(TC{k}.vein.GG,1,2),std(TC{k}.artery.GG,1,2)]/sqrt(numdataset);
for j = 1:2
    Color.Shade = COLOR(j,:);
    Color.ShadeAlpha=0.25;
    Color.Line=COLOR(j,:);
    ShadedErrorbar(datam(:,j)', datastd(:,j)', tCoor, Color);
end
hold on;
plot(0:trial.nStim,ones(1+trial.nStim,1)*95,'-k','LineWidth',2);% vein
legend({'','Vein','','Artery'})
xlabel('Time [s]')
ylabel('rCBFv [%]')
set(gca, 'FontSize', 10);
title(ROIlib{k});
set(gca, 'YLim', [95, 125])
end
%% averaged Time Course: v/a rCBV in different ROIs
tCoor = [1:trial.nlength]-trial.nRest-1;
% COLOR = ['r','b','k'];
for k = 1:3
fig = figure;
set(fig, 'Position', [500 500 400 300])
datam = [mean(TC{k}.vein.GGV,2),mean(TC{k}.artery.GGV,2)];
datastd = [std(TC{k}.vein.GGV,1,2),std(TC{k}.artery.GGV,1,2)]/sqrt(numdataset);
for j = 1:2
    Color.Shade = COLOR(j,:);
    Color.ShadeAlpha=0.25;
    Color.Line=COLOR(j,:);
    ShadedErrorbar(datam(:,j)', datastd(:,j)', tCoor, Color);
end
hold on;
plot(0:trial.nStim,ones(1+trial.nStim,1)*95,'-k','LineWidth',2);% vein
legend({'','Vein','','Artery'})
xlabel('Timee [s]')
ylabel('rCBV [%]')
set(gca, 'FontSiz', 10);
title(ROIlib{k});
set(gca, 'YLim', [95, 125])
end
%% averaged Time Course: v/a rCBF in different ROIs
tCoor = [1:trial.nlength]-trial.nRest-1;
% COLOR = ['r','b','k'];
for k = 1:3
fig = figure;
set(fig, 'Position', [500 500 400 300])
datam = [mean(TC{k}.vein.GGF,2),mean(TC{k}.artery.GGF,2)];
datastd = [std(TC{k}.vein.GGF,1,2),std(TC{k}.artery.GGF,1,2)]/sqrt(numdataset);
for j = 1:2
    Color.Shade = COLOR(j,:);
    Color.ShadeAlpha=0.25;
    Color.Line=COLOR(j,:);
    ShadedErrorbar(datam(:,j)', datastd(:,j)', tCoor, Color);
end
hold on;
plot(0:trial.nStim,ones(1+trial.nStim,1)*95,'-k','LineWidth',2);% vein
legend({'','Vein','','Artery'})
xlabel('Time [s]')
ylabel('rCBF [%]')
set(gca, 'FontSize', 10);
title(ROIlib{k});
set(gca, 'YLim', [95, 130])
end

%% CNR
% CNRs.GG(:,:,4) = fliplr(CNRs.GG(:,:,4));
% CNRs.GGV(:,:,4) = fliplr(CNRs.GGV(:,:,4));
% CNRs.GGFV(:,:,4) = fliplr(CNRs.GGFV(:,:,4));
% CNRs.PDI(:,:,4) = fliplr(CNRs.PDI(:,:,4));
% CNRs.GGF(:,:,4) = fliplr(CNRs.GGF(:,:,4));

Fig = figure;
set(Fig, 'Position', [500 500 900 450])
subplot(231);
imagesc(mean(CNRs.GG,3)); title('G1(rCBFv)'); axis image; caxis([0.5,2]);colorbar;axis off;
subplot(232)
imagesc(mean(CNRs.GGV,3));title('G1(rCBV)'); axis image; caxis([0.5,2]);colorbar;axis off;
subplot(233)
imagesc(mean(CNRs.GGF,3));title('G1(rCBF)'); axis image; caxis([0.5,2]);colorbar;axis off;
subplot(234)
imagesc(mean(CNRs.GGFV,3));title('G1(optimal)'); axis image; caxis([0.5,2]);colorbar;axis off;
subplot(235)
imagesc(mean(CNRs.PDI,3));title('PDI'); axis image; caxis([0.5,2]);colorbar;axis off;