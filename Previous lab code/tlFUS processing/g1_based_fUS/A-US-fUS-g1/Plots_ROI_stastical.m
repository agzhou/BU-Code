
R = (max(movmean(ROIGG.m, 3)/100)*max(movmean(ROIGGV.m, 3)/100)-1)/(max(movmean(ROIGGV.m, 3)/100)-1);
disp(['Ratio of Peak rCBF / Peak rCBV: ', num2str(R)])

Fig = figure; set(Fig,'Position',[400 600 450 300])
plot(ROIGG.m); 
hold on; plot(ROIGGV.m);
hold on; plot(ROIPDI.m);
hold on; plot(ROIGGFV.m);
hold on; plot(ROIGGF.m);
hold on; 
plot(trial.nRest+1:trial.nRest+trial.nStim,ones(trial.nStim,1)*min(ROIPDI.m(:)),'-k','LineWidth',2);
legend({'g1-CBFv','g1-CBV','pdi','novel','g1-CBF'})
title(['rCBFlow/rCBV = ', num2str(R)])
xlabel('Time[s]');
ylabel('%');
ylim([90 130]);


[vpGG,~] = max(movmean(ROIGG.vein.m(window(1):window(2)), 3)/100);
[vpGGV,~] = max(movmean(ROIGGV.vein.m(window(1):window(2)), 3)/100);
disp(['Peak rCBFv and Peak rCBV in veins: ', num2str(vpGG), ' ', num2str(vpGGV)])

[apGG,~] = max(movmean(ROIGG.artery.m(window(1):window(2)), 3)/100);
[apGGV,~] = max(movmean(ROIGGV.artery.m(window(1):window(2)), 3)/100);
disp(['Peak rCBFv and Peak rCBV in arteries: ', num2str(apGG), ' ', num2str(apGGV)])

interpx = interp1(1:trial.nlength, ROIGG.vein.m, 1:0.2:trial.nlength, "spline");
[~,vtpGG0] = max(interpx(window(1)*5:window(2)*5));
vtpGG = vtpGG0/5;
interpx = interp1(1:trial.nlength, ROIGGV.vein.m, 1:0.2:trial.nlength, "spline");
[~,vtpGGV0] = max(interpx(window(1)*5:window(2)*5));
vtpGGV = vtpGGV0/5;
disp(['Time to Peak rCBFv and Time to Peak rCBV in veins: ', num2str(vtpGG), ' ', num2str(vtpGGV)])

interpx = interp1(1:trial.nlength, ROIGG.artery.m, 1:0.2:trial.nlength, "spline");
[~,atpGG0] = max(interpx(window(1)*5:window(2)*5));
atpGG = atpGG0/5;
interpx = interp1(1:trial.nlength, ROIGGV.artery.m, 1:0.2:trial.nlength, "spline");
[~,atpGGV0] = max(interpx(window(1)*5:window(2)*5));
atpGGV = atpGGV0/5;
disp(['Time to Peak rCBFv and Time to Peak rCBV in arteries: ', num2str(atpGG), ' ', num2str(atpGGV)])

figure; plot(ROIGG.vein.m); hold on;  plot(ROIGG.artery.m);
figure; plot(ROIGGV.vein.m); hold on;  plot(ROIGGV.artery.m);

% switch ROItype
%     case 'rCBFv_based'
% 
% ROI_Fv.GG = [ROI_Fv.GG,ROIGG.m];
% ROI_Fv.GGV = [ROI_Fv.GGV,ROIGGV.m];
% ROI_Fv.GGF = [ROI_Fv.GGF,ROIGGF.m];
% ROI_Fv.GGFV = [ROI_Fv.GGFV,ROIGGFV.m];
% ROI_Fv.PDI = [ROI_Fv.PDI,ROIPDI.m];
% 
% 
% Rs_Fv = [Rs_Fv, R];
% 
% apGGs_Fv = [apGGs_Fv, apGG];
% apGGVs_Fv = [apGGVs_Fv, apGGV];
% vpGGs_Fv = [vpGGs_Fv, vpGG];
% vpGGVs_Fv = [vpGGVs_Fv, vpGGV];
% 
% atpGGs_Fv = [atpGGs_Fv, atpGG];
% atpGGVs_Fv = [atpGGVs_Fv, atpGGV];
% vtpGGs_Fv = [vtpGGs_Fv, vtpGG];
% vtpGGVs_Fv = [vtpGGVs_Fv, vtpGGV];
% 
%     case 'rCBV_based'
% ROI_V.GG = [ROI_V.GG, ROIGG.m];
% ROI_V.GGV = [ROI_V.GGV,ROIGGV.m];
% ROI_V.GGF = [ROI_V.GGF,ROIGGF.m];
% 
% ROI_V.GGFV = [ROI_V.GGFV,ROIGGFV.m];
% ROI_V.PDI = [ROI_V.PDI,ROIPDI.m];
% 
% Rs_V = [Rs_V, R];
% 
% apGGs_V = [apGGs_V, apGG];
% apGGVs_V = [apGGVs_V, apGGV];
% vpGGs_V = [vpGGs_V, vpGG];
% vpGGVs_V = [vpGGVs_V, vpGGV];
% 
% atpGGs_V = [atpGGs_V, atpGG];
% atpGGVs_V = [atpGGVs_V, atpGGV];
% vtpGGs_V = [vtpGGs_V, vtpGG];
% vtpGGVs_V = [vtpGGVs_V, vtpGGV];
%  case 'rCBF_based'
%      ROI_F.GG = [ROI_F.GG,ROIGG.m];
% ROI_F.GGV = [ROI_F.GGV,ROIGGV.m];
% ROI_F.GGF = [ROI_F.GGF,ROIGGF.m];
% 
% ROI_F.GGFV = [ROI_F.GGFV,ROIGGFV.m];
% ROI_F.PDI = [ROI_F.PDI,ROIPDI.m];
% 
% Rs_F = [Rs_F, R];
% 
% apGGs_F = [apGGs_F, apGG];
% apGGVs_F = [apGGVs_F, apGGV];
% vpGGs_F = [vpGGs_F, vpGG];
% vpGGVs_F = [vpGGVs_F, vpGGV];
% 
% atpGGs_F = [atpGGs_F, atpGG];
% atpGGVs_F = [atpGGVs_F, atpGGV];
% vtpGGs_F = [vtpGGs_F, vtpGG];
% vtpGGVs_F = [atpGGVs_F, vtpGGV];
% end
