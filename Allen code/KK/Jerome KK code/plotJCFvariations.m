CPP% Plot JCF variation data

%% Load Data
clearvars
close all

load JCFvariations2.mat

%% Plot results
figure
tiledlayout(2,4)
nexttile(1,[2,1])
plotGammaScaleImage(pS.xCoord*1e3,pS.zCoord*1e3,JCF1,0.1)
axis image
title('Raised to 1')
nexttile(2,[2,1])
plotGammaScaleImage(pS.xCoord*1e3,pS.zCoord*1e3,JCF2,0.1)
axis image
title('Raised to 2')
nexttile(3,[2,1])
plotGammaScaleImage(pS.xCoord*1e3,pS.zCoord*1e3,JCF3,0.1)
axis image
title('Raised to 3')
nexttile(4,[2,1])
plotGammaScaleImage(pS.xCoord*1e3,pS.zCoord*1e3,JCF4,0.1)
axis image
title('Raised to 4')

figure
tiledlayout(2,4)
nexttile(1,[2,1])
plotGammaScaleImage(pS.xCoord*1e3,pS.zCoord*1e3,JCF5,0.1)
axis image
title('No Conj Raised to 1')
nexttile(2,[2,1])
plotGammaScaleImage(pS.xCoord*1e3,pS.zCoord*1e3,JCF6,0.1)
axis image
title('No Conj Raised to 2')
nexttile(3,[2,1])
plotGammaScaleImage(pS.xCoord*1e3,pS.zCoord*1e3,JCF7,0.1)
axis image
title('No Conj Raised to 3')
nexttile(4,[2,1])
plotGammaScaleImage(pS.xCoord*1e3,pS.zCoord*1e3,JCF8,0.1)
axis image
title('No Conj Raised to 4')
figure
plotGammaScaleImage(pS.xCoord*1e3,pS.zCoord*1e3,JCF9,0.1)
axis image
title('Default CPP')