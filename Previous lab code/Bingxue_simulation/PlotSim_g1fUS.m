%% CBFv vs vtot agl diameter
%%%%%%%%%%%%%%%%
% addpath D:\g1_based_fUS\NumericalSimulation\
% addpath D:\g1_based_fUS\A-US-fUS-g1\SubFunctions\
%%
% Vtot = -[3,6,9,12,15,18]*1e-3;%,15,18]*1e-3;
Vtot = -0.03;
% agl = [0,45,90];
agl = 0;
% WidthVes = [10,30,50,70,90];%um
WidthVes = 30;
nRpt = 5;
snr = 5;
WidthVes0 = 50;
%%
CBFvs = zeros(length(Vtot),length(agl), nRpt);
CBVs = zeros(length(Vtot),length(agl), nRpt);
PDIs = zeros(length(Vtot),length(agl), nRpt);
nparticles = zeros(length(Vtot),length(agl), nRpt);
%% generate CBFv index data
% tic
%get IQ0 for calculating noise
[IQ0, FIQ, g1, Numer, tCoor, fCoor, tauCoor, Noise, nparticle] = Numerical_g1fUS_bx(Vtot(1), WidthVes0, snr, 0);
for j = 1:length(agl)
for k = 1: length(Vtot)
    disp(num2str(k));
for i = 1:nRpt
    %generate noise with snr level each run
    IQn = awgn(squeeze(IQ0),snr,'measured','db');
    Noise = IQn - squeeze(IQ0);
    %calculation
    [IQ, FIQ, g1, Numer, tCoor, fCoor, tauCoor, Noise, nparticle] = Numerical_g1fUS_bx(Vtot(k), WidthVes0, Noise, agl(j));
    absg1 = abs(squeeze(g1));
    CBFvs(k,j,i) = sqrt(abs(log(absg1(1,:)./absg1(15,:))));
    CBVs(k,j,i) = ((absg1(1))./(1-(absg1(1))));
    PDIs(k,j,i) = sum(abs(squeeze(IQ)).^2);
    nparticles(k,j,i) = nparticle;
    g1s(k,j,i,:) = squeeze(g1);
    IQs(k,j,i,:) = squeeze(IQ);
end
end
end
% toc

%%
mCBFv = mean(CBFvs,3); stdCBFv = std(CBFvs,1,3);
mCBV = mean(CBVs,3); stdCBV = std(CBVs,1,3);
mnptcl = mean(nparticles,3); stdnptcl = std(nparticles,1,3);
mPDI = mean(PDIs,3); stdPDI = std(PDIs,1,3);
%%
figure; errorbar(repmat(Vtot/Vtot(1),[3,1])', mCBFv./mCBFv(1,:), stdCBFv/mCBFv(1), 'LineWidth',2); 
hold on; plot(Vtot/Vtot(1),Vtot/Vtot(1),'-.k');
xlabel('True rCBFv');
ylabel('Measured rCBFv');
legend('0^{o}','30^{o}','60^{o}','90^{o}')

figure; errorbar(mnptcl/mnptcl(1), mCBV/mCBV(1), stdCBV/mCBV(1), 'LineWidth',2); 
hold on; plot(mnptcl/mnptcl(1), mnptcl/mnptcl(1), '-.k');
xlabel('True rCBV');
ylabel('Measured rCBV');
legend('0^{o}','30^{o}','60^{o}','90^{o}')
figure; errorbar(mnptcl/mnptcl(1), mPDI/mPDI(1), stdPDI/mPDI(1), 'LineWidth',2);
hold on; plot(mnptcl/mnptcl(1), mnptcl/mnptcl(1), '-.k');
xlabel('True rCBV');
ylabel('Measured rCBV(PDI)');
legend('0^{o}','45^{o}','90^{o}')

%%
fig = figure(14); set(fig,'Position',[800 400 800 300]);
hold on;
for j = 1: (length(agl))
rVP = abs(Vtot)*1e3;
[p11, S] = polyfit(rVP(1:end)',mCBFv(1:end,j),1);
y11 = polyval(p11,rVP(1:end));
subplot(121); hold on;
errorbar(rVP(1:end), mCBFv(1:end,j), stdCBFv(1:end,j),'.','MarkerSize',10,'LineWidth',1); 
hold on; box off; %axis equal tight;
xlabel(['Velocity [mm/s]']);
ylabel(['CBFv_{index} [a.u.]']); 
%h = legend({'','Fit: k = 1.007, b = 0.039'});set(h,'box','off');
% title(['k=',num2str(p11(1)),' b=',num2str(p11(2))])
legend('0^{\circ}','45^{\circ}','90^{\circ}')
% legend('10\mum','30\mum','50\mum','70\mum','90\mum')
fitline(j,:) = p11;
ylim([0,0.8])
end
T_rVP = []; T_mCBFv = [];
jend = [6,6,6];
for j = 1:length(agl)
rVP = abs(Vtot)*1e3;
T_rVP = [T_rVP,rVP(1:jend(j))];
T_mCBFv  = [T_mCBFv;mCBFv(1:jend(j),j)];
end
[p11, S] = polyfit(T_rVP',T_mCBFv,1);
y11 = polyval(p11,[0,rVP]);
% figure(30); 
hold on;
plot([0,rVP], y11,'--k','LineWidth',1.2);ylim([0,0.8])

%%
Color = [0 0.4470 0.7410;
    0.8500 0.3250 0.0980;
    0.9290 0.6940 0.1250];
fig = figure(14); set(fig,'Position',[600 200 400 300]);

for j = 1: (length(agl))
rVP = abs(Vtot)*1e3;
[p11, S] = polyfit(rVP(1:end)',mCBFv0(1:end,j),1);
y11 = polyval(p11,[0,rVP(1:end)]);
hold on;
errorbar(rVP(1:end), mCBFv0(1:end,j), stdCBFv0(1:end,j),'.','MarkerSize',10,'LineWidth',1,'Color',Color(j,:)); 
hold on; box off; %axis equal tight;
xlabel(['Velocity [mm/s]']);
ylabel(['CBFv_{index} [a.u.]']); 
% h = legend({'','Fit: k = 1.007, b = 0.039'});set(h,'box','off');
% title(['k=',num2str(p11(1)),' b=',num2str(p11(2))])

% legend('10\mum','30\mum','50\mum','70\mum','90\mum')
fitline(j,:) = p11;
ylim([0,0.8])

hold on;
plot([0,rVP], y11,'-.','LineWidth',1.2,'Color',Color(j,:));

end
legend({'0^{\circ}','Fit','45^{\circ}','Fit','90^{\circ}','Fit'})%:k=0.042,b=2.1e-3
T_rVP = []; T_mCBFv = [];
jend = [6,6,6];
for j = 1:length(agl)
rVP = abs(Vtot)*1e3;
T_rVP = [T_rVP,rVP(1:jend(j))];
T_mCBFv  = [T_mCBFv;mCBFv0(1:jend(j),j)];
end
% [p11, S] = polyfit(T_rVP',T_mCBFv,1);
% y11 = polyval(p11,[0,rVP]);
% figure(30); hold on;
% plot([0,rVP], y11,'-.k','LineWidth',1.2);ylim([0,1.5])


%% CBFv index vs. diameter
%%%%%%%%%%%%%%%%%%%%%%%
fig = figure; set(fig,'Position',[600 200 400 300]);
for j = 1%:(length(agl)+2)
rVP = WidthVes;%abs(Vtot)*1e3;
[p11, S] = polyfit(rVP(1:end),mCBFv(3,1:end),1);
y11 = polyval(p11,rVP(1:end));
 hold on;
errorbar(rVP(1:end), mCBFv(3,1:end), stdCBFv(3,1:end),'MarkerSize',10,'LineWidth',1.5,'Color','b'); 
hold on; box off; %axis equal tight;
xlabel(['Diameter [\mum]']);
ylabel(['CBF-speed_{index} [a.u.]']); 
%h = legend({'','Fit: k = 1.007, b = 0.039'});set(h,'box','off');
% title(['k=',num2str(p11(1)),' b=',num2str(p11(2))])
% legend('0^{\circ}','45^{\circ}','90^{\circ}')
legend('Velocity: 9 mm/s. Angle: 45^{\circ}')
fitline(j,:) = p11;
ylim([0.2,1])
end
%% CBV index vs. velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = figure; set(fig,'Position',[600 200 400 300]);
for j = 1%:(length(agl)+2)
rVP = abs(Vtot)*1e3;
[p11, S] = polyfit(rVP(1:end)',mCBV(1:end,2),1);
y11 = polyval(p11,rVP(1:end));
hold on;
errorbar(rVP(1:end), mCBV(1:end,1)+0.5, stdCBV(1:end,1),'MarkerSize',10,'LineWidth',1.5,'Color','b'); 
hold on; box off; %axis equal tight;
xlabel(['Velocity [mm/s]']);
ylabel(['CBV_{index} [a.u.]']); 
%h = legend({'','Fit: k = 1.007, b = 0.039'});set(h,'box','off');
% title(['k=',num2str(p11(1)),' b=',num2str(p11(2))])
% legend('0^{\circ}','45^{\circ}','90^{\circ}')
legend('Diameter: 50\mum. Angle: 45^{\circ}')
fitline(j,:) = p11;
ylim([0.2,1])
end
%%
nt = 20000; nTau = 100;
for j = 1:length(agl)
for k = 1: length(Vtot)
    disp(num2str(k));
for i = 1:nRpt
%     IQ0(1,1,:) = IQs(k,j,i,:);
%     IQ0(1,1,:) = IQ0 - mean(IQ0,'all');
%     [g1, Numer] = IQ2g1(IQ0,1,nt,nTau);
%     g1s0(k,j,i,:) = g1;
%     absg1 = abs(squeeze(g1));
    g1 = g1s0(k,j,i,:);
    absg1 = abs(squeeze(g1));
    CBFvs0(k,j,i) = sqrt(abs(log(absg1(1,:)./absg1(10,:))));
    CBVs0(k,j,i) = ((absg1(1))./(1-(absg1(1))));
%     PDIs0(k,j,i) = sum(abs(squeeze(IQ)).^2);
end
end
end


%% CBV vs. diameter/psf
%%%%%%%%%%%%%%%%% generate CBV data
%%
rV = [10,20,30,40,50,60,70,80,90,100,110,120];%
Vtot0 = -10e-3;
WidthVes0 = 10;%24
nRpt = 5;
snr = 1;
agl = [0,15,30,45,60,75,90];
psf = [75,100,125];
%%
CBVs = zeros(length(rV),3, nRpt);
PDIs = zeros(length(rV),3, nRpt);
nparticles = zeros(length(rV),3, nRpt);

% tic
%get IQ0 for calculating noise
[IQ0, FIQ, g1, Numer, tCoor, fCoor, tauCoor, Noise, nparticle] = Numerical_g1fUS_bx(100, WidthVes0, snr);
for j = 1:length(psf)
for k = 1: length(rV)
    disp(num2str(k));
for i = 1:nRpt
    %generate noise with snr level each run
    IQn = awgn(squeeze(IQ0),snr,'measured','db');
    Noise = IQn - squeeze(IQ0);
    %calculation
    [IQ, FIQ, g1, Numer, tCoor, fCoor, tauCoor, Noise, nparticle] = Numerical_g1fUS_bx(psf(j), rV(k), Noise);
    absg1 = abs(squeeze(g1));
    CBVs(k,j,i) = ((absg1(1))./(1-(absg1(1))));
    PDIs(k,j,i) = sum(abs(squeeze(IQ)).^2);
    nparticles(k,j,i) = nparticle;
end
end
end
% toc

%%
mCBV = mean(CBVs,2); stdCBV = std(CBVs,1,2);
mnptcl = mean(nparticles,2); stdnptcl = std(nparticles,1,2);
mPDI = mean(PDIs,2); stdPDI = std(PDIs,1,2);

figure; errorbar(mnptcl/mnptcl(1), mCBV/mCBV(1), stdCBV/mCBV(1),'b', 'LineWidth',2); 
hold on; plot(mnptcl/mnptcl(1), mnptcl/mnptcl(1), '-.k');
xlabel('True rCBV');
ylabel('Measured rCBV');
figure; errorbar(mnptcl/mnptcl(1), mPDI/mPDI(1), stdPDI/mPDI(1),'b', 'LineWidth',2);
hold on; plot(mnptcl/mnptcl(1), mnptcl/mnptcl(1), '-.k');
xlabel('True rCBV');
ylabel('Measured rCBV(PDI)');

%%
rVP = rV/75; 
[p11, S] = polyfit(rVP(1:end-8),mCBV(1:end-8),1);
y11 = polyval(p11,rVP(1:end));
fig = figure; set(fig,'Position',[800 400 800 300]);
subplot(121);
errorbar(rVP(1:end), mCBV(1:end), stdCBV(1:end), '.b','MarkerSize',10,'LineWidth',1); hold on; plot(rVP(1:end), y11,'-b','LineWidth',0.75)
hold on; box off; %axis equal tight;
xlabel(['Diameter/PSF(fixed)']);
ylabel(['CBV_{index}']); xlim([0,1])
%h = legend({'','Fit: k = 1.007, b = 0.039'});set(h,'box','off');
title(['k=',num2str(p11(1)),' b=',num2str(p11(2))])


%% CBV index plots
%%%%%%%% diameter/psf different psf
%% load data
CBVs_ = NaN(12,3,5); nparticles_ = NaN(12,3,5); PDIs_ = NaN(12,3,5);
load(['D:\g1_based_fUS\NumericalSimulation\','agl0wdth10-10-100psf75.mat'])
CBVs_(1:10,1,:) = CBVs; nparticles_(1:10,1,:) = nparticles; PDIs_(1:10,1,:) = PDIs;
load(['D:\g1_based_fUS\NumericalSimulation\','agl0wdth10-10-100psf100.mat'])
CBVs_(1:10,2,:) = CBVs; nparticles_(1:10,2,:) = nparticles; PDIs_(1:10,2,:) = PDIs;
load(['D:\g1_based_fUS\NumericalSimulation\','agl0wdth10-10-120psf125.mat'])
CBVs_(:,3,:) = CBVs; nparticles_(:,3,:) = nparticles; PDIs_(:,3,:) = PDIs;
clear CBVs nparticles PDIs
CBVs = CBVs_; nparticles = nparticles_; PDIs = PDIs_;
clear CBVs_ nparticles_ PDIs_
%%
mCBV = mean(CBVs,3); stdCBV = std(CBVs,1,3);
mnptcl = mean(nparticles,3); stdnptcl = std(nparticles,1,3);
mPDI = mean(PDIs,3); stdPDI = std(PDIs,1,3);

figure; errorbar(mnptcl./mnptcl(1,:), mCBV./mCBV(1,:), stdCBV./mCBV(1,:), 'LineWidth',2); 
hold on; plot(mnptcl./mnptcl(1,:), mnptcl./mnptcl(1,:), '-.k');
xlabel('True rCBV');
ylabel('Measured rCBV');
legend('PSF=75\mum','PSF=100\mum','PSF=125\mum')
figure; errorbar(mnptcl./mnptcl(1,:), mPDI./mPDI(1,:), stdPDI./mPDI(1,:), 'LineWidth',2);
hold on; plot(mnptcl./mnptcl(1,:), mnptcl./mnptcl(1,:), '-.k');
xlabel('True rCBV');
ylabel('Measured rCBV(PDI)');
legend('PSF=75\mum','PSF=100\mum','PSF=125\mum')
%%
psf = [75,100,125]; rV = 10:10:120;
jend = [4,6,6];
nptcl = [0.75,1,1.25].^1; %mean(mnptcl(1:12,:),1)  [0.7389, 0.9907, 1.2933];
fig = figure; set(fig,'Position',[600 200 400 300]);
for j = 1: length(psf)
rVP = rV/psf(j);
[p11, S] = polyfit(rVP(1:jend(j)),mCBV(1:jend(j),j)/nptcl(j),1);
y11 = polyval(p11,[0,rVP(1:end)]);
hold on;
errorbar(rVP(1:end), mCBV(1:end,j)/nptcl(j), stdCBV(1:end,j)/nptcl(j),'MarkerSize',10,'LineWidth',1); %hold on; plot([0,rVP(1:end)], y11,'-.','LineWidth',0.75) 
hold on; box off; %axis equal tight;
xlabel(['Diameter/PSF ']);%
ylabel(['Normalized CBV_{index} [a.u.]']);xlim([0,1])
%h = legend({'','Fit: k = 1.007, b = 0.039'});set(h,'box','off');
% title(['k=',num2str(p11(1)),' b=',num2str(p11(2))])
legend('\sigma_x=75\mum','\sigma_x=100\mum','\sigma_x=125\mum')
fitline(j,:) = p11;
end
%%
T_rVP = []; T_mCBV = [];
jend = [4,6,7];
for j = 1:length(psf)
rVP = rV/psf(j);
T_rVP = [T_rVP,rVP(1:jend(j))];
T_mCBV  = [T_mCBV;mCBV(1:jend(j),j)/nptcl(j)];
end
[p11, S] = polyfit(T_rVP,T_mCBV,1);
y11 = polyval(p11,[0,rVP]);
figure(9); hold on;
plot([0,rVP], y11,'-.k','LineWidth',1.2)
% figure(41); hold on;
% plot([0,rVP], y11*0.9,'-k','LineWidth',1.2)

%% diameter/psf in different angles
rV = [10,20,30,40,50,60,70,80,90,100];%
Vtot0 = -10e-3;
WidthVes0 = 10;%24
nRpt = 5;
snr = 1;
agl = [0,15,30,45,60,65,70,75,90];
%%
CBVs = zeros(length(rV), nRpt);
PDIs = zeros(length(rV), nRpt);
nparticles = zeros(length(rV), nRpt);

% tic
%get IQ0 for calculating noise
[IQ0, FIQ, g1, Numer, tCoor, fCoor, tauCoor, Noise, nparticle] = Numerical_g1fUS_bx(0, WidthVes0, snr);
for j = 1:length(agl)
for k = 1: length(rV)
    disp(num2str(k));
for i = 1:nRpt
    %generate noise with snr level each run
    IQn = awgn(squeeze(IQ0),snr,'measured','db');
    Noise = IQn - squeeze(IQ0);
    %calculation
    [IQ, FIQ, g1, Numer, tCoor, fCoor, tauCoor, Noise, nparticle] = Numerical_g1fUS_bx(agl(j), rV(k), Noise);
    absg1 = abs(squeeze(g1));
    CBVs(k,j,i) = ((absg1(1))./(1-(absg1(1))));
    PDIs(k,j,i) = sum(abs(squeeze(IQ)).^2);
    nparticles(k,j,i) = nparticle;
end
end
end
% toc

%%
% load(['D:\g1_based_fUS\NumericalSimulation\','agl0-15-90wdth10-10-100psf100.mat']);

mCBV = mean(CBVs,3); stdCBV = std(CBVs,1,3);
mnptcl = mean(nparticles,3); stdnptcl = std(nparticles,1,3);
mPDI = mean(PDIs,3); stdPDI = std(PDIs,1,3);

figure; errorbar(mnptcl./mnptcl(1,:), mCBV./mCBV(1,:), stdCBV./mCBV(1,:), 'LineWidth',2); 
hold on; plot(mnptcl./mnptcl(1,:), mnptcl./mnptcl(1,:), '-.k');
xlabel('True rCBV');
ylabel('Measured rCBV');
legend('0^{o}','15^{o}','30^{o}','45^{o}','60^{o}','65^{o}','70^{o}','75^{o}','90^{o}')
figure; errorbar(mnptcl./mnptcl(1,:), mPDI./mPDI(1,:), stdPDI./mPDI(1,:), 'LineWidth',2);
hold on; plot(mnptcl./mnptcl(1,:), mnptcl./mnptcl(1,:), '-.k');
xlabel('True rCBV');
ylabel('Measured rCBV(PDI)');
legend('0^{o}','15^{o}','30^{o}','45^{o}','60^{o}','65^{o}','70^{o}','75^{o}','90^{o}')
%%
rVP = rV/(100*cos(agl(4)/180*pi)); 
[p11, S] = polyfit(rVP(1:end-7)',mCBV(1:end-7,4),1);
y11 = polyval(p11,rVP(1:end));
fig = figure; set(fig,'Position',[800 400 800 300]);
subplot(121);
errorbar(rVP(1:end), mCBV(1:end,4), stdCBV(1:end,4), '.b','MarkerSize',10,'LineWidth',1); hold on; plot(rVP(1:end), y11,'-b','LineWidth',0.75)
hold on; box off; %axis equal tight;
xlabel(['Diameter/PSF(fixed)']);
ylabel(['CBV_{index}']); xlim([0,1])
%h = legend({'','Fit: k = 1.007, b = 0.039'});set(h,'box','off');
title(['k=',num2str(p11(1)),' b=',num2str(p11(2))])
%%
fig = figure; set(fig,'Position',[800 400 800 300]);
COLOR = [0 0.4470 0.7410;
0.8500 0.3250 0.0980;
0.9290 0.6940 0.1250;
0.4940 0.1840 0.5560;
0.4660 0.6740 0.1880;
0.05 0.92 0.9;
0.92 0.05 0.9;
0. 0. 0.9;
0.3010 0.7450 0.9330;
0.456 0.89 0.27];

for j = 1: (length(agl)-4)
rVP = rV/(100*(cos(agl(j)/180*pi)^0.2));
[p11, S] = polyfit(rVP(1:end-7)',mCBV(1:end-7,j),1);
y11 = polyval(p11,rVP(1:end));
subplot(121); hold on;
errorbar(rVP(1:end), mCBV(1:end,j), stdCBV(1:end,j),'MarkerSize',10,'LineWidth',1,'Color',COLOR(j,:)); %hold on; plot(rVP(1:end), y11,'-b','LineWidth',0.75)
hold on; box off; %axis equal tight;
xlabel(['Diameter/(PSF*cos(\theta)^{0.2})']);
ylabel(['Normalized CBV_{index} [a.u.]']); xlim([0,1])
%h = legend({'','Fit: k = 1.007, b = 0.039'});set(h,'box','off');
% title(['k=',num2str(p11(1)),' b=',num2str(p11(2))])
legend('0^{\circ}','15^{\circ}','30^{\circ}','45^{\circ}','60^{\circ}','65^{o}','70^{o}','75^{\circ}')
fitline(j,:) = p11;
ylim([0,6]);
end
% figure; plot(agl,(cos(agl/180*pi))); hold on;  plot(agl,(cos(agl/180*pi)).^0.5); hold on; plot(agl,(cos(agl/180*pi)).^0.2);
% legend('cos(\theta)','cos(\theta)^{0.5}','cos(\theta)^{0.2}'); grid on;

%%
T_rVP = []; T_mCBV = [];%
jend = 6*ones(5,1);
for j = 1:(length(agl)-4)
rVP = rV/(100*(cos(agl(j)/180*pi)^0.2));
T_rVP = [T_rVP,rVP(1:jend(j))];
T_mCBV  = [T_mCBV;mCBV(1:jend(j),j)];
end
[p11, S] = polyfit(T_rVP,T_mCBV',1);
y11 = polyval(p11,[0,rVP]);
figure(5); hold on;
plot([0,rVP], y11,'-.k','LineWidth',1.2)
%%
figure(5); hold on;
linestyle = {':','- -'};
for j = 6: (length(agl)-2)
rVP = rV/(100*(cos(agl(j)/180*pi)^0.2));
[p11, S] = polyfit(rVP(1:end-7)',mCBV(1:end-7,j),1);
y11 = polyval(p11,rVP(1:end));
yyaxis right;
errorbar(rVP(1:end), mCBV(1:end,j), stdCBV(1:end,j),'MarkerSize',10,'LineWidth',1,'Color',0.7*ones(1,3),'LineStyle',linestyle{j-5}); 
hold on; box off; %axis equal tight;
xlabel(['Diameter/(PSF*cos(\theta)^{0.2})']);
xlim([0,1])
%legend('65^{o}','70^{o}','75^{\circ}')
fitline(j,:) = p11;
end
ax = gca;
ax.YAxis(2).Color = 0.7*ones(1,3);