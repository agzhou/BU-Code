clear all;
%% figure1 simulation:
addpath D:\g1_based_fUS\NumericalSimulation\
addpath D:\g1_based_fUS\A-US-fUS-g1\SubFunctions\

% snr = 20;
% Vtot = -10e-3;
% nScatter = 10;

Vtot = [-5e-3, -10e-3];
WidthVes = [10,12];
%%
Vtot = [-10e-3, -15e-3; 10e-3, 15e-3];
WidthVes = [12,24;12,24];

%% CBFv
[IQ, FIQ, g1, Numer, tCoor, fCoor, tauCoor] = Numerical_g1fUS_bx(Vtot(1), WidthVes(1),30,0);
IQdata(:,1) = squeeze(IQ);
FIQdata(:,1) = squeeze(FIQ);
absg1(:,1) = abs(squeeze(g1));
rg1(:,1) = real(squeeze(g1));
ig1(:,1) = imag(squeeze(g1));
clear IQ FIQ g1
[IQ, FIQ, g1, Numer, tCoor, fCoor, tauCoor] = Numerical_g1fUS_bx(Vtot(2), WidthVes(1),30,0);
IQdata(:,2) = squeeze(IQ);
FIQdata(:,2) = squeeze(FIQ);
absg1(:,2) = abs(squeeze(g1));
rg1(:,2) = real(squeeze(g1));
ig1(:,2) = imag(squeeze(g1));
clear IQ FIQ g1

fig = figure; 
set(fig,'Position',[500,500,300,300])
subplot(2,1,1), h1 = plot(tCoor, abs(IQdata(:,1)),'b'); h1.Color(4) = 0.3;
hold on; h2 = plot(tCoor, abs(IQdata(:,2)),'r'); h2.Color(4) = 0.3;
xlim([0 600])
set(gca, 'FontSize', 8);
xlabel('Time [ms]')
ylabel('sIQ')
subplot(2,1,2),h1 = plot(fCoor,abs(FIQdata(:,1)),'b'); h1.Color(4) = 0.3;
hold on; h2 = plot(fCoor, abs(FIQdata(:,2)),'r'); h2.Color(4) = 0.3;
xlim([-600 600])
xlabel('Frequency [Hz]')
ylabel('Power');
legend({['V = ',num2str(Vtot(1)*1e3),'mm/s'],['V = ',num2str(Vtot(2)*1e3),'mm/s']});
legend('boxoff');
set(gca, 'FontSize', 8);

fig = figure;
set(fig,'Position',[300,300,250,600])
subplot(311);plot(tauCoor,absg1(:,1),'-b','LineWidth',1);xlabel('\tau [ms]');ylabel('abs(g_{1})');
hold on;plot(tauCoor,absg1(:,2),'-r','LineWidth',1);
%  legend({'v_{x}=5mm/s, v_{z}=5mm/s','v_{x}=10mm/s, v_{z}=5mm/s'},'Location','bestoutside','Orientation','horizontal');
legend({['V = ',num2str(Vtot(1)*1e3),'mm/s'],['V = ',num2str(Vtot(2)*1e3),'mm/s']});
legend('boxoff');
hold on;subplot(312);plot(tauCoor,rg1(:,1),'-b','LineWidth',1);xlabel('\tau [ms]');ylabel('Re(g_{1})');
hold on;plot(tauCoor,rg1(:,2),'-r','LineWidth',1);
hold on;subplot(313);plot(tauCoor,ig1(:,1),'-b','LineWidth',1);xlabel('\tau [ms]');ylabel('Im(g_{1})');
hold on;plot(tauCoor,ig1(:,2),'-r','LineWidth',1);
hold off

%% CBFv 2 vessels
[IQ, FIQ, g1, Numer, tCoor, fCoor, tauCoor] = Numerical2vessels_g1fUS(Vtot(:,1), WidthVes(:,1));
IQdata(:,1) = squeeze(IQ);
FIQdata(:,1) = squeeze(FIQ);
absg1(:,1) = abs(squeeze(g1));
rg1(:,1) = real(squeeze(g1));
ig1(:,1) = imag(squeeze(g1));
clear IQ FIQ g1

[IQ, FIQ, g1, Numer, tCoor, fCoor, tauCoor] = Numerical2vessels_g1fUS(Vtot(:,2), WidthVes(:,1));
IQdata(:,2) = squeeze(IQ);
FIQdata(:,2) = squeeze(FIQ);
absg1(:,2) = abs(squeeze(g1));
rg1(:,2) = real(squeeze(g1));
ig1(:,2) = imag(squeeze(g1));
clear IQ FIQ g1

fig = figure; 
set(fig,'Position',[500,500,300,300])
subplot(2,1,1), h1 = plot(tCoor, abs(IQdata(:,1)),'b'); h1.Color(4) = 0.3;
hold on; h2 = plot(tCoor, abs(IQdata(:,2)),'r'); h2.Color(4) = 0.3;
xlim([0 600])
set(gca, 'FontSize', 8);
xlabel('Time [ms]')
ylabel('sIQ')
subplot(2,1,2),h1 = plot(fCoor,abs(FIQdata(:,1)),'b'); h1.Color(4) = 0.3;
hold on; h2 = plot(fCoor, abs(FIQdata(:,2)),'r'); h2.Color(4) = 0.3;
xlim([-600 600])
xlabel('Frequency [Hz]')
ylabel('Power');
legend({['V = ',num2str(Vtot(1,1)*1e3),'mm/s'],['V = ',num2str(Vtot(1,2)*1e3),'mm/s']});
legend('boxoff');
set(gca, 'FontSize', 8);

fig = figure;
set(fig,'Position',[300,300,250,600])
subplot(311);plot(tauCoor,absg1(:,1),'-b','LineWidth',1);xlabel('\tau [ms]');ylabel('abs(g_{1})');
hold on;plot(tauCoor,absg1(:,2),'-r','LineWidth',1);
%  legend({'v_{x}=5mm/s, v_{z}=5mm/s','v_{x}=10mm/s, v_{z}=5mm/s'},'Location','bestoutside','Orientation','horizontal');
legend({['V = ',num2str(Vtot(1,1)*1e3),'mm/s'],['V = ',num2str(Vtot(1,2)*1e3),'mm/s']});
legend('boxoff');
hold on;subplot(312);plot(tauCoor,rg1(:,1),'-b','LineWidth',1);xlabel('\tau [ms]');ylabel('Re(g_{1})');
hold on;plot(tauCoor,rg1(:,2),'-r','LineWidth',1);
hold on;subplot(313);plot(tauCoor,ig1(:,1),'-b','LineWidth',1);xlabel('\tau [ms]');ylabel('Im(g_{1})');
hold on;plot(tauCoor,ig1(:,2),'-r','LineWidth',1);
hold off

%% CBV
[IQ, FIQ, g1,Numer, tCoor, fCoor, tauCoor, Noise, nparticle] = Numerical_g1fUS_bx(Vtot(2),WidthVes(1),1);
IQdata(:,1) = squeeze(IQ);
FIQdata(:,1) = squeeze(FIQ);
absg1(:,1) = abs(squeeze(g1));
rg1(:,1) = real(squeeze(g1));
ig1(:,1) = imag(squeeze(g1));
G1(:,1) = squeeze(Numer);
disp(['G1(1):', num2str(squeeze(g1(1,1,1)))])
disp(['number of particles: ', num2str(nparticle)])
clear IQ FIQ g1

[IQ, FIQ, g1, Numer, tCoor, fCoor, tauCoor, Noise, nparticle] = Numerical_g1fUS_bx(Vtot(2),WidthVes(2), Noise);%WidthVes(2)
IQdata(:,2) = squeeze(IQ);
FIQdata(:,2) = squeeze(FIQ);
absg1(:,2) = abs(squeeze(g1));
rg1(:,2) = real(squeeze(g1));
ig1(:,2) = imag(squeeze(g1));
G1(:,2) = squeeze(Numer);
disp(['G1(1):', num2str(squeeze(g1(1,1,1)))])
disp(['number of particles: ', num2str(nparticle)])
clear IQ FIQ g1

fig = figure; 
set(fig,'Position',[500,500,300,300])
subplot(2,1,1), h1 = plot(tCoor, abs(IQdata(:,1)),'b'); h1.Color(4) = 0.3;
hold on; h2 = plot(tCoor, abs(IQdata(:,2)),'r'); h2.Color(4) = 0.3;
xlim([0 600])
set(gca, 'FontSize', 8);
xlabel('Time [ms]')
ylabel('sIQ')
subplot(2,1,2),h1 = plot(fCoor,abs(FIQdata(:,1)),'b'); h1.Color(4) = 0.3;
hold on; h2 = plot(fCoor, abs(FIQdata(:,2)),'r'); h2.Color(4) = 0.3;
xlim([-600 600])
xlabel('Frequency [Hz]')
ylabel('Power');
legend({['WidthVes = ',num2str(WidthVes(1,1))],['WidthVes = ',num2str(WidthVes(1,2))]});
legend('boxoff');
set(gca, 'FontSize', 8);

fig = figure;
set(fig,'Position',[300,300,250,600])
subplot(311);plot(tauCoor,absg1(:,1),'-b','LineWidth',1);xlabel('\tau [ms]');ylabel('abs(g_{1})');
hold on;plot(tauCoor,absg1(:,2),'-r','LineWidth',1);
%  legend({'v_{x}=5mm/s, v_{z}=5mm/s','v_{x}=10mm/s, v_{z}=5mm/s'},'Location','bestoutside','Orientation','horizontal');
legend({['WidthVes = ',num2str(WidthVes(1,1))],['WidthVes = ',num2str(WidthVes(1,2))]});
legend('boxoff');
hold on;subplot(312);plot(tauCoor,rg1(:,1),'-b','LineWidth',1);xlabel('\tau [ms]');ylabel('Re(g_{1})');
hold on;plot(tauCoor,rg1(:,2),'-r','LineWidth',1);
hold on;subplot(313);plot(tauCoor,ig1(:,1),'-b','LineWidth',1);xlabel('\tau [ms]');ylabel('Im(g_{1})');
hold on;plot(tauCoor,ig1(:,2),'-r','LineWidth',1);
hold off
%% measured 
CBFv_ind = sqrt(abs(log(absg1(1,:)./absg1(21,:))));
rCBFv = CBFv_ind(2)/CBFv_ind(1);
disp(['measured rCBFv: ', num2str(rCBFv*100), '%'])

CBV_ind = (absg1(1,:))./(1-(absg1(1,:)));
rCBV = CBV_ind(2)/CBV_ind(1);
disp(['measured rCBV: ', num2str(rCBV*100), '%'])

rPDI = sum(abs(IQdata(:,2)).^2)/sum(abs(IQdata(:,1)).^2);
disp(['measured rPDI: ', num2str(rPDI*100), '%'])
%% CBV 2 vessels
[IQ, FIQ, g1, Numer, tCoor, fCoor, tauCoor, Noise] = Numerical2vessels_g1fUS(Vtot(:,1), WidthVes(:,1));
IQdata(:,1) = squeeze(IQ);
FIQdata(:,1) = squeeze(FIQ);
absg1(:,1) = abs(squeeze(g1));
rg1(:,1) = real(squeeze(g1));
ig1(:,1) = imag(squeeze(g1));
clear IQ FIQ g1
[IQ, FIQ, g1, Numer, tCoor, fCoor, tauCoor] = Numerical2vessels_g1fUS(Vtot(:,1), WidthVes(:,2), Noise);
IQdata(:,2) = squeeze(IQ);
FIQdata(:,2) = squeeze(FIQ);
absg1(:,2) = abs(squeeze(g1));
rg1(:,2) = real(squeeze(g1));
ig1(:,2) = imag(squeeze(g1));
clear IQ FIQ g1

fig = figure; 
set(fig,'Position',[500,500,300,300])
subplot(2,1,1), h1 = plot(tCoor, abs(IQdata(:,1)),'b'); h1.Color(4) = 0.3;
hold on; h2 = plot(tCoor, abs(IQdata(:,2)),'r'); h2.Color(4) = 0.3;
xlim([0 600])
set(gca, 'FontSize', 8);
xlabel('Time [ms]')
ylabel('sIQ')
subplot(2,1,2),h1 = plot(fCoor,abs(FIQdata(:,1)),'b'); h1.Color(4) = 0.3;
hold on; h2 = plot(fCoor, abs(FIQdata(:,2)),'r'); h2.Color(4) = 0.3;
xlim([-600 600])
xlabel('Frequency [Hz]')
ylabel('Power');
legend({['WidthVes = ',num2str(WidthVes(1,1))],['WidthVes = ',num2str(WidthVes(1,2))]});
legend('boxoff');
set(gca, 'FontSize', 8);

fig = figure;
set(fig,'Position',[300,300,250,600])
subplot(311);plot(tauCoor,absg1(:,1),'-b','LineWidth',1);xlabel('\tau [ms]');ylabel('abs(g_{1})');
hold on;plot(tauCoor,absg1(:,2),'-r','LineWidth',1);
%  legend({'v_{x}=5mm/s, v_{z}=5mm/s','v_{x}=10mm/s, v_{z}=5mm/s'},'Location','bestoutside','Orientation','horizontal');
legend({['WidthVes = ',num2str(WidthVes(1,1))],['WidthVes = ',num2str(WidthVes(1,2))]});
legend('boxoff');
hold on;subplot(312);plot(tauCoor,rg1(:,1),'-b','LineWidth',1);xlabel('\tau [ms]');ylabel('Re(g_{1})');
hold on;plot(tauCoor,rg1(:,2),'-r','LineWidth',1);
hold on;subplot(313);plot(tauCoor,ig1(:,1),'-b','LineWidth',1);xlabel('\tau [ms]');ylabel('Im(g_{1})');
hold on;plot(tauCoor,ig1(:,2),'-r','LineWidth',1);
hold off



%% CBF vs. noise
nRpt = 20;
noise = [3, 1, -1, -3, -5, -7];%dB
%%
tic
for k = 1: length(noise)
for i = 1:nRpt
%     [IQ, FIQ, g1, Numer, tCoor, fCoor, tauCoor] = Numerical_g1fUS(Vtot(1), WidthVes(1), noise(k));
%     g1s(:,i,1) = g1;
%     absg1 = abs(squeeze(g1));
%     CBFv1(i) = sqrt(log(absg1(1)./absg1(10)));
%     [IQ, FIQ, g1, Numer, tCoor, fCoor, tauCoor] = Numerical_g1fUS(Vtot(2), WidthVes(1), noise(k));
%     g1s(:,i,2) = g1;
%     absg1 = abs(squeeze(g1));
%     CBFv2(i) = sqrt(log(absg1(1)./absg1(10)));
%     rCBF(i) = CBFv2(i)./CBFv1(i);
    
    [IQ, FIQ, g1, Numer, tCoor, fCoor, tauCoor, Noise] = Numerical_g1fUS(Vtot(1), WidthVes(1), noise(k));
    g1s(:,i,3) = g1;
    absg1 = abs(squeeze(g1));
    CBV1(i) = absg1(1)./(1-absg1(1));
    [IQ, FIQ, g1, Numer, tCoor, fCoor, tauCoor] = Numerical_g1fUS(Vtot(1), WidthVes(2), Noise);
    g1s(:,i,4) = g1;
    absg1 = abs(squeeze(g1));
    CBV2(i) = absg1(1)./(1-absg1(1));
    rCBV(i) = CBV2(i)./CBV1(i);
end
% mrCBF(:,k) = mean(rCBF);
% stdCBF(:,k) = std(rCBF);

mrCBV(:,k) = mean(rCBV);
stdCBV(:,k) = std(rCBV);
end
toc
% mrCBF1 = mrCBF;
% stdCBF1 = stdCBF;

mrCBV1 = mrCBV;
stdCBV1 = stdCBV;
%%

figure; errorbar(noise, mrCBF1, stdCBF1,'b', 'LineWidth',2);
hold on; plot(noise, Vtot(2)/Vtot(1)*ones(size(noise)),'-.k');
xlim([-8,4]); xlabel('Noise level [dB]');
ylim([1,2]); ylabel('Measured rCBFv');
set(gca,'XDir','reverse');

figure; errorbar(noise, mrCBV1, stdCBV1,'b', 'LineWidth',2);
hold on; plot(noise, WidthVes(2)/WidthVes(1)*ones(size(noise)),'-.k');
xlim([-8,4]); xlabel('Noise level [dB]');
ylim([1,2.5]); ylabel('Measured rCBV');
set(gca,'XDir','reverse');

%% CBF vs. v0
rV = 1.2;
Vtot0 = [21, 23, 25, 27, 29]*(-1e-3)-20e-3;
%%
tic
for k = 1: length(Vtot0)
for i = 1:nRpt
    [IQ, FIQ, g1, Numer, tCoor, fCoor, tauCoor] = Numerical_g1fUS(Vtot0(k), WidthVes(1), -1);
    g1s(:,i,1) = g1;
    absg1 = abs(squeeze(g1));
    CBFv1(i) = sqrt(log(absg1(1)./absg1(10)));
    [IQ, FIQ, g1, Numer, tCoor, fCoor, tauCoor] = Numerical_g1fUS(Vtot0(k)*rV, WidthVes(1), -1);
    g1s(:,i,2) = g1;
    absg1 = abs(squeeze(g1));
    CBFv2(i) = sqrt(log(absg1(1)./absg1(10)));
    rCBF(i) = CBFv2(i)./CBFv1(i);
end
mrCBF(:,k) = mean(rCBF);
stdCBF(:,k) = std(rCBF);
end
toc
mrCBF2 = mrCBF(1:end-1);
stdCBF2 = stdCBF(1:end-1);
%%

figure; errorbar(-Vtot0*1e3, mrCBF2, stdCBF2,'b', 'LineWidth',2);
hold on; plot(-Vtot0*1e3, rV*ones(size(Vtot0)), '-.k');
xlim([40,50]); xlabel('CBFv baseline [mm/s]');
ylim([0.9,1.5]); ylabel('Measured rCBFv');
% set(gca,'XDir','reverse');

%% CBF vs.time lag
clear rV Vtot0
rV = 1.2;
Vtot0 = -25*1.2e-3;
taulag = [5,6,7,8,9,10];
%%
tic
for k = 1: length(taulag)
for i = 1:nRpt
    [IQ, FIQ, g1, Numer, tCoor, fCoor, tauCoor] = Numerical_g1fUS(Vtot0, WidthVes(1), -1);
    g1s(:,i,1) = g1;
    absg1 = abs(squeeze(g1));
    CBFv1(i) = sqrt(log(absg1(1)./absg1(taulag(k))));
    [IQ, FIQ, g1, Numer, tCoor, fCoor, tauCoor] = Numerical_g1fUS(Vtot0*rV, WidthVes(1), -1);
    g1s(:,i,2) = g1;
    absg1 = abs(squeeze(g1));
    CBFv2(i) = sqrt(log(absg1(1)./absg1(taulag(k))));
    rCBF(i) = CBFv2(i)./CBFv1(i);
end
mrCBF(:,k) = mean(rCBF);
stdCBF(:,k) = std(rCBF);
end
toc
mrCBF4 = mrCBF;
stdCBF4 = stdCBF;
%%
figure; errorbar(taulag, mrCBF4, stdCBF4,'b', 'LineWidth',2);
hold on; plot(taulag,rV*ones(size(taulag)), '-.k');
xlim([4,11]); xlabel('Selected time lag');
ylim([0.8,1.4]); ylabel('Measured rCBFv');
% set(gca,'XDir','reverse');

%% CBF vs.ratio
clear rV Vtot0
rV = [1.05, 1.1, 1.2, 1.3, 1.4, 1.5];
Vtot0 = -10e-3;
%%
tic
for k = 1: length(rV)
for i = 1:nRpt
    [IQ, FIQ, g1, Numer, tCoor, fCoor, tauCoor] = Numerical_g1fUS(Vtot0, WidthVes(1), 1);
    g1s(:,i,1) = g1;
    absg1 = abs(squeeze(g1));
    CBFv1(i) = sqrt(log(absg1(1)./absg1(10)));
    [IQ, FIQ, g1, Numer, tCoor, fCoor, tauCoor] = Numerical_g1fUS(Vtot0*rV(k), WidthVes(1), 1);
    g1s(:,i,2) = g1;
    absg1 = abs(squeeze(g1));
    CBFv2(i) = sqrt(log(absg1(1)./absg1(10)));
    rCBF(i) = CBFv2(i)./CBFv1(i);
end
mrCBF(:,k) = mean(rCBF);
stdCBF(:,k) = std(rCBF);
end
toc
mrCBF3 = mrCBF;
stdCBF3 = stdCBF;
%%
figure; errorbar(rV, mrCBF3, stdCBF3,'b', 'LineWidth',2);
hold on; plot(rV, rV, '-.k');
xlim([0.9,1.7]); xlabel('True rCBFv');
ylim([0.9,1.7]); ylabel('Measured rCBFv');
% set(gca,'XDir','reverse');
%%
[p11, S] = polyfit(rV,mrCBF,1);
y11 = polyval(p11,rV);
fig = figure; set(fig,'Position',[800 400 800 300]);
subplot(121);
errorbar(rV, mrCBF, stdCBF, '.b','MarkerSize',10,'LineWidth',1); hold on; plot(rV, y11,'-b','LineWidth',0.75)
hold on; box off; axis equal tight;
xlim([0.9,1.7]);xlabel(['Relative V_{set}']);
ylim([0.9,1.7]);ylabel(['Relative CBFv_{index}']); 
%h = legend({'','Fit: k = 1.007, b = 0.039'});set(h,'box','off');
title(['k=',num2str(p11(1)),' b=',num2str(p11(2))])


%% CBV vs.ratio
clear rV Vtot0
rV = [10,11,12,13,14,15]*1;%[10,11,12,13,14,15]*15;%[4,6,8,10,12];[27,54,81,108];*20
Vtot0 = -10e-3;
WidthVes0 = 10;%24
nRpt = 10;
%%
CBVs = zeros(length(rV), nRpt, 2);
PDIs = zeros(length(rV), nRpt, 2);
nparticles = zeros(length(rV), nRpt, 2);
tic
for k = 1: length(rV)
    disp(num2str(k));
for i = 1:nRpt
    [IQ, FIQ, g1, Numer, tCoor, fCoor, tauCoor, Noise, nparticle] = Numerical_g1fUS_bx(Vtot0, WidthVes0, 1);
    absg1 = abs(squeeze(g1));
%     CBV1(i) = (sqrt(absg1(1)-1)./(1-sqrt(absg1(1)-1)));
    CBVs(k,i,1) = ((absg1(1))./(1-(absg1(1))));
    PDIs(k,i,1) = sum(abs(squeeze(IQ)).^2);
    nparticles(k,i,1) = nparticle;
    [IQ, FIQ, g1, Numer, tCoor, fCoor, tauCoor, Noise, nparticle] = Numerical_g1fUS_bx(Vtot0, rV(k), Noise);
%     g1s(:,i,k) = g1;
    absg1 = abs(squeeze(g1));
%     CBV2(i) = (sqrt(absg1(1)-1)./(1-sqrt(absg1(1)-1)));
    CBVs(k,i,2) = ((absg1(1))./(1-(absg1(1))));
    PDIs(k,i,2) = sum(abs(squeeze(IQ)).^2);
    nparticles(k,i,2) = nparticle;
end
end
toc
rCBV = CBVs(:,:,2)./CBVs(:,:,1);
rPDI = PDIs(:,:,2)./PDIs(:,:,1);
rnptcl = nparticles(:,:,2)./nparticles(:,:,1);
mrCBV = mean(rCBV,2); stdCBV = std(rCBV,1,2);
mrnptcl = mean(rnptcl,2); stdrnptcl = std(rnptcl,1,2);
%%
save(['D:\g1_based_fUS\NumericalSimulation\','agl0wdth50.mat'],'CBVs','PDIs','nparticles');
%%
% rV0 = rV/rV(1);
figure; errorbar(mrnptcl, mrCBV/mrCBV(1), stdCBV,'b', 'LineWidth',2); %hold on;errorbar(rV0, mean(CBV0./CBV0(1,:),2), std(CBV0./CBV0(1,:),1,2),'r', 'LineWidth',2);
hold on; plot(mrnptcl, mrnptcl, '-.k');
xlim([0.8,1.8]); xlabel('True rCBV');
ylim([0.8,2]); ylabel('Measured rCBV');
figure; errorbar(mrnptcl, median(rPDI,2), std(rPDI,1,2),'b', 'LineWidth',2);
hold on; plot(mrnptcl, mrnptcl, '-.k');
xlim([0.8,1.8]); xlabel('True rCBV');
ylim([0.8,1.8]); ylabel('Measured rCBV(PDI)');
%%
[p11, S] = polyfit(mrnptcl(1:end),mrCBV(1:end),1);
y11 = polyval(p11,mrnptcl(1:end));
fig = figure; set(fig,'Position',[800 400 800 300]);
subplot(121);
errorbar(mrnptcl(1:end), mrCBV(1:end), stdCBV(1:end), '.b','MarkerSize',10,'LineWidth',1); hold on; plot(mrnptcl(1:end), y11,'-b','LineWidth',0.75)
hold on; box off; axis equal tight;
xlim([0.7,2.3]);xlabel(['Relative Concentration_{set}']);
ylim([0.7,2.3]);ylabel(['Relative CBV_{index}']); 
%h = legend({'','Fit: k = 1.007, b = 0.039'});set(h,'box','off');
title(['k=',num2str(p11(1)),' b=',num2str(p11(2))])

