
%%%% load X:\G1based fUS
%%%% data\vUS\RESULT-WhiskerStim_100321_BL4_ForepawUnsynch-vUS %%%%%%% run
%%%% MAIN_vUS_invivo.m

cropX = 30:190;%40:205;
cropY = 10:142-10;%6:142-10;
% CR mask
mCR = double(sum(CR(cropY,cropX,:),3)>0);

addpath D:\CODE\DataProcess\A-US-fUS-g1\SubFunctions
[VzCmap, VzCmapDn,VzCmapUp]=Colormaps_fUS;

%% compare colordoppler(no directional filter)Vd and improved colordoppler(with mask)Vcz
Vd = ColorDoppler(sIQ(cropY,cropX,:), PRSSinfo);

figure; subplot(221);imagesc(Vd);caxis([0,30]);axis image; title('colorDoppler'); colorbar;
subplot(222); imagesc(Vcz(cropY,cropX,2));caxis([0,30]);axis image;title('improved colorDoppler');colorbar;

subplot(223);imagesc(Vd);caxis([-30,0]);axis image;
cmp = colormap; cmp = flipud(cmp); colormap(cmp); title('colorDoppler'); colorbar;
subplot(224); imagesc(Vcz(cropY,cropX,1));caxis([-30,0]);axis image;title('improved colorDoppler');colorbar;

%% compare colordoppler(no directional filter) and improved colordoppler(without mask)Vcz_
[nz0,nx0,nt]=size(sIQ(cropY,cropX,:));
fIQ=sysNoiseRemove(sIQ(cropY,cropX,:),PRSSinfo.rFrame);
%fIQ = fft(sIQ,nt,3);
Vcz0 = [];Vcz_ = [];
for iNP=1:2
    iFIQ=zeros(size(fIQ));
    iFIQ(:,:,(floor(nt/2))*(iNP-1)+1:floor(nt/2)*iNP)=fIQ(:,:,(floor(nt/2))*(iNP-1)+1:floor(nt/2)*iNP);
    iIQ=(ifft(iFIQ,nt,3));
    iVcz=(ColorDoppler(iIQ,PRSSinfo));
    Vcz0(:,:,iNP)=-1*iVcz*1e3;%mm/s
    Vcz_(:,:,iNP)=imresize(Vcz0(:,:,iNP),[nz0,nx0]*PRSSinfo.rfnScale,'nearest');
end

figure; subplot(221);imagesc(Vd);caxis([0,30]);axis image; title('colorDoppler'); colorbar;
cmp = colormap; colormap(cmp);
subplot(222); imagesc(Vcz_(:,:,2));caxis([0,30]);axis image;title('improved colorDoppler'); colorbar;

subplot(223);imagesc(Vd);caxis([-30,0]);axis image;
cmp = colormap; cmp = flipud(cmp); colormap(cmp); title('colorDoppler'); colorbar;
subplot(224); imagesc(Vcz_(:,:,1));caxis([-30,0]);axis image;title('improved colorDoppler'); colorbar;

%% no directional filter: compare colordoppler(freq)(Vd) and 'standard' doppler velocimetry(Vstd)
Vstd = ColorDoppler_std(sIQ(cropY,cropX,:), PRSSinfo);
figure; subplot(221);imagesc(Vd);caxis([0,30]);axis image; title('colorDoppler'); colorbar;
cmp = colormap; colormap(cmp);
subplot(222); imagesc(Vstd);caxis([0,30]);axis image;title('std colorDoppler'); colorbar;
subplot(223);imagesc(Vd);caxis([-30,0]);axis image;
cmp = colormap; cmp = flipud(cmp); colormap(cmp); title('colorDoppler'); colorbar;
subplot(224); imagesc(Vstd);caxis([-30,0]);axis image;title('std colorDoppler'); colorbar;

%% with directional filter: compare colordoppler(freq)(Vd_) and 'standard' doppler velocimetry(Vstd_)
[nz0,nx0,nt]=size(sIQ(cropY,cropX,:));
fIQ=sysNoiseRemove(sIQ(cropY,cropX,:),PRSSinfo.rFrame);
%fIQ = fft(sIQ,nt,3);
Vd0_ = [];Vd = [];Vstd0_=[];Vstd_=[];
for iNP=1:2
    iFIQ=zeros(size(fIQ));
    iFIQ(:,:,(floor(nt/2))*(iNP-1)+1:floor(nt/2)*iNP)=fIQ(:,:,(floor(nt/2))*(iNP-1)+1:floor(nt/2)*iNP);
    iIQ=(ifft(iFIQ,nt,3));
    iVd_= ColorDoppler(iIQ,PRSSinfo);
    Vd0_(:,:,iNP)=-1*iVd_*1e3;%mm/s
    Vd_(:,:,iNP)=imresize(Vd0_(:,:,iNP),[nz0,nx0]*PRSSinfo.rfnScale,'nearest');
    
    iVstd_= ColorDoppler_std(iIQ,PRSSinfo);
    Vstd0_(:,:,iNP)=-1*iVstd_*1e3;%mm/s
    Vstd_(:,:,iNP)=imresize(Vstd0_(:,:,iNP),[nz0,nx0]*PRSSinfo.rfnScale,'nearest');
end

figure; subplot(221);imagesc(Vd_(:,:,2));caxis([0,30]);axis image; title('colorDoppler'); colorbar;
cmp = colormap; colormap(cmp);
subplot(222); imagesc(Vstd_(:,:,2));caxis([0,30]);axis image;title('std colorDoppler'); colorbar;

subplot(223);imagesc(Vd_(:,:,1));caxis([-30,0]);axis image;
cmp = colormap; cmp = flipud(cmp); colormap(cmp); title('colorDoppler'); colorbar;
subplot(224); imagesc(Vstd_(:,:,1));caxis([-30,0]);axis image;title('std colorDoppler'); colorbar;

%%
figure; subplot(221);imagesc(V(cropY,cropX,2));caxis([0,30]);axis image; title('vUS'); colorbar;
cmp = colormap; colormap(cmp);
subplot(222); imagesc(sumGGV0(cropY,cropX,2));caxis([0,10*0.8]);axis image;title('CBF-speed index'); colorbar;
subplot(223);imagesc(V(cropY,cropX,1));caxis([-30,0]);axis image;
cmp = colormap; cmp = flipud(cmp); colormap(cmp); title('vUS'); colorbar;
subplot(224); imagesc(-sumGGV0(cropY,cropX,1));caxis([-10*0.8,0]);axis image;title('CBF-speed index'); colorbar;

%% plots
msk(:,:,1) = abs(V(cropY,cropX,1))>0;
msk(:,:,2) = abs(V(cropY, cropX,2))>0;
msk(:,:,3) = (msk(:,:,1)+msk(:,:,2))>0;

fig = figure;
set(fig, 'Position', [100, 100, 400, 300])
imagesc(-Vstd.*msk(:,:,3));%
caxis([-30,30]*2.5/3);axis image;title('std colorDoppler'); colorbar;
axis off;
colormap(VzCmap)
%%
npmsk = Vstd>0;
fig = figure;
set(fig, 'Position', [100, 100, 400, 300])
imagesc((-sumGG0(cropY,cropX,3).*npmsk+sumGG0(cropY,cropX,3).*(1-npmsk)).*msk(:,:,3));% 
caxis([-2,2]);
axis image;title('CBF-speed index'); colorbar;
colormap(VzCmap)
axis off;
%%
fig = figure;
set(fig, 'Position', [100, 100, 400, 300])
imagesc(-Vd);
caxis([-30,30]);
axis image;title('ColorDoppler'); colorbar;
colormap(VzCmap)
axis off;

%% after directional filter
cRange = [-80, 80];
fig = figure;
set(fig, 'Position', [100, 100, 400, 300])
h1 = axes;
imagesc(Vstd_(:,:,1));%.*CR(cropY,cropX,1)
caxis(cRange*1e3);axis image;title('std colorDoppler'); colorbar;
axis off;
colormap(VzCmap)
hold on
h2 = axes;
imagesc(Vstd_(:,:,2).*CR(cropY,cropX,2));
caxis(cRange*1e3);axis image;title('std colorDoppler'); colorbar;
alpha(h2,double((Vstd_(:,:,2)<7e4)&(CR(cropY,cropX,2)>0)));
axis off;
colormap(VzCmap)
linkaxes([h1, h2])

cRange = [-80, 80];
fig = figure;
set(fig, 'Position', [100, 100, 400, 300])
h1 = axes;
imagesc(Vstd_(:,:,1));%.*CR(cropY,cropX,1)
caxis(cRange*1e3);axis image;title('std colorDoppler'); colorbar;
axis off;
colormap(VzCmap)

fig = figure;
set(fig, 'Position', [100, 100, 400, 300])
h2 = axes;
imagesc(Vstd_(:,:,2));%.*CR(cropY,cropX,2)
caxis(cRange*1e3);axis image;title('std colorDoppler'); colorbar;
%alpha(h2,double((Vstd_(:,:,2)<7e4)&(CR(cropY,cropX,2)>0)));
axis off;
colormap(VzCmap)
linkaxes([h1, h2])

%%
cRange = [-80, 80];
fig = figure;
set(fig, 'Position', [100, 100, 400, 300])
h1 = axes;
imagesc(Vd_(:,:,1).*CR(cropY,cropX,1));
caxis(cRange*1e3);axis image;title('ColorDoppler'); colorbar;
axis off;
colormap(VzCmap)
hold on
h2 = axes;
imagesc(Vd_(:,:,2).*CR(cropY,cropX,2));
caxis(cRange*1e3);axis image;title('ColorDoppler'); colorbar;
alpha(h2,double((Vd_(:,:,2)<7e4)&(CR(cropY,cropX,2)>0)));
axis off;
colormap(VzCmap)
linkaxes([h1, h2])

%%
fig = figure;
set(fig, 'Position', [100, 100, 400, 300])
h1 = axes;
imagesc(Vcz(cropY,cropX,1));
caxis([-30,30]*1e0);axis image;title('improved ColorDoppler'); colorbar;
axis off;
colormap(VzCmap)
hold on
h2 = axes;
imagesc(Vcz(cropY,cropX,2));
caxis([-30,30]*1e0);axis image;title('improved ColorDoppler'); colorbar;
alpha(h2,double(Vcz(cropY,cropX,2)>0.5));
axis off;
linkaxes([h1, h2])
colormap(VzCmap)

% fig = figure;
% set(fig, 'Position', [100, 100, 400, 300])
% h1 = axes;
% imagesc(-sumGG0(cropY,cropX,1));
% caxis([-5,5]*1e0);axis image;title('CBF-speed index'); colorbar;
% axis off;
% colormap(VzCmap)
% hold on
% h2 = axes;
% imagesc(sumGG0(cropY,cropX,2));
% caxis([-5,5]*1e0);axis image;title('CBF-speed index'); colorbar;
% alpha(h2,double(Vcz(cropY,cropX,2)>0.5));
% axis off;
% linkaxes([h1, h2])
% colormap(VzCmap)

%% CBF-speed index %%% overlap
%VzCmap_ = [flipud(VzCmap(1:32,:)); flipud(VzCmap(32:end,:))];
msk(:,:,1) = abs(V(cropY,cropX,1))>0;
msk(:,:,2) = abs(V(cropY, cropX,2))>0;
msk(:,:,3) = (msk(:,:,1)+msk(:,:,2))>0;

cRange = [-1,1]*1.5;
fig = figure;
set(fig,'Position',[100, 100, 400, 300])
h1 = axes;
imagesc(-sumGG0(cropY, cropX,1).*msk(:,:,1)); % veins 
colorbar;axis equal tight;colormap(h1,(VzCmap));
%alpha(h1,double(-sumGGV0(:,:,1)<-0.5));
caxis(cRange)
axis(h1,'off')
h2 = axes;
imagesc(sumGG0(cropY, cropX,2).*msk(:,:,2)); % arteries
colorbar;colormap(h2,(VzCmap));axis equal tight;
alpha(h2,(double(sumGG0(cropY, cropX,2)>0.3).*(1-(msk(:,:,1)-msk(:,:,2))>0)));
caxis(cRange)
axis(h2,'off')
linkaxes([h1,h2]);

cRange = [-1,1]*1.5;
fig = figure;
set(fig,'Position',[100, 100, 400, 300])
h1 = axes;
imagesc(sumGG0(cropY, cropX,1).*msk(:,:,1)); % veins 
colorbar;axis equal tight;colormap(h1,(flipud(VzCmapUp)));
%alpha(h1,double(-sumGGV0(:,:,1)<-0.5));
caxis(cRange)
axis(h1,'off')

fig = figure;
set(fig,'Position',[100, 100, 400, 300])
h2 = axes;
imagesc(sumGG0(cropY, cropX,2).*msk(:,:,2)); % arteries 
colorbar;colormap(h2,(VzCmapDn));axis equal tight;
%alpha(h2,(double(msumGG(cropY, cropX,2)>0.3).*(1-(msk(:,:,1)-msk(:,:,2))>0)));
caxis(cRange)
axis(h2,'off')

%% plot spectrum & GG
pixX = 77+9; pixY = 144+29;%(25, 140)77,144; 30,148
pixX = 31; pixY = 146;
fCoor = linspace(-2500,2500,1000);
tau = linspace(0,20,100);
pix = squeeze(sIQ(pixX,pixY,:));
fpix = abs(fftshift(fft(pix)));
fig = figure;
set(fig, 'Position', [100, 100, 300, 200])
plot(fCoor, fpix,'k');
xlim([-1000,1000]*1);
xlabel('Frequency [Hz]');
ylabel('[a.u.]')
ylim([0,max(fpix)])

GGpix = sIQ2GG(sIQ(pixX,pixY,:),PRSSinfo);
fig = figure;
set(fig, 'Position', [100, 100, 300, 200])
plot(tau,squeeze(abs(GGpix)),'k');
xlabel('\tau [ms]');
ylabel('abs(g_1(\tau))');

fIQpix = sysNoiseRemove(sIQ(pixX,pixY,:),PRSSinfo.rFrame);
for iNP=1:2
    iFIQ=zeros(size(fIQpix));
    iFIQ(:,:,(floor(nt/2))*(iNP-1)+1:floor(nt/2)*iNP)=fIQpix(:,:,(floor(nt/2))*(iNP-1)+1:floor(nt/2)*iNP);
    iIQ=(ifft(iFIQ,nt,3));
    iGGpix(:,:,:,iNP)=sIQ2GG(iIQ,PRSSinfo);
end

fig = figure;
set(fig, 'Position', [100, 100, 300, 200])
plot(fCoor(end/2+1:end), fpix(end/2+1:end))
hold on;
plot(fCoor(1:end/2), fpix(1:end/2))
xlabel('Frequency [Hz]');
ylabel('[a.u.]')
ylim([0,max(fpix)])
xlim([-1000,1000]*1)

fig = figure;
set(fig, 'Position', [100, 100, 300, 200])
plot(tau,movmean(squeeze(abs(iGGpix(:,:,:,1))),2));
hold on
plot(tau,movmean(squeeze(abs(iGGpix(:,:,:,2))),2));
xlabel('\tau [ms]');
ylabel('abs(g_1(\tau))');

%% compare GGpix 1,2,3 at different SNR
fig = figure;
set(fig, 'Position', [100, 100, 300, 200])
plot(tau,movmean(GGpixs1,2));
hold on
plot(tau,movmean(GGpixs2,4));
hold on
plot(tau,movmean(GGpixs4,8));
hold on
plot(tau,movmean(GGpixs3,1));
xlabel('\tau [ms]');
ylabel('abs(g_1(\tau))');


%% compare improved colordoppler and tl-fUS CBF-speed index (sumGGV)

function Vstd = ColorDoppler_std(sIQ, PRSSinfo)
[nz, nx, nxRpt] = size(sIQ) ;
if PRSSinfo.g1nT>nxRpt-PRSSinfo.g1nTau
    PRSSinfo.g1nT=nxRpt-PRSSinfo.g1nTau-PRSSinfo.g1StartT+1;
%     disp(['Warning: nt is larger than nxRpt-ntau, and is modified to be nxRpt-ntau=',num2str(PRSSinfo.g1nT),'!']);
end
int = 1;
Numer = mean((conj(sIQ(:,:,PRSSinfo.g1StartT:PRSSinfo.g1StartT-1+PRSSinfo.g1nT))).*(sIQ(:,:,int+PRSSinfo.g1StartT:int+PRSSinfo.g1StartT-1+PRSSinfo.g1nT)),3);
phase = angle(Numer);
Vstd = phase*(PRSSinfo.C/PRSSinfo.f0)/(4*pi*(int/PRSSinfo.rFrame))*1e3;%mm/s
end

