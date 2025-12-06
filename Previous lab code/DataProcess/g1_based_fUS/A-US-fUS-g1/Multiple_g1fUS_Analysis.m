
clear all; 
addpath D:\g1_based_fUS\A-US-fUS-g1\SubFunctions\

%% I.load GG data
[FileName,FilePath]=uigetfile('V:\G1based fUS data');
fileInfo=strsplit(FileName(1:end-4),'-');
myFile=matfile([FilePath,FileName]);
P=myFile.P;
prompt={'Start Repeat', 'Number of Repeats'};
name='File info'; 
defaultvalue={'1','480','1'};
numinput=inputdlg(prompt,name, 1, defaultvalue);
startRpt=str2num(numinput{1});
nRpt=str2num(numinput{2});          % number of repeat for each coronal plane

indSkipped=1; k = 1; 
for iRpt=startRpt:startRpt+nRpt-1
 iFileInfo=fileInfo;
 iFileInfo{8}=num2str(iRpt);
 iFileName=[strjoin(iFileInfo,'-'),'.mat'];
            %     load([FilePath,iFileName]);
 if exist([FilePath,iFileName],'file')
  myFile=matfile([FilePath,iFileName]);
   GG = myFile.GG0;
%     VV = myFile.V;
%   Numer = myFile.Numer;
%   Denom = myFile.Denom;
%   IQDenom = myFile.IQDenom;
%   GIQ = myFile.GIQ;
%   GGn0 = myFile.GGnp(:,:,:,1);
%   GGp0 = myFile.GGnp(:,:,:,2);
  PDIHP = myFile.PDIHHP;
  eqNoise = myFile.eqNoise;
  disp([iFileName,' was loaded!'])
  else
   disp([iFileName, ' skipped!'])
  SkipFile(indSkipped)=iRpt;
   indSkipped=indSkipped+1;  
 end
 GG0(:,:,:,k) = GG;
% VV0(:,:,:,k) = VV;
%  Numer0(:,:,:,k) = Numer;
%  Denom0(:,:,:,k) = Denom;
%  IQDenom0(:,:,:,k) = IQDenom;
%  GIQ0(:,:,:,k) = GIQ;
%  GGn(:,:,:,k) = GGn0;
%  GGp(:,:,:,k) = GGp0;
 PDIHP0(:,:,:,k) = PDIHP;
 eqNoise0(:,:,k) = eqNoise;
 k = k+1;
end
eqNoise = mean(eqNoise0,3);
PDI = PDIHP0./eqNoise.^1.8;
ImgBG = mean(squeeze(PDI(:,:,3,:)),3).^0.35;
[nz,nx,ntau,~] = size(GG0);

% trial defination 
trial.nBase = 30;%25;
trial.n = 10;
trial.nStim = 15;%5;
trial.nRecover = 30;%25;
trial.nRest = 5; % included in nRecover
trial.nlength = trial.nStim + trial.nRecover;
% trial.parttern = Stim(trial.nBase-10+1:trial.nBase-10+(trial.nStim+trial.nRecover)); 

stim = zeros(nRpt,1);
for k = 1: trial.n % # of trials
stim(trial.nBase+(k-1)*trial.nlength+1: trial.nBase+trial.nStim+(k-1)*trial.nlength) = 1;
end
trial.stim = stim;
xCoor = (1:ntau)*1e3/P.CCFR; % tau axis[ms]

% Global HRF 
% hrf = hemodynamicResponse(1,[2 16 0.5 1 20 0]);
hrf = hemodynamicResponse(1,[1.5 10 0.5 1 20 0 16]);
Stim = zeros(trial.nlength,1);
Stim(5+1:5+trial.nStim,:)=1;
stimhrf = filter(hrf,1,Stim);
figure; plot(stimhrf);

%% II.single pixel analysis
% %manually select datapoint to check gg
% prompt={'New Select(0: No; 1: Yes)'};
% name='Select';
% defaultvalue={'1'};
% numinput=inputdlg(prompt,name, 1, defaultvalue);
% newSlect=str2num(numinput{1});
% 
% if newSlect==1
%     figure(2); imagesc(ImgBG); axis image;colormap(gray);
%     [slt(2), slt(1)]=ginput(1); % [x z]
%     slt=round(slt);
% end
% 
% gg = squeeze(GG0(slt(1), slt(2),:,:));
% 
% %make plots: 2D gg
% Fig=figure;
% set(Fig,'Position',[400 200 800 800])
% subplot(131);imagesc(abs(gg')); title(['abs(gg) [', num2str(slt),']']);axis image; 
% subplot(132);imagesc(real(gg')); title(['real(gg) [', num2str(slt),']']);axis image; 
% subplot(133);imagesc(imag(gg')); title(['imag(gg) [', num2str(slt),']']);axis image; 
% 
% % Fig=figure;
% % set(Fig,'Position',[400 200 800 800])
% % subplot(131);imagesc(abs(gg')); title(['abs(gg) [', num2str(slt),']']);axis image; hold on; h2 = imagesc(Stim); axis image; set(h2,'AlphaData',Stim*0.5);
% % subplot(132);imagesc(real(gg')); title(['real(gg) [', num2str(slt),']']);axis image; hold on; h2 = imagesc(Stim); axis image; set(h2,'AlphaData',Stim*0.5);
% % subplot(133);imagesc(imag(gg')); title(['imag(gg) [', num2str(slt),']']);axis image; hold on; h2 = imagesc(Stim); axis image; set(h2,'AlphaData',Stim*0.5);

slt = [140,55];
gg = squeeze(GG0(slt(1),slt(2),:,:)); 
ggStim = gg(:,logical(stim));
ggRecover = gg(:,logical(1-(stim)));
pdi = squeeze(squeeze(PDI(slt(1),slt(2),3,:)));

% %fft analysis
% fCoor = linspace(-P.CCFR/2,P.CCFR/2,100);
% Fgg = fftshift(fft(gg));
% FggStim = Fgg(:,logical(stim));
% FggRecover = Fgg(:,logical(1-(stim)));
% Fig=figure;
% set(Fig,'Position',[400 200 800 800])
% subplot(131);imagesc(fCoor,[1:nRpt],abs(Fgg')); title(['abs(gg) [', num2str(slt),']']);%axis equal tight; 
% subplot(132);imagesc(fCoor,[1:nRpt],real(Fgg')); title(['real(gg) [', num2str(slt),']']);%axis equal tight; 
% subplot(133);imagesc(fCoor,[1:nRpt],imag(Fgg')); title(['imag(gg) [', num2str(slt),']']);%axis equal tight;
% COLOR=[1 0 0;
%     0.08 0.17 0.55
%     0.31 0.31 0.31];
% Color1.ShadeAlpha = 0.15;
% Color1.Shade = COLOR(1,:);
% Color1.Line = COLOR(1,:);
% Color2.ShadeAlpha = 0.15;
% Color2.Shade = COLOR(2,:);
% Color2.Line = COLOR(2,:);
% Fig=figure;
% set(Fig,'Position',[400 200 600 800])
% subplot(311); ShadedErrorbar(mean(abs(FggRecover),2)', std(abs(FggRecover),1,2)', fCoor, Color2); 
% hold on; ShadedErrorbar(mean(abs(FggStim),2)', std(abs(FggStim),1,2)', fCoor, Color1); title(['abs(gg) [', num2str(slt),']']);xlabel('Freq[Hz]');
% subplot(312); ShadedErrorbar(mean(real(FggRecover),2)', std(real(FggRecover),1,2)', fCoor, Color2); 
% hold on; ShadedErrorbar(mean(real(FggStim),2)', std(real(FggStim),1,2)', fCoor, Color1); title(['real(gg) [', num2str(slt),']']);xlabel('Freq[Hz]');
% subplot(313); ShadedErrorbar(mean(imag(FggRecover),2)', std(imag(FggRecover),1,2)', fCoor, Color2); 
% hold on; ShadedErrorbar(mean(imag(FggStim),2)', std(imag(FggStim),1,2)', fCoor, Color1); title(['imag(gg) [', num2str(slt),']']);xlabel('Freq[Hz]');
% 
% %fft: real imag abs gg
% Fig=figure;
% set(Fig,'Position',[400 200 600 800])
% subplot(411); ShadedErrorbar(mean(fftshift(fft(abs(ggRecover)),1),2)', std(fftshift(fft(abs(ggRecover)),1),1,2)', xCoor, Color2); 
% hold on; ShadedErrorbar(mean(fftshift(fft(abs(ggStim)),1),2)', std(fftshift(fft(abs(ggStim)),1),1,2)', xCoor, Color1); title(['abs(gg) [', num2str(slt),']']);xlabel('tau[ms]');
% subplot(412); ShadedErrorbar(mean(fftshift(fft(real(ggRecover)),1),2)', std(fftshift(fft(real(ggRecover)),1),1,2)', xCoor, Color2); 
% hold on; ShadedErrorbar(mean(fftshift(fft(real(ggStim)),1),2)', std(fftshift(fft(real(ggStim)),1),1,2)', xCoor, Color1); title(['real(gg) [', num2str(slt),']']);xlabel('tau[ms]');
% subplot(413); ShadedErrorbar(mean(fftshift(fft(imag(ggRecover)),1),2)', std(fftshift(fft(imag(ggRecover)),1),1,2)', xCoor, Color2); 
% hold on; ShadedErrorbar(mean(fftshift(fft(imag(ggStim)),1),2)', std(fftshift(fft(imag(ggStim)),1),1,2)', xCoor, Color1); title(['imag(gg) [', num2str(slt),']']);xlabel('tau[ms]');
% subplot(414); ShadedErrorbar(mean(abs(fftshift(fft(angle(ggRecover)),1)),2)', std(abs(fftshift(fft(angle(ggRecover)),1)),1,2)', xCoor, Color2); 
% hold on; ShadedErrorbar(mean(abs(fftshift(fft(angle(ggStim)),1)),2)', std(abs(fftshift(fft(angle(ggStim)),1)),1,2)', xCoor, Color1); title(['phase(gg) [', num2str(slt),']']);xlabel('tau[ms]');
% 
% ggn = gg./(gg(1,:));
% Fggn = abs(fftshift(fft(angle(ggn)),1));
% [~,pk] = max(Fggn(51:end,:),[],1);

%% make plots: stim vs recover gg
COLOR=[1 0 0;
    0.08 0.17 0.55
    0.31 0.31 0.31];
Color1.ShadeAlpha = 0.15;
Color1.Shade = COLOR(1,:);
Color1.Line = COLOR(1,:);
Color2.ShadeAlpha = 0.15;
Color2.Shade = COLOR(2,:);
Color2.Line = COLOR(2,:);
Fig=figure;
set(Fig,'Position',[400 200 600 800])
subplot(411); ShadedErrorbar(mean(abs(ggRecover),2)', std(abs(ggRecover),1,2)', xCoor, Color2); 
hold on; ShadedErrorbar(mean(abs(ggStim),2)', std(abs(ggStim),1,2)', xCoor, Color1); title(['abs(gg) [', num2str(slt),']']);xlabel('tau[ms]');
subplot(412); ShadedErrorbar(mean(real(ggRecover),2)', std(real(ggRecover),1,2)', xCoor, Color2); 
hold on; ShadedErrorbar(mean(real(ggStim),2)', std(real(ggStim),1,2)', xCoor, Color1); title(['real(gg) [', num2str(slt),']']);xlabel('tau[ms]');
subplot(413); ShadedErrorbar(mean(imag(ggRecover),2)', std(imag(ggRecover),1,2)', xCoor, Color2); 
hold on; ShadedErrorbar(mean(imag(ggStim),2)', std(imag(ggStim),1,2)', xCoor, Color1); title(['imag(gg) [', num2str(slt),']']);xlabel('tau[ms]');
subplot(414); ShadedErrorbar(mean(angle(ggRecover),2)', std(angle(ggRecover),1,2)', xCoor, Color2); 
hold on; ShadedErrorbar(mean(angle(ggStim),2)', std(angle(ggStim),1,2)', xCoor, Color1); title(['phase(gg) [', num2str(slt),']']);xlabel('tau[ms]');
 
% % loge
% Fig=figure;
% set(Fig,'Position',[400 200 600 800])
% subplot(311); ShadedErrorbar(mean(log(abs(ggRecover)),2)', std(log(abs(ggRecover)),1,2)', xCoor, Color2); 
% hold on; ShadedErrorbar(mean(log(abs(ggStim)),2)', std(log(abs(ggStim)),1,2)', xCoor, Color1); title(['abs(gg) [', num2str(slt),']']);xlabel('tau[ms]');
% subplot(312); ShadedErrorbar(mean(real(ggRecover),2)', std(real(ggRecover),1,2)', xCoor, Color2); 
% hold on; ShadedErrorbar(mean(real(ggStim),2)', std(real(ggStim),1,2)', xCoor, Color1); title(['real(gg) [', num2str(slt),']']);xlabel('tau[ms]');
% subplot(313); ShadedErrorbar(mean(imag(ggRecover),2)', std(imag(ggRecover),1,2)', xCoor, Color2); 
% hold on; ShadedErrorbar(mean(imag(ggStim),2)', std(imag(ggStim),1,2)', xCoor, Color1); title(['imag(gg) [', num2str(slt),']']);xlabel('tau[ms]');
% 
% Fig=figure;
% set(Fig,'Position',[400 200 600 800])
% subplot(311);  plot(abs(gg(:,logical(stim))),'-r'); hold on;plot(abs(gg(:,logical(1-(stim)))),'.b');title(['abs(gg) [', num2str(slt),']']);
% subplot(312);  plot(real(gg(:,logical(stim))),'-r');hold on;plot(real(gg(:,logical(1-(stim)))),'.b');title(['real(gg) [', num2str(slt),']']);
% subplot(313); plot(imag(gg(:,logical(stim))),'-r');hold on;plot(imag(gg(:,logical(1-(stim)))),'.b'); title(['imag(gg) [', num2str(slt),']']);

%% make plots: stim vs recover gg (for conference)
COLOR=[1 0 0;
    0.08 0.17 0.55
    0.31 0.31 0.31];
Color1.ShadeAlpha = 0.2;
Color1.Shade = COLOR(1,:);
Color1.Line = COLOR(1,:);
Color2.ShadeAlpha = 0.15;
Color2.Shade = COLOR(2,:);
Color2.Line = COLOR(2,:);
Fig=figure;
set(Fig,'Position',[400 200 250 600])
subplot(411); ShadedErrorbar(mean(abs(ggRecover),2)', std(abs(ggRecover),1,2)', xCoor, Color2); 
hold on; ShadedErrorbar(mean(abs(ggStim),2)', std(abs(ggStim),1,2)', xCoor, Color1); xlabel('\tau (ms)');ylabel('abs(g_{1})');
subplot(412); ShadedErrorbar(mean(real(ggRecover),2)', std(real(ggRecover),1,2)', xCoor, Color2); 
hold on; ShadedErrorbar(mean(real(ggStim),2)', std(real(ggStim),1,2)', xCoor, Color1); xlabel('\tau (ms)');ylabel('Re(g_{1})');
subplot(413); ShadedErrorbar(mean(imag(ggRecover),2)', std(imag(ggRecover),1,2)', xCoor, Color2); 
hold on; ShadedErrorbar(mean(imag(ggStim),2)', std(imag(ggStim),1,2)', xCoor, Color1); xlabel('\tau (ms)');ylabel('Im(g_{1})');
subplot(414); ShadedErrorbar(mean(angle(ggRecover),2)', std(angle(ggRecover),1,2)', xCoor, Color2); 
hold on; ShadedErrorbar(mean(angle(ggStim),2)', std(angle(ggStim),1,2)', xCoor, Color1); xlabel('\tau (ms)');ylabel('Phase(g_{1})');

%% III.single pixel g1 based fUS
% % interpolate gg
% [tau0,t0] = meshgrid(1:400,1:20);
% [tau1,t1] = meshgrid(1:400,0+20/100:20/100:20);
% gg1 = interp2(tau0,t0,gg,tau1,t1,'spline');

igg = imag(gg);
rgg = real(gg);
absgg = abs(gg);
absngg = absgg./absgg(1,:);
% ratiogg = absgg(1,:)./absgg(10,:);

% % % Smooth gg
% % igg = movmean(igg,3);
% % rgg = movmean(rgg,3);
% % absgg = movmean(absgg,3);
% % 
% % igg = sgolayfilt(double(igg),3,5);
% % rgg = sgolayfilt(double(rgg),3,5);
% % absgg = sgolayfilt(double(absgg),3,5);

% % 1. use time lag of first valley rgg
% peakrange = [1,50];
% [locmin,Index] = findpeakIndex(rgg,peakrange);
% 
% % 2. abs(gg(1)-gg(10))
% decayrate = absgg(1,:)-absgg(10,:);
% 
% % 3. sum absgg all
% sumrange = [1,100];
% sumgg = sum(abs(absngg(sumrange(1):sumrange(2),:)),1);
% 
% % 4. sum absngg 1:index
% for k = 1: nRpt
%     index = squeeze(squeeze(Index(1,k)));
%     sum0 = sum(abs(absngg(sumrange(1): index, k)),1);
%     sumgg(1,k) = sum0;
% end
% % Index = sum(abs(absgg(sumrange(1):sumrange(2),:).*repmat(log10([sumrange(1):sumrange(2)])',[1,nRpt])),1);
% 
% % plot g1-fUS vs. pdi 
% Fig=figure;
% set(Fig,'Position',[500 500 600 400])
% subplot(211);plot(-Index); title(['real(gg) locmin(',num2str(peakrange),')  [', num2str(slt),']']);
% subplot(212);plot(pdi);title(['PDI [', num2str(slt),']']);
% 
% % Fig=figure;
% % set(Fig,'Position',[500 500 600 400])
% % subplot(211);plot(-movmean(Index,3)); title(['imag(gg) locmin(',num2str(minrange),')  [', num2str(slt),']']);
% % subplot(212);plot(movmean(pdi,3));title(['PDI [', num2str(slt),']']);
% Fig=figure;
% set(Fig,'Position',[500 500 600 400])
% subplot(211);plot(-Index); title(['abs(gg) sum(',num2str(sumrange),')  [', num2str(slt),']']);
% subplot(212);plot(pdi);title(['PDI [', num2str(slt),']']);

%% check PC for single pixel
% trialgg = reshape(real(gg(1,trial.nBase-5+1: end-5)), [trial.nlength,trial.n]);
% mtrialgg = median(trialgg,2);

% fgg = absgg(1,trial.nBase-5+1: end-5)-absgg(10,trial.nBase-5+1: end-5);
% trialgg = reshape(fgg, [trial.nlength,trial.n]);
% mtrialgg = median(trialgg,2);

% [B,A]=butter(4,200/(1000/2),'low');    %coefficients for the high pass filter
% absgg1(101:100+ntau,:)=absgg;
% absgg1(1:100,:)=flip(absgg1(101:200,:),1);
% absgg2=filter(B,A,absgg1,[],1);    
% absggHP=[absgg2(103:end,:); absgg2(1:2,:)]; 

% absggHP = medfilt1(absgg, 3, [], 1);
absggHP = movmean(absgg,3,1);

fgg = sqrt(abs(-(log(max(absggHP(1:10,:),[],1))-log(min(absggHP(1:10,:),[],1)))));%
trialgg = reshape(fgg(1,trial.nBase-5+1: end-5), [trial.nlength,trial.n]);
mtrialgg = median(trialgg,2);

CorrMap(mtrialgg,stimhrf)

trialpdi = reshape(medfilt1(pdi(trial.nBase-5+1: end-5,1),3), [trial.nlength,trial.n]);
mtrialpdi = median(trialpdi,2);

CorrMap(mtrialpdi,stimhrf)

ratiotrialgg = trialgg./mean(trialgg(1:5,:))*100;
ratiotrialpdi = trialpdi./mean(trialpdi(1:5,:))*100;
mratiotrialgg = median(ratiotrialgg,2);
mratiotrialpdi = median(ratiotrialpdi,2);

Fig=figure;
set(Fig,'Position',[500 500 600 400])
subplot(211);yyaxis left; plot(1:length(pdi),abs(fgg)); title(['gg  [', num2str(slt),']']);yyaxis right; plot(stim)
subplot(212);yyaxis left;plot(1:length(pdi),pdi);title(['PDI [', num2str(slt),']']);yyaxis right; plot(stim)
Fig=figure;
set(Fig,'Position',[500 500 600 400])
subplot(211);plot(abs(trialgg)); title(['gg  [', num2str(slt),']']);
subplot(212);plot(trialpdi);title(['PDI [', num2str(slt),']']);
figure; plot(mratiotrialgg); hold on; plot(mratiotrialpdi); 
hold on; plot(6:10,ones(5,1)*min(mratiotrialpdi(:)),'-k','LineWidth',2);legend({'gg','pdi'})

% figure; yyaxis left; plot(mtrialgg); hold on; yyaxis right; plot(mtrialpdi); 
% hold on; plot(6:10,ones(5,1)*min(mtrialpdi(:)),'-k','LineWidth',2);legend({'gg','pdi'})

%% IV.whole image g1 based fUS correlation
% % interploate GG 
% GG1 = zeros(nz,nx,100,nRpt);
% tic
% for i = 1:nz
%     for j = 1:nx
%         GG1(i,j,:,:) = interp2(tau0,t0,squeeze(squeeze(GG0(i,j,:,:))),tau1,t1,'spline');
%     end
% end
% toc

% % average trials first
% GG1 = reshape(GG0(:,:,:,trial.nBase-5+1: end-5), [nz,nx,ntau,trial.nlength,trial.n]);
% GG2 = median(GG1,5);
% clear GG1;

peakrange = [1,100];
sumrange = [1,100];

iGG = imag(GG0);
rGG = real(GG0);
absGG = abs(GG0);
phaseGG = angle(GG0);
absnGG = absGG./absGG(:,:,1,:);
rnGG = rGG./rGG(:,:,1,:);

% absGGHP = medfilt1(absGG, 3, [], 3);
absGGHP = movmean(absGG, 3, 3);
% [B,A]=butter(4,200/(1000/2),'low');    %coefficients for the high pass filter
% absGG1(:,:,101:100+ntau,:)=absGG;
% absGG1(:,:,1:100,:)=flip(absGG1(:,:,101:200,:),3);
% absGG2=filter(B,A,absGG1,[],3);    
% absGGHP=cat(3, absGG2(:,:,103:end,:), absGG2(:,:,1:2,:)); 

mphaseGG = mean(angle(GG0),4);
% figure; imagesc(max(abs(mphaseGG),[],3)); axis image; colorbar; 
vthre = 0.8;%0.5;
vmsk = max(abs(mphaseGG),[],3)>vthre;
% figure; imagesc(vmsk);axis image;title(['Threshold: >',num2str(vthre)]);

% % %% 1. use time lag of first valley rGG
% % IndexGG = findpeakIndex(rGG,peakrange);
% % IndexGG = squeeze(-IndexGG);
% % 
% % % [locmin,IndexGG] = findpeakIndex(iGG,peakrange); % find pi/2 max
% % % DecayGG = squeeze(abs(iGG(:,:,1,:)-locmin)./(IndexGG-0));
% % 
% % %% 2. sum absnGG 1:index 
% % [locmin,IndexGG] = findpeakIndex(iGG,peakrange);% find pi/2 max iGG
% % % IndexGG = round(movmean(IndexGG,3,4));
% % [locmin,IndexGG] = findpeakIndex(iGG,peakrange);% find 3*pi/2 min iGG
% % % IndexGG = round(movmean(IndexGG,3,4));
% % % IndexGG = round(medfilter(IndexGG,3,4));
% % tic
% % for m = 1:nz
% %     for n = 1:nx
% %         for k = 1: nRpt
% %             index = squeeze(squeeze(IndexGG(m,n,1,k)));
% %             sum0 = sum(abs(absnGG(m,n,sumrange(1): index, k)),3);
% %             sumGG(m,n,k) = sum0;
% %         end
% %     end
% % end
% % toc
% % 
% % %% 3. sum absnGG/ decayGG (index based on max phase)
% % dir = sum(mphaseGG(:,:,1:10),3);
% % dir(dir>0) = 1;
% % dir(dir<0) = -1;
% % mphaseGG_bar = mphaseGG.*dir;
% % [Locmin, Index] = max(mphaseGG_bar,[],3);
% % tic
% % for m = 1:nz
% %     for n = 1:nx
% %         index = squeeze(squeeze(Index(m,n)));
% %         sum0 = sum(abs(absnGG(m,n,sumrange(1): index,:)),3);
% %         sumGG(m,n,:) = sum0;
% %     end
% % end
% % toc

%% 4. sum absngg (index based on imag pi)
dir = sum(mphaseGG(:,:,1:10),3);
dir(dir>0) = 1;
dir(dir<0) = -1;
mphaseGG_bar = mphaseGG.*dir;
[Locmin, Index] = max(mphaseGG_bar,[],3);
zeroGG = iGG(:,:,4:99,:).*iGG(:,:,5:100,:);
zeropGG = find(zeroGG<0);
tic
for m = 1:nz
    for n = 1:nx
         [zr,zc] = find(squeeze(squeeze(zeroGG(m,n,:,:)))<0);
         zc1 = diff([0;zc]);
         zr1 = zr(find(zc1>0))+3;
         zr1 = interp1(zc(find(zc1>0)),zr1,[1:nRpt]);
         clear zc;
         vp = max(abs(mphaseGG(m,n,:)),[],3);
%          if vp < vthre
%          for k = 1: nRpt             
%          sum0 = sum(abs(absnGG(m,n,sumrange(1): round(zr1(k)),k)),3);
%          sumGG(m,n,k) = -sum0; % rCBFv
%          sumGGV(m,n,k) = -sum0; % rCBV
%          sumGGFV(m,n,k) = -sum0; % mixture of rCBFv, rCBV
%          end
%          else
             index = squeeze(squeeze(Index(m,n)));
             % sum0 = sum(abs(absnGG(m,n,sumrange(1): index,:)),3);
            
             % sum0 = squeeze(sqrt(-log(absGG(m,n,1,:)))-sqrt(-log(absGG(m,n,10,:))));
             % sum0 = squeeze(absGG(m,n,1,:)./(1-absGG(m,n,1,:)));
             % sum0 = squeeze(sqrt(abs(-log(medfilt1((absGG(m,n,1,:))./(absGG(m,n,5,:)),3,[],4)))));
             sum0 = squeeze(sqrt(abs(-log(max(absGGHP(m,n,5,:),[],3)./min(absGGHP(m,n,1,:),[],3)))));
             sum1 = squeeze(absGG(m,n,1,:)./(1-absGG(m,n,1,:)));
             sum2 = squeeze((rGG(m,n,1,:)-min(rGG(m,n,10,:),[],3)));
             sumGG(m,n,:) = sum0;
             sumGGV(m,n,:) = sum1;
             sumGGFV(m,n,:) = sum2;
%          end
    end
end
toc

% %% 5. GG(1)-GG(10)
% DecayGG = squeeze(absGG(:,:,1,:)-absGG(:,:,10,:));
% DecayGG = squeeze(absGG(:,:,1,:)-min(absGG(:,:,1:10,:),[],3)); % gg(1)-gg(10)
% 
% %% 6. sum absGG all
% sumGG = sum(abs(absGG(:,:,sumrange(1):sumrange(2),:)),3);
% sumGG = squeeze(-sumGG); % sum GG / nGG
% 
% ratiotau0 = ([sumrange(1):sumrange(2)]').^2;
% ratiotau0(50:sumrange(2)) = 1;
% ratiotau = permute(repmat(repmat(ratiotau0,[1,nRpt]),[1,1,nz,nx]),[3,4,1,2]);
% sumGG = sum(abs(absGG(:,:,sumrange(1):sumrange(2),:).*ratiotau),3);
% sumGG = squeeze(-sumGG); % weighted sum
%%
save([FilePath,'GG_invivo_results.mat'],'sumGG', 'sumGGV', 'sumGGFV', 'PDI', 'eqNoise', 'trial', 'dir', 'vmsk', 'nx', 'nz', 'ntau');

% %% Correlation map
% % trialGG = reshape(DecayGG(:,:,trial.nBase-5+1: end-5), [nz,nx,trial.nlength,trial.n]);
% % mtrialGG = median(trialGG,4);
% % [actmap, coefmap] = CoorCoeffMap(mtrialGG, stimhrf', 0); 
% 
% trialsumGG = reshape(sumGG(:,:,trial.nBase-5+1: end-5), [nz,nx,trial.nlength,trial.n]);
% mtrialsumGG = median(trialsumGG,4);
% [actmap_GG0, coefmap_GG0] = CoorCoeffMap(-mtrialsumGG, stimhrf', 0);
% 
% trialsumGGV = reshape(sumGGV(:,:,trial.nBase-5+1: end-5), [nz,nx,trial.nlength,trial.n]);
% mtrialsumGGV = median(trialsumGGV,4);
% [actmap_GGV0, coefmap_GGV0] = CoorCoeffMap(-mtrialsumGGV, stimhrf', 0);
% 
% trialsumGGFV = reshape(sumGGFV(:,:,trial.nBase-5+1: end-5), [nz,nx,trial.nlength,trial.n]);
% mtrialsumGGFV = median(trialsumGGFV,4);
% [actmap_GGFV0, coefmap_GGFV0] = CoorCoeffMap(-mtrialsumGGFV, stimhrf', 0);
% %% ROI averaged time course
% newSlect = 1;
% if newSlect==1
%     figure; imagesc(coefmap_GG0); axis image;colormap(jet);caxis([0,1]);
%     [loc_x,loc_y]=ginput(6); % [x z]
%     BW=roipoly(coefmap_GG0,loc_x,loc_y);
% end
% Fig = figure; set(Fig,'Position',[600 600 500 350])
% imagesc(BW.*coefmap_GG0); axis image;caxis([-1,1]); colorbar;colormap('jet'); title('ROI');
% 
% roiGG = squeeze(sum(sum(sumGG.*repmat(BW.*vmsk,[1, 1,nRpt]),1),2));
% % Fig = figure;set(Fig,'Position',[400 400 800 170]) 
% % yyaxis left;plot(-roiGG, 'b', 'LineWidth', 1);ylabel('G1-fUS(CBF)'); xlabel('Time[s]'); 
% % yyaxis right; plot(stim);grid on;
% 
% newSlect = 1;
% if newSlect==1
%     figure; imagesc(coefmap_GGV0); axis image;colormap(jet);caxis([0,1]);
%     [loc_x,loc_y]=ginput(6); % [x z]
%     BW=roipoly(coefmap_GGV0,loc_x,loc_y);
% end
% Fig = figure; set(Fig,'Position',[600 600 500 350])
% imagesc(BW.*coefmap_GGV0); axis image;caxis([-1,1]); colorbar;colormap('jet'); title('ROI');
% 
% roiGGV = squeeze(sum(sum(sumGGV.*repmat(BW.*vmsk,[1, 1,nRpt]),1),2));
% % Fig = figure;set(Fig,'Position',[400 400 800 170]) 
% % yyaxis left;plot(-roiGGV, 'b', 'LineWidth', 1);ylabel('G1-fUS(CBV)'); xlabel('Time[s]'); 
% % yyaxis right; plot(stim);grid on;
% 
% newSlect = 1;
% if newSlect==1
%     figure; imagesc(coefmap_GGFV0); axis image;colormap(jet);caxis([0,1]);
%     [loc_x,loc_y]=ginput(6); % [x z]
%     BW=roipoly(coefmap_GGFV0,loc_x,loc_y);
% end
% Fig = figure; set(Fig,'Position',[600 600 500 350])
% imagesc(BW.*coefmap_GGFV0); axis image;caxis([-1,1]); colorbar;colormap('jet'); title('ROI')
% 
% roiGGFV = squeeze(sum(sum(sumGGFV.*repmat(BW.*vmsk,[1, 1,nRpt]),1),2));
% % Fig = figure;set(Fig,'Position',[400 400 800 170]) 
% % yyaxis left;plot(-roiGGFV, 'b', 'LineWidth', 1);ylabel('G1-fUS(novel)'); xlabel('Time[s]'); 
% % yyaxis right; plot(stim);grid on;
% %% coefmap use ROI based HRF
% % figure(2); imagesc(coefmap); axis image;colormap(jet);caxis([0,1]);
% HRF = calROIHRF(-sumGG,trial,BW.*vmsk);
% [actmap_GG, coefmap_GG] = CoorCoeffMap(-mtrialsumGG, HRF, 0);
% 
% HRF = calROIHRF(-sumGGV,trial,BW.*vmsk);
% [actmap_GGV, coefmap_GGV] = CoorCoeffMap(-mtrialsumGGV, HRF, 0);
% 
% HRF = calROIHRF(-sumGGFV,trial,BW.*vmsk);
% [actmap_GGFV, coefmap_GGFV] = CoorCoeffMap(-mtrialsumGGFV, HRF, 0);
% 
% %% V.PDI comparison
% trialPDI = reshape(squeeze(PDI(:,:,3,trial.nBase-5+1: end-5)), [nz,nx,trial.nlength,trial.n]);
% mtrialPDI = median(trialPDI,4);
% [actmap_PDI0, coefmap_PDI0] = CoorCoeffMap(mtrialPDI, stimhrf', 0); 
% 
% % ROI averaged time course
% roiPDI = squeeze(sum(sum(squeeze(PDI(:,:,3,:)).*repmat(BW,[1, 1,nRpt]),1),2));
% % yyaxis left;plot(roiPDI, 'b', 'LineWidth', 1);ylabel('PDI'); xlabel('Time[s]'); 
% % yyaxis right; plot(stim);grid on;
% 
% % coefmap use ROI based HRF
% HRF = calROIHRF(squeeze(PDI(:,:,3,:)),trial,BW.*vmsk);
% [actmap_PDI, coefmap_PDI] = CoorCoeffMap(mtrialPDI, HRF, 0); 
% 
% 
% %% trial averaged time course GG vs. PDI
% trialroiGG = reshape(medfilt1(-roiGG(trial.nBase-trial.nRest+1: end-trial.nRest),3), [trial.nlength,trial.n]);
% trialroiPDI = reshape(medfilt1(roiPDI(trial.nBase-trial.nRest+1: end-trial.nRest),3), [trial.nlength,trial.n]);
% ratiotrialroiGG = trialroiGG./mean(trialroiGG(1:trial.nRest,:))*100;
% ratiotrialroiPDI = trialroiPDI./mean(trialroiPDI(1:trial.nRest,:))*100;
% mratiotrialroiGG = median(ratiotrialroiGG,2);
% stdratiotrialroiGG = std(ratiotrialroiGG,1,2);
% semratiotrialroiGG = std(ratiotrialroiGG,1,2)./trial.n;
% mratiotrialroiPDI = median(ratiotrialroiPDI,2);
% stdratiotrialroiPDI = std(ratiotrialroiPDI,1,2);
% semratiotrialroiPDI = std(ratiotrialroiPDI,1,2)./trial.n;
%     
% trialroiGGV = reshape(medfilt1(-roiGGV(trial.nBase-trial.nRest+1: end-trial.nRest),3), [trial.nlength,trial.n]);
% ratiotrialroiGGV = trialroiGGV./mean(trialroiGGV(1:trial.nRest,:))*100;
% mratiotrialroiGGV = median(ratiotrialroiGGV,2);
% stdratiotrialroiGGV = std(ratiotrialroiGGV,1,2);
% semratiotrialroiGGV = std(ratiotrialroiGGV,1,2)./trial.n;
% 
% trialroiGGFV = reshape(medfilt1(-roiGGFV(trial.nBase-trial.nRest+1: end-trial.nRest),3), [trial.nlength,trial.n]);
% ratiotrialroiGGFV = trialroiGGFV./mean(trialroiGGFV(1:trial.nRest,:))*100;
% mratiotrialroiGGFV = median(ratiotrialroiGGFV,2);
% stdratiotrialroiGGFV = std(ratiotrialroiGGFV,1,2);
% semratiotrialroiGGFV = std(ratiotrialroiGGFV,1,2)./trial.n;
% %% plots
% % stimhrf correlation map 
% Fig = figure;set(Fig,'Position',[200 200 1000 800])
% subplot(421)
% imagesc(coefmap_GG0); axis image; caxis([-1,1]); colormap('jet');colorbar;title('Correlation map (g1-CBFv)');
% subplot(422)
% imagesc(actmap_GG0); axis image; caxis([-1,1]); colormap('jet');colorbar;title('Activation map (Zscore>1.6)');
% subplot(423)
% imagesc(coefmap_GGV0); axis image; caxis([-1,1]); colormap('jet');colorbar;title('Correlation map (g1-CBV)');
% subplot(424)
% imagesc(actmap_GGV0); axis image; caxis([-1,1]); colormap('jet');colorbar;
% subplot(425)
% imagesc(coefmap_PDI0); axis image; caxis([-1,1]); colormap('jet');colorbar;title('Correlation map (PDI)');
% subplot(426)
% imagesc(actmap_PDI0); axis image; caxis([-1,1]); colormap('jet');colorbar;
% subplot(427)
% imagesc(coefmap_GGFV0); axis image; caxis([-1,1]); colormap('jet');colorbar;title('Correlation map (novel)');
% subplot(428)
% imagesc(actmap_GGFV0); axis image; caxis([-1,1]); colormap('jet');colorbar;
% 
% 
% % roihrf correlation map
% Fig = figure;set(Fig,'Position',[200 200 1000 800])
% subplot(421)
% imagesc(coefmap_GG); axis image; caxis([-1,1]); colormap('jet');colorbar;title('Correlation map (g1-CBF)');
% subplot(422)
% imagesc(actmap_GG); axis image; caxis([-1,1]); colormap('jet');colorbar;title('Activation map (Zscore>1.6)');
% subplot(423)
% imagesc(coefmap_GGV); axis image; caxis([-1,1]); colormap('jet');colorbar;title('Correlation map (g1-CBV)');
% subplot(424)
% imagesc(actmap_GGV); axis image; caxis([-1,1]); colormap('jet');colorbar;
% subplot(425)
% imagesc(coefmap_PDI); axis image; caxis([-1,1]); colormap('jet');colorbar;title('Correlation map (g1-PDI)');
% subplot(426)
% imagesc(actmap_PDI); axis image; caxis([-1,1]); colormap('jet');colorbar;
% subplot(427)
% imagesc(coefmap_GGFV); axis image; caxis([-1,1]); colormap('jet');colorbar;title('Correlation map (novel)');
% subplot(428)
% imagesc(actmap_GGFV); axis image; caxis([-1,1]); colormap('jet');colorbar;
% 
% % time course
% Fig = figure;set(Fig,'Position',[400 400 700 500])
% subplot(411)
% yyaxis left;plot(-roiGG, 'b', 'LineWidth', 1);ylabel('g1-fUS(CBF)'); xlabel('Time[s]'); 
% yyaxis right; plot(stim);grid on;
% subplot(412)
% yyaxis left;plot(-roiGGV, 'b', 'LineWidth', 1);ylabel('g1-fUS(CBV)'); xlabel('Time[s]'); 
% yyaxis right; plot(stim);grid on;
% subplot(413)
% yyaxis left;plot(roiPDI, 'b', 'LineWidth', 1);ylabel('PDI'); xlabel('Time[s]'); 
% yyaxis right; plot(stim);grid on;
% subplot(414)
% yyaxis left;plot(-roiGGFV, 'b', 'LineWidth', 1);ylabel('g1-fUS(novel)'); xlabel('Time[s]'); 
% yyaxis right; plot(stim);grid on;
% 
% %% averaged time course compare
% disp('Ratio of Peak rCBF / Peak rCBV: ')
% R = (max(movmean(mratiotrialroiGG, 3)/100)*max(movmean(mratiotrialroiGGV, 3)/100)-1)/(max(movmean(mratiotrialroiGGV, 3)/100)-1)
% 
% Fig = figure; set(Fig,'Position',[400 600 450 300])
% % plot(movmean(mratiotrialroiGG,3)); 
% % hold on; plot(movmean(mratiotrialroiGGV,3));
% % hold on; plot(movmean(mratiotrialroiPDI,3));  
% plot(mratiotrialroiGG); 
% hold on; plot(mratiotrialroiGGV);
% hold on; plot(mratiotrialroiPDI);
% hold on; plot(mratiotrialroiGGFV);
% hold on; 
% plot(trial.nRest:trial.nRest+trial.nStim-1,ones(trial.nStim,1)*min(mratiotrialroiPDI(:)),'-k','LineWidth',2);
% legend({'g1-CBF','g1-CBV','pdi','novel'})
% title(['rCBFlow/rCBV = ', num2str(R)])
% xlabel('Time[s]');
% ylabel('%');
% ylim([90 120])

