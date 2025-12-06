
clear all; 
addpath D:\g1_based_fUS\A-US-fUS-g1\SubFunctions\

%% I.load GG data
[FileName,FilePath]=uigetfile('V:\R-Functional1-20190503-PMP9-Awake\DATA\20190503-Whisker-CP6-Awake-Airpuff-Left\');
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
  GGn0 = myFile.GGnp(:,:,:,1);
  GGp0 = myFile.GGnp(:,:,:,2);
  PDIHP = myFile.PDIHP;
  eqNoise = myFile.eqNoise;
  disp([iFileName,' was loaded!'])
  else
   disp([iFileName, ' skipped!'])
  SkipFile(indSkipped)=iRpt;
   indSkipped=indSkipped+1;  
 end
 GG0(:,:,:,k) = GG;
 GGn(:,:,:,k) = GGn0;
 GGp(:,:,:,k) = GGp0;
 PDIHP0(:,:,:,k) = PDIHP;
 eqNoise0(:,:,k) = eqNoise;
 k = k+1;
end
eqNoise = mean(eqNoise0,3);
PDI = PDIHP0./eqNoise.^1.8;
ImgBG = mean(squeeze(PDI(:,:,3,:)),3).^0.35;
[nz,nx,ntau,~] = size(GG0);

trial.nBase = 25;
trial.n = 10;
trial.nStim = 5;
trial.nRecover = 25;
trial.nRest = 5; % included in nRecover
trial.nlength = trial.nStim + trial.nRecover;
% trial.parttern = Stim(trial.nBase-10+1:trial.nBase-10+(trial.nStim+trial.nRecover)); 

stim = zeros(nRpt,1);
for k = 1: trial.n % # of trials
stim(trial.nBase+(k-1)*trial.nlength+1: trial.nBase+trial.nStim+(k-1)*trial.nlength) = 1;
end
trial.stim = stim;
xCoor = (1:ntau)*1e3/P.CCFR; % tau axis[ms]

%% II.single pixel analysis
prompt={'New Select(0: No; 1: Yes)'};
name='Select';
defaultvalue={'1'};
numinput=inputdlg(prompt,name, 1, defaultvalue);
newSlect=str2num(numinput{1});

if newSlect==1
    figure(2); imagesc(ImgBG); axis image;colormap(gray);
    [slt(2), slt(1)]=ginput(1); % [x z]
    slt=round(slt);
end

gg = squeeze(GG0(slt(1), slt(2),:,:));

%make plots: 2D gg
Fig=figure;
set(Fig,'Position',[400 200 800 800])
subplot(131);imagesc(abs(gg')); title(['abs(gg) [', num2str(slt),']']);axis image; 
subplot(132);imagesc(real(gg')); title(['real(gg) [', num2str(slt),']']);axis image; 
subplot(133);imagesc(imag(gg')); title(['imag(gg) [', num2str(slt),']']);axis image; 

% Fig=figure;
% set(Fig,'Position',[400 200 800 800])
% subplot(131);imagesc(abs(gg')); title(['abs(gg) [', num2str(slt),']']);axis image; hold on; h2 = imagesc(Stim); axis image; set(h2,'AlphaData',Stim*0.5);
% subplot(132);imagesc(real(gg')); title(['real(gg) [', num2str(slt),']']);axis image; hold on; h2 = imagesc(Stim); axis image; set(h2,'AlphaData',Stim*0.5);
% subplot(133);imagesc(imag(gg')); title(['imag(gg) [', num2str(slt),']']);axis image; hold on; h2 = imagesc(Stim); axis image; set(h2,'AlphaData',Stim*0.5);

slt = [95,73];
gg = squeeze(GG0(slt(1),slt(2),:,:)); 
ggStim = gg(:,logical(stim));
ggRecover = gg(:,logical(1-(stim)));
pdi = squeeze(squeeze(PDI(slt(1),slt(2),3,:)));

%fft analysis
fCoor = linspace(-P.CCFR/2,P.CCFR/2,100);
Fgg = fftshift(fft(gg));
FggStim = Fgg(:,logical(stim));
FggRecover = Fgg(:,logical(1-(stim)));
Fig=figure;
set(Fig,'Position',[400 200 800 800])
subplot(131);imagesc(fCoor,[1:nRpt],abs(Fgg')); title(['abs(gg) [', num2str(slt),']']);%axis equal tight; 
subplot(132);imagesc(fCoor,[1:nRpt],real(Fgg')); title(['real(gg) [', num2str(slt),']']);%axis equal tight; 
subplot(133);imagesc(fCoor,[1:nRpt],imag(Fgg')); title(['imag(gg) [', num2str(slt),']']);%axis equal tight;
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
subplot(311); ShadedErrorbar(mean(abs(FggRecover),2)', std(abs(FggRecover),1,2)', fCoor, Color2); 
hold on; ShadedErrorbar(mean(abs(FggStim),2)', std(abs(FggStim),1,2)', fCoor, Color1); title(['abs(gg) [', num2str(slt),']']);xlabel('Freq[Hz]');
subplot(312); ShadedErrorbar(mean(real(FggRecover),2)', std(real(FggRecover),1,2)', fCoor, Color2); 
hold on; ShadedErrorbar(mean(real(FggStim),2)', std(real(FggStim),1,2)', fCoor, Color1); title(['real(gg) [', num2str(slt),']']);xlabel('Freq[Hz]');
subplot(313); ShadedErrorbar(mean(imag(FggRecover),2)', std(imag(FggRecover),1,2)', fCoor, Color2); 
hold on; ShadedErrorbar(mean(imag(FggStim),2)', std(imag(FggStim),1,2)', fCoor, Color1); title(['imag(gg) [', num2str(slt),']']);xlabel('Freq[Hz]');

%make plots: stim vs recover gg
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

% loge
Fig=figure;
set(Fig,'Position',[400 200 600 800])
subplot(311); ShadedErrorbar(mean(log(abs(ggRecover)),2)', std(log(abs(ggRecover)),1,2)', xCoor, Color2); 
hold on; ShadedErrorbar(mean(log(abs(ggStim)),2)', std(log(abs(ggStim)),1,2)', xCoor, Color1); title(['abs(gg) [', num2str(slt),']']);xlabel('tau[ms]');
subplot(312); ShadedErrorbar(mean(real(ggRecover),2)', std(real(ggRecover),1,2)', xCoor, Color2); 
hold on; ShadedErrorbar(mean(real(ggStim),2)', std(real(ggStim),1,2)', xCoor, Color1); title(['real(gg) [', num2str(slt),']']);xlabel('tau[ms]');
subplot(313); ShadedErrorbar(mean(imag(ggRecover),2)', std(imag(ggRecover),1,2)', xCoor, Color2); 
hold on; ShadedErrorbar(mean(imag(ggStim),2)', std(imag(ggStim),1,2)', xCoor, Color1); title(['imag(gg) [', num2str(slt),']']);xlabel('tau[ms]');

Fig=figure;
set(Fig,'Position',[400 200 600 800])
subplot(311);  plot(abs(gg(:,logical(stim))),'-r'); hold on;plot(abs(gg(:,logical(1-(stim)))),'.b');title(['abs(gg) [', num2str(slt),']']);
subplot(312);  plot(real(gg(:,logical(stim))),'-r');hold on;plot(real(gg(:,logical(1-(stim)))),'.b');title(['real(gg) [', num2str(slt),']']);
subplot(313); plot(imag(gg(:,logical(stim))),'-r');hold on;plot(imag(gg(:,logical(1-(stim)))),'.b'); title(['imag(gg) [', num2str(slt),']']);

%% make plots for conference
COLOR=[1 0 0;
    0.08 0.17 0.55
    0.31 0.31 0.31];
Color1.ShadeAlpha = 0.2;
Color1.Shade = COLOR(1,:);
Color1.Line = COLOR(1,:);
Color2.ShadeAlpha = 0.2;
Color2.Shade = COLOR(2,:);
Color2.Line = COLOR(2,:);
Fig=figure;
set(Fig,'Position',[400 200 250 600])
subplot(311); ShadedErrorbar(mean(abs(ggRecover),2)', std(abs(ggRecover),1,2)', xCoor, Color2); 
hold on; ShadedErrorbar(mean(abs(ggStim),2)', std(abs(ggStim),1,2)', xCoor, Color1); xlabel('\tau (ms)');ylabel('abs(g_{1})');
subplot(312); ShadedErrorbar(mean(real(ggRecover),2)', std(real(ggRecover),1,2)', xCoor, Color2); 
hold on; ShadedErrorbar(mean(real(ggStim),2)', std(real(ggStim),1,2)', xCoor, Color1); xlabel('\tau (ms)');ylabel('Re(g_{1})');
subplot(313); ShadedErrorbar(mean(imag(ggRecover),2)', std(imag(ggRecover),1,2)', xCoor, Color2); 
hold on; ShadedErrorbar(mean(imag(ggStim),2)', std(imag(ggStim),1,2)', xCoor, Color1); xlabel('\tau (ms)');ylabel('Im(g_{1})');


%% III.single pixel g1 based fUS
[tau0,t0] = meshgrid(1:400,1:20);
[tau1,t1] = meshgrid(1:400,0+20/100:20/100:20);
gg1 = interp2(tau0,t0,gg,tau1,t1,'spline');
igg = imag(gg);
rgg = real(gg);
absgg = abs(gg);

% igg = movmean(igg,3);
% rgg = movmean(rgg,3);
% absgg = movmean(absgg,3);
% 
% igg = sgolayfilt(double(igg),3,5);
% rgg = sgolayfilt(double(rgg),3,5);
% absgg = sgolayfilt(double(absgg),3,5);

peakrange = [1,50];
[locmin,Index] = findpeakIndex(rgg,peakrange);
% decayrate = rgg

decayrate = absgg(1,:)-absgg(10,:);

[locmin,Index] = findpeakIndex(igg,peakrange);
decayrate = squeeze(abs(igg(1,:)-locmin)./(Index-0));

sumrange = [1,100];
Index = sum(abs(absgg(sumrange(1):sumrange(2),:)),1);
Index = sum(abs(absgg(sumrange(1):sumrange(2),:).*repmat(log10([sumrange(1):sumrange(2)])',[1,nRpt])),1);

% Fig=figure;
% set(Fig,'Position',[500 500 600 400])
% subplot(211);plot(-Index); title(['imag(gg) first peak  [', num2str(slt),']']);
% subplot(212);plot(pdi);title(['PDI [', num2str(slt),']']);

Fig=figure;
set(Fig,'Position',[500 500 600 400])
subplot(211);plot(-Index); title(['real(gg) locmin(',num2str(peakrange),')  [', num2str(slt),']']);
subplot(212);plot(pdi);title(['PDI [', num2str(slt),']']);

% Fig=figure;
% set(Fig,'Position',[500 500 600 400])
% subplot(211);plot(-movmean(Index,3)); title(['imag(gg) locmin(',num2str(minrange),')  [', num2str(slt),']']);
% subplot(212);plot(movmean(pdi,3));title(['PDI [', num2str(slt),']']);
Fig=figure;
set(Fig,'Position',[500 500 600 400])
subplot(211);plot(-Index); title(['abs(gg) sum(',num2str(sumrange),')  [', num2str(slt),']']);
subplot(212);plot(pdi);title(['PDI [', num2str(slt),']']);

%% IV.whole image g1 based fUS correlation
% GG1 = zeros(nz,nx,100,nRpt);
% tic
% for i = 1:nz
%     for j = 1:nx
%         GG1(i,j,:,:) = interp2(tau0,t0,squeeze(squeeze(GG0(i,j,:,:))),tau1,t1,'spline');
%     end
% end
% toc
GG1 = reshape(GG0(:,:,:,trial.nBase-5+1: end-5), [nz,nx,ntau,trial.nlength,trial.n]);
GG2 = median(GG1,5);
clear GG1;

iGG = imag(GG0);
rGG = real(GG0);
absGG = abs(GG0);
phaseGG = angle(GG0);
absnGG = absGG./absGG(:,:,1,:);

mphaseGG = mean(angle(GG0),4);
figure; imagesc(max(abs(mphaseGG),[],3)); axis image; colorbar; 
vthre = 0.5;
vmsk = max(abs(mphaseGG),[],3)>vthre;
figure; imagesc(vmsk);axis image;title(['Threshold: >',num2str(vthre)]);

IndexGG = findpeakIndex(rGG,peakrange);
IndexGG = squeeze(-IndexGG);

% [locmin,IndexGG] = findpeakIndex(iGG,peakrange); % find pi/2 max
% DecayGG = squeeze(abs(iGG(:,:,1,:)-locmin)./(IndexGG-0));

peakrange = [1,100];
sumrange = [1,100];
[locmin,IndexGG] = findpeakIndex(iGG,peakrange);% find pi/2 max iGG
% IndexGG = round(movmean(IndexGG,3,4));
[locmin,IndexGG] = findpeakIndex(iGG,peakrange);% find 3*pi/2 min iGG
IndexGG = round(movmean(IndexGG,3,4));
IndexGG = round(medfilter(IndexGG,3,4));
tic
for m = 1:nz
    for n = 1:nx
        for k = 1: nRpt
            index = squeeze(squeeze(IndexGG(m,n,1,k)));
        sum0 = sum(abs(absnGG(m,n,sumrange(1): index, k)),3);
        sumGG(m,n,k) = sum0;
        end
    end
end
toc

DecayGG = squeeze(absGG(:,:,1,:)-absGG(:,:,10,:));
DecayGG = squeeze(absGG(:,:,1,:)-min(absGG(:,:,1:10,:),[],3)); % gg(1)-gg(10)

[locmin,IndexangGG] = findpeakIndex(angGG,peakrange); % find first peak
IndexangGG = round(movmean(IndexangGG,3,4));
DecayrGG = squeeze(rGG(:,:,1,:)-rGG(IndexangGG))./(IndexangGG-1);% dectect

sumGG = sum(abs(absGG(:,:,sumrange(1):sumrange(2),:)),3);
sumGG = squeeze(-sumGG); % sum GG / nGG

ratiotau0 = ([sumrange(1):sumrange(2)]').^2;
ratiotau0(50:sumrange(2)) = 1;
ratiotau = permute(repmat(repmat(ratiotau0,[1,nRpt]),[1,1,nz,nx]),[3,4,1,2]);
sumGG = sum(abs(absGG(:,:,sumrange(1):sumrange(2),:).*ratiotau),3);
sumGG = squeeze(-sumGG); % weighted sum



trialGG = reshape(DecayGG(:,:,trial.nBase-5+1: end-5), [nz,nx,trial.nlength,trial.n]);
mtrialGG = median(trialGG,4);

trialsumGG = reshape(sumGG(:,:,trial.nBase-5+1: end-5), [nz,nx,trial.nlength,trial.n]);
mtrialsumGG = median(trialsumGG,4);

% correlation map calculation
% hrf = hemodynamicResponse(1,[2 16 0.5 1 20 0]);
hrf = hemodynamicResponse(1,[1.5 10 0.5 1 20 0 16]);
Stim = zeros(trial.nlength,1);
Stim(5+1:5+trial.nStim,:)=1;
stimhrf = filter(hrf,1,Stim);
figure; plot(stimhrf);

HRF = calROIHRF(IndexGG,trial);
[actmap, coefmap] = CoorCoeffMap(mtrialGG, stimhrf', 0); 
[actmap, coefmap] = CoorCoeffMap(mtrialsumGG, stimhrf', 0); 

%% V.PDI comparison
trialPDI = reshape(squeeze(PDI(:,:,3,trial.nBase-5+1: end-5)), [nz,nx,trial.nlength,trial.n]);
mtrialPDI = median(trialPDI,4);

coefmap = CoorCoeffMap(mtrialPDI, stimhrf', 0); 




