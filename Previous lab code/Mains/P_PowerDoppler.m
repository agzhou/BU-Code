clear all;
tic
load('D:\CODE\Mains\DAQParameters.mat');
defaultpath=DAQInfo.savepath;
addpath('D:\CODE\Functions');
[FileName,FilePath]=uigetfile(defaultpath);  % read data of a small part of the brain cortex (IQR matrix)  
disp('Loading data...');
load ([FilePath, FileName]);
disp('Data loaded!');
% IQ=IQData{1};
dataRAW=IQ;
clear IQ0 IQData;
[nz, nx, nt]=size(dataRAW);
fileinfo=strsplit(FileName(1:end-4),'-');
Agl=str2num(fileinfo{2});
nAgl=str2num(fileinfo{3});
fCC=str2num(fileinfo{4});
nCC=str2num(fileinfo{5});
iSupFrame=str2num(fileinfo{6});
toc
%% data processing parameters
prompt={'SVD Rank (low):', ['SVD Rank (High):(Max Rank: ',num2str(nCC),')'],'High pass cutoff frequency (Hz)',...
    ['nCC_process (nCC total: ',num2str(nCC),')'], 'Image refine scale','Transducer center frequency (Mhz)','Element Pitch (mm)'};
name='Power Doppler data processing';
% defaultvalue={'15', '110', '40',num2str(nCC),'5', num2str(P.frequency)};
defaultvalue={'20', num2str(nCC), '60',num2str(nCC),'5', '15.625','0.1'};
numinput=inputdlg(prompt,name, 1, defaultvalue);
RankLow=str2num(numinput{1});
RankHigh=str2num(numinput{2});
DefCutFreq=str2num(numinput{3});
nCC_proc=str2num(numinput{4});
RefScale=str2num(numinput{5});     % image refine scale
fCenter=str2num(numinput{6});     % transducer center frequency
dx=str2num(numinput{7});          % transducer element pitch
dz=1.540/fCenter/2;             % axial samppling step, in mm, the default is 4 samppling points for one wavelength  
%% x,z coordinates after image refined
% xCoor=linspace(0,nx*dx,nx*RefScale);
% zCoor=linspace(0,nz*dz,nz*RefScale);
xCoor=P.xCoor;
zCoor=P.zCoor;

%% SVD process 1 (direct SVD use MATLAB)
% [nz,nx,nt]=size(IQR);
% S=reshape(IQR,[nz*nx,nt]);
% [UU,DD,VV]=svd(S);
% DD0=zeros(size(DD));
% DD0(:,21:end)=DD(:,21:end);
% sBlood=reshape(UU*DD0*VV',[nz,nx,nt]);
tic
%% SVD process 2 (eigen-to-SVD use MATLAB)
IQR=dataRAW(:,:,1:nCC_proc); % IQ data used for data processing
[nz,nx,nt]=size(IQR);
rank=[RankLow:RankHigh];
S=reshape(IQR,[nz*nx,nt]);
S_COVt=(S'*S);
[V,D]=eig(S_COVt); % V is the right singular Vector of S/eigenvector; D is the eigenvalue/square of Singular value
for it=1:nt 
    Ddiag(it)=abs(sqrt(D(it,it)));
end
Ddiag=20*log10(Ddiag/max(Ddiag)); % singular value in db
[Ddesc, Idesc]=sort(Ddiag,'descend');

for it=1:nt
    Vdesc(:,it)=V(:,Idesc(it)); % Vdesc is the right singluar matrix in SVD (has little numerical error)
end
Vrank=zeros(size(Vdesc));
Vrank(:,rank)=Vdesc(:,rank);
Vnoise=zeros(size(Vdesc));
Vnoise(:,end)=Vdesc(:,end);
UDelta=S*Vdesc;
sBlood0=reshape(UDelta*Vrank',[nz,nx,nt]);
% sBlood=sBlood0./repmat(std(abs(sBlood0),1,2),[1,nx]);
%%%% Noise equalization 
sNoise=reshape(UDelta*Vnoise',[nz,nx,nt]);
B=ones(100, 100);%[50,50]
% sNoiseMed=convn(abs(squeeze(mean(sNoise,3))),B,'same');
sNoiseMed=medfilt2(abs(squeeze(mean(sNoise,3))),[30,30],'symmetric');
sNoiseMedNorm=sNoiseMed/min(sNoiseMed(:));
sBlood=sBlood0./repmat(sNoiseMedNorm,[1,1,nt]);

%%

% C = corrcoef(abs(UU(:,1:600)));
% figure; imagesc(abs(C)); colorbar; colormap('jet');
% axis off
% %% dIQ Noise equalization
% Noisefit = fit([1:nz]',PDISVD(:,40),'poly2');
% fittedNoise = Noisefit([1:nz]');
% fittedNoise = fittedNoise./min(fittedNoise);
% eqNoise = repmat(fittedNoise,[1,nx]);
% figure; imagesc(eqNoise); axis image; colorbar;
% figure; plot(PDISVD(:,40));hold on; plot(PDISVD(:,40)./fittedNoise);
% sBlood=sBlood0./repmat(eqNoise,[1,1,nt]);
% %% SVD analysis %%%     
% % added by Bingxue Liu, 0409 2021
% figure;
% yyaxis left; plot(Ddesc);
% DDesc = cumsum(sort(sqrt(diag(D)),'descend'));
% yyaxis right; plot(20*log10(DDesc/max(DDesc)));
% hold on; %xline(rank(1));
% 
% figure;
% iVector = [1, 20, 50, 100, 150, 200];
% subplot(211);
% plot(real(Vdesc(:,iVector)),'LineWidth',1);title('Temperal Singular Vectors Real Part');
% legend(strsplit(num2str(iVector)),'Location','eastoutside');
% subplot(212);
% plot(imag(Vdesc(:,iVector)),'LineWidth',1);title('Temperal Singular Vectors Imag Part');
% legend(strsplit(num2str(iVector)),'Location','eastoutside');
% 
% % [UU,DD,VV]=svd(S);
% for i = 1:length(iVector)
%     Ui(:,:,i) = reshape(UU(:,iVector(i)),[nz,nx]);
%     figure;
%     imagesc((abs(Ui(:,:,i))).^0.4); axis image; colormap(gray);colorbar;
%     title([num2str(iVector(i)),' th Spatial Singular Vector']);
%     figure;
%     subplot(211);
%     plot(real(Vdesc(:,iVector(i))),'LineWidth',1);title('Temperal Singular Vectors Real Part');
%     legend(strsplit(num2str(iVector(i))),'Location','eastoutside');
%     subplot(212);
%     plot(imag(Vdesc(:,iVector(i))),'LineWidth',1);title('Temperal Singular Vectors Imag Part');
%     legend(strsplit(num2str(iVector(i))),'Location','eastoutside');
% end
% 
% % save spatial singular vectors as video
% figure;
% Moviie = VideoWriter([FilePath,strcat('svd-',FileName),'-AVI.avi']);
% Moviie.Quality = 100;
% Moviie.FrameRate = 5;
% open(Moviie);
% 
% for i = 1 : length(iVector)
%     Ui(:,:,i) = reshape(UU(:,iVector(i)),[nz,nx]);
% end
% cmin = min(min((log(abs(Ui(:,:,i))))));
% cmax = max(max((log(abs(Ui(:,:,i))))));
% for i = 1 : length(iVector)
%     % make a movie
%    imagesc(log(abs(Ui(:,:,i)))); axis image; colormap(gray);colorbar;
%     caxis([cmin, cmax]); 
%     hColorbar = colorbar;
%     set(hColorbar, 'Ticks', sort([hColorbar.Limits, hColorbar.Ticks]));
%     title([num2str(iVector(i)),' th Spatial Singular Vector']);
% %         text(20,20,['Time = ' num2str((i)*ImgInfo.t/nf),'s'],...
% %             'Position',[150,20],...
% %             'Units','pixels',...
% %             'FontSize',12,'Color',[0,1,0])
%         hold off;
% %    pause (1)
%     vframe = getframe(gcf); 
%     writeVideo(Moviie, vframe);
%    figure(gcf)
% end
% close(Moviie);
% %%%%%%%%
% 
% for i = 1:nt
% %     Vdesc0(:,i) = Vdesc(:,i).*conj(mean(Vdesc(:,i)))./abs(mean(Vdesc(:,i))).^2;
% Vdesc0(:,i) = (Vdesc(:,i)-min(Vdesc(:,i)))./(max(Vdesc(:,i))-min(Vdesc(:,i)));
% end
% FVdesc = fftshift(fft(Vdesc0,[],1),1);
% fCoor = linspace(-fCC/2,fCC/2,nt);
% 
% figure; imagesc([1:nt],fCoor, 20*log10(abs(FVdesc)/max(max(abs(FVdesc)))));colormap(jet);axis square;
% hold on; line([rank(1),rank(1)],[fCoor(1),fCoor(end)],'Color','black','LineWidth',1);
% hold on; line([1,nt],[-DefCutFreq, -DefCutFreq],'Color','black','LineWidth',1);
% hold on; line([1,nt],[DefCutFreq, DefCutFreq],'Color','black','LineWidth',1);
% xlabel('iNumber of Singular Vectors'); ylabel('Frequency [Hz]'); title('Spectrum of Temperal Singular Vectors');

%% high pass filter
[B,A]=butter(4,DefCutFreq/fCC*2,'high');    %coefficients for the high pass filter
IQR1(:,:,21:20+nCC_proc)=IQR-repmat(IQR(:,:,1),[1,1,nCC_proc]);
for iCC=1:20
    IQR1(:,:,iCC)=IQR1(:,:,41-iCC);
end
sBloodHP=filter(B,A,IQR1,[],3);    % blood signal (filtering in the time dimension)
sBloodHP=sBloodHP(:,:,21:end)./repmat(sNoiseMedNorm,[1,1,nt]);           % the first 4 temporal samples are eliminates (filter oscilations)
%% power doppler
PDISVD=mean(abs(sBlood).^2,3); 
PDISVDdb=10*log10(PDISVD./max(PDISVD(:))); % SVD-based PD image in dB
PDIHP=mean(abs(sBloodHP).^2,3); 
PDIHPdb=10*log10(PDIHP./max(PDIHP(:)));% High pass filter-based PD image in dB
% PDISVDRef=imresize(PDISVDdb,RefScale);
% PDIHPRef=imresize(PDIHPdb,RefScale);
toc
% %% Frequency Analysis
% % modified by Bingxue Liu 04/10/2021
%     figure; imagesc(log(PDISVD)); axis image; colorbar;colormap(gray); 
%     [slt(2), slt(1)]=ginput(1); % [x z]
%     slt=round(slt);
% 
%   
%     sIQ= sBlood;
%     
%     
%     IQAVG=squeeze(IQ(slt(1),slt(2),:));
%     nt=length(IQAVG);
%     fIQ=fftshift(fft(IQAVG));
%     
%     sIQAVG=squeeze(sBlood(slt(1),slt(2),:));
%     nt=length(sIQAVG);
%     FsIQ=fftshift(fft(sIQAVG));
%     
%     Pneg=sum(abs(FsIQ(1:floor(end/2)))-median(abs(FsIQ)));
%     Ppos=sum(abs(FsIQ(floor(end/2)+1:end))-median(abs(FsIQ)));
%     rNeg=Pneg/(Pneg+Ppos);
%     rPos=Ppos/(Pneg+Ppos);
%     tCoor=linspace(1,nt,nt)/fCC*1000;
%     
%     figure;
%     subplot(2,1,1),plot(fCoor,abs(fIQ),'k')
%     hold on; plot(fCoor,abs(FsIQ),'r')
% %     xlim([-800 800])
%     ylim([0 mean(abs(fIQ))])
%     xlabel('Frequency [Hz]')
%     ylabel('Power Spectrum Density')
%     legend({'Raw','SVD'})
%     title(['Pneg=',num2str(rNeg),'; Ppos=',num2str(rPos)])
%     
%     subplot(2,1,2),plot(tCoor,(angle(IQAVG)),'k');
%     hold on; plot(tCoor,(angle(sIQAVG)),'r');
%     title('phase(IQ)');% ylim([0,1])
%     xlabel('time lag, [ms]')
%     legend({'Raw','SVD'})
%     

%% figure plot
% figure,
% imagesc(PDIdb);
% caxis([-22 0]);
% colormap(hot);
% colorbar;
fig=figure;
set(fig,'Position',[300 400 900 300]);
subplot(1,2,1);imagesc(xCoor, zCoor, PDISVDdb);
caxis(MyCaxis(PDISVDdb));
colormap(hot);
colorbar;
title ('SVD-based PDI')
xlabel('X [mm]'); ylabel('Z [mm]');
axis equal tight
subplot(1,2,2);imagesc(xCoor, zCoor, PDIHPdb);
caxis(MyCaxis(PDIHPdb));
colormap(hot);
colorbar;
title ('HP-based PDI')
xlabel('X [mm]'); ylabel('Z [mm]');
axis equal tight

