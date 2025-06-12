function G1fUSmain(msktype, GGthre, GGFVthre, mskthre) 

if nargin<2
    GGthre = 5;
    GGFVthre = 10;
    mskthre = 0.8;
end

addpath D:\g1_based_fUS\A-US-fUS-g1\SubFunctions\

%% I.load GG data
[FileName,FilePath]=uigetfile('V:\G1based fUS data\');
fileInfo=strsplit(FileName(1:end-4),'-');
myFile=matfile([FilePath,FileName]);
P=myFile.P;
prompt={'Start Repeat', 'Number of Repeats'};
name='File info'; 
defaultvalue={'1','325','1'};
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

%% default trial defination 
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

%% IV.whole image g1 based fUS correlation
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

mphaseGG = mean(angle(GG0),4);
% figure; imagesc(max(abs(mphaseGG),[],3)); axis image; colorbar; 
vthre = mskthre;%0.5;
vmsk = max(abs(mphaseGG),[],3)>vthre;
% figure; imagesc(vmsk);axis image;title(['Threshold: >',num2str(vthre)]);


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
         switch msktype
             case "rigid"
                 sum0 = squeeze(sqrt(abs(-log(max(absGGHP(m,n,GGthre,:),[],3)./min(absGGHP(m,n,1,:),[],3)))));
                 sum1 = squeeze(absGG(m,n,1,:)./(1-absGG(m,n,1,:)));
                 sum2 = squeeze((rGG(m,n,1,:)-min(rGG(m,n,GGFVthre,:),[],3)));
                 sumGG(m,n,:) = sum0;
                 sumGGV(m,n,:) = sum1;
                 sumGGFV(m,n,:) = sum2;
             case "soft"
                 if vp < vthre
                     for k = 1: nRpt             
                         sum0 = sum(abs(absnGG(m,n,sumrange(1): round(zr1(k)),k)),3);
                         sumGG(m,n,k) = -sum0; % rCBFv
                         sumGGV(m,n,k) = -sum0; % rCBV
                         sumGGFV(m,n,k) = -sum0; % mixture of rCBFv, rCBV
                     end
                 else
                     index = squeeze(squeeze(Index(m,n)));
                     sum0 = squeeze(sqrt(abs(-log(max(absGGHP(m,n,GGthre,:),[],3)./min(absGGHP(m,n,1,:),[],3)))));
                     sum1 = squeeze(absGG(m,n,1,:)./(1-absGG(m,n,1,:)));
                     sum2 = squeeze((rGG(m,n,1,:)-min(rGG(m,n,GGFVthre,:),[],3)));
                     sumGG(m,n,:) = sum0;
                     sumGGV(m,n,:) = sum1;
                     sumGGFV(m,n,:) = sum2;
                 end
         end
    end
end
toc

%%
save([FilePath,'GG_invivo_results.mat'],'sumGG', 'sumGGV', 'sumGGFV', 'PDI', 'eqNoise', 'trial', 'dir', 'vmsk', 'nx', 'nz', 'ntau');
