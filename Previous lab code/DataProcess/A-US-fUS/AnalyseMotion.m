% clear all;
addpath 'D:\CODE\Mains\A-US-fUS\SubFunctions'
addpath 'D:\CODE\DataProcess\BingxueLiu\Code fUS DATA analysis\SubFunctions'

[FileName,FilePath]=uigetfile('G:\0622_BL3_sedated\');%'H:\test_pupil_accel\15mintest');
load([FilePath,FileName]);
ImgInfo.savePathName=[FilePath,FileName];
ImgInfo.t = 500/1000;
fUSmovie(log(abs(IQ)),ImgInfo);
mIQ1 = squeeze(mean(mean(IQ,1),2));
 
SignalRank = [21,200];
[sIQ, Noise]=SVDfilter(IQ,SignalRank);
ImgInfo.savePathName=[FilePath,strcat('s',FileName)];
% fUSmovie(log(abs(sIQ)),ImgInfo);
msIQ  = squeeze(mean(mean(sIQ,1),2));
figure; imagesc(log(mean(abs(sIQ).^2,3)./Noise.^1.8));axis image; colorbar;colormap(hot);
title(['SVD: ',num2str(SignalRank),' Noise: end-50']);