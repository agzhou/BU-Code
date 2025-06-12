% use this program right after Power Doppler Processing toward last frame
% use the last frame's noise to equialize all the frames;
% Transform sIQ to PDI across frames;
%
% Last modified: 1/5/2021, Bingxue Liu

clear sIQ PDIm PDImdb;

%% load file path and name
load('D:\CODE\Mains\DAQParameters.mat');
defaultpath=DAQInfo.savepath;
addpath('.\SubFunctions');
[FileName,FilePath]=uigetfile(defaultpath);
FileInfo=strsplit(FileName(1:end-4),'-');
myFileName = ['StackedPDI-',strjoin(FileInfo(2:end-1),'-'),'.mat'];

%% save P
myFile=matfile([FilePath,FileName]);
P=myFile.P;

%% identify frames to be process
prompt = {'Start Reapt', 'Number of Repeats'};
name = 'Power Doppler Data Post Processing';
defaultvalue={'1','1'};
numinput=inputdlg(prompt,name, 1, defaultvalue);
nf = str2num(numinput{2})-str2num(numinput{1})+1;

%% noise equalization and form PDI
for i = 1: nf
    iFileInfo=FileInfo;         
    iFileInfo{8}=num2str(i);
    iFileName = [strjoin(iFileInfo,'-'),'.mat'];
    Y = load([FilePath,iFileName]);
    sIQ = Y.sIQ;
    PDIm(:,:,i) = squeeze(mean(abs(sIQ./sNoiseMedNorm).^2,3)); 
    PDImdb(:,:,i)= log10(PDIm(:,:,i)); %./max(max(PDIm(:,:,i))));
%     PDImdb(:,:,i) = (PDImdb0(:,:,i)-min(min(PDImdb0(:,:,i))))/(max(max(PDImdb0(:,:,i)))-min(min(PDImdb0(:,:,i))));% normalized (0,1)
end
%% save
save([FilePath,myFileName],'P','PDImdb');

%%
PDImdb0 = (PDImdb-min(min(min(PDImdb))))/(max(max(max(PDImdb)))-min(min(min(PDImdb))));
figure;[PDIroi, rect] = imcrop(PDImdb0(:,:,end));
rect = floor(rect);

clear Roi Roiv;
for i = 1: nf
    Roi(:,:) = PDImdb(rect(2):(rect(2)+rect(4)),rect(1):(rect(1)+rect(3)),i);
    Roiv(:,i) = Roi(:);
end
mRoiv = mean(Roiv,1);
figure(1);
% wv1 = 700:1:900;
% mRoiv1 = interp1(wv,mRoiv,wv1,'linear');
plot(mRoiv','.-'); %ylim([0.2,0.5]);%wv1,mRoiv1,'-');
%%

fUSmovie(fliplr(PDImdb0),15); % 15min

