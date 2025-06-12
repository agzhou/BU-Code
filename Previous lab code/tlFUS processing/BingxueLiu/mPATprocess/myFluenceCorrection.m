clear all;
load D:\CODE\Mains\PAT\PATAQParameters_L22_14v.mat;
addpath D:\CODE\DataProcess\BingxueLiu\mPATprocess\SubFunctions;

defaultpath=PATAQInfo.savepath;
[FileName,FilePath]=uigetfile(defaultpath);
fileInfo=strsplit(FileName(1:end-4),'-');
prompt={'Start Repeat', 'Number of Repeats', 'Start Wavelength(nm)', 'End Wavelength(nm)', 'Wavelength Interval(nm)'};
name='File info';
defaultvalue={'1', '17', '690', '850', '10'};
numinput=inputdlg(prompt,name, 1, defaultvalue);
nStart = str2num(numinput{1});
nEnd = str2num(numinput{2});
mWL = str2num(numinput{3}): str2num(numinput{5}): str2num(numinput{4});
nWL = size(mWL,2);
dWL = str2num(numinput{5});
nAver = (nEnd-nStart+1)/nWL;
disp('Loading PA images...')
for iRpt = nStart: nEnd
    iFileInfo=fileInfo(1:end-1);
    iFileInfo{8}=num2str(iRpt); % or7
    iFileName=[strjoin(iFileInfo,'-'),'.mat'];
    Y = load ([FilePath, iFileName]);
    Xm(:,:,iRpt) = jitCor(Y.IQ); %% for jitter correction for laser master data;    
end
disp('Data loaded!')

%% load fluence map
disp('Load fluence map...')
load 'D:\CODE\DataProcess\BingxueLiu\mPATprocess\fluence_map'; % load fluence map for fluence correction
[nx0, ny0] = size(ratiophi_map);
[nx, ny] = size(Xm(:,:,1));
[xcoor0, ycoor0] = meshgrid(1:1/ny0:2-1/ny0, 1:1/nx0:2-1/nx0);
[xcoor, ycoor] = meshgrid(1:1/ny:2-1/ny, 1:1/nx:2-1/nx);
ratiophi_map = interp2(xcoor0, ycoor0, ratiophi_map, xcoor, ycoor,'same'); % interpolate fluence map;

%% fluence equalization
filename = '881549_02.txt';
[Pmeter,delimiterOut]=importdata([FilePath,filename]);
Amp = Pmeter.data(:,2); % energy per pulse recorded by starlab
for i = 1: nWL
    Ampm(i) = sum(Amp(nAver*(i-1)+1:nAver*i))/nAver; % averaged pulse energy per wl
%     Xm(:,:,i) = sum(X(:,:,nAver*(i-1)+1:nAver*i),3)/nAver; % averaged PA image per wl
    Xmn(:,:,i) = Xm(:,:,i)./Ampm(i); % fluence equalization
    Xmnexp(:,:,i) = Xmn(:,:,i)./ratiophi_map;  % fluence correction
end
Ampm = Ampm./min(Ampm);

myFileName = 'FluenceCorrectionData.mat';
save ([FilePath, myFileName],'Xmn','Xmnexp','mWL');

%% make plots

% Xn1 = (X(:,:,1)-min(min(X(:,:,1))))/(max(max(X(:,:,1)))-min(min(X(:,:,1))));
% figure; [roi, rect] = imcrop(Xn1);
% rect = floor(rect);
% Xm = zeros(160,256,21);
% clear Roi Roiv;
% 
% for i = 1: 21
%   
%     Roi(:,:) = Xmnexp(rect(2):(rect(2)+rect(4)),rect(1):(rect(1)+rect(3)),i);
%     Roiv(:,i) = Roi(:);
% end
% mRoiv = mean(Roiv,1);
% figure(2);hold on;
% wv = 700:10:900;
% wv1 = 700:1:900;
% mRoiv1 = interp1(wv,mRoiv,wv1,'linear');
% plot(wv,mRoiv,'.',wv1,mRoiv1,'-');
% 
% MakeAvimPAT(Xmn, [700:10:900], FilePath);



