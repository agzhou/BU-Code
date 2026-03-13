% Dataset examples. 
% NOTE: THIS FILE NEEDS MATLAB 2019B OR LATER TO RUN CORRECTLY

%% Load Datasets
clearvars
close all

% if(isMATLABReleaseOlderThan("R2019b"))
%     error("This example needs matlab release 2019b or later to work correctly")
% end

varname1 = 'images';
varname2 = 'pLarge';

% Load abodminal data
% imgAbd = load("Abdominal_28-Nov-2023.mat",varname1).(varname1);
% pAbd = load("Abdominal_28-Nov-2023.mat",varname2).(varname2);

% Load abodminal data updated
imgAbd = load("Abdominal_01-Dec-2023.mat",varname1).(varname1);
pAbd = load("Abdominal_01-Dec-2023.mat",varname2).(varname2);

% Load Calf data using Boas lab transducer
% imgCalfBoas = load("CalfData_Boas_27-Nov-2023.mat",varname1).(varname1);
% pCalfBoas = load("CalfData_Boas_27-Nov-2023.mat",varname2).(varname2);


% Load Calf data using Boas lab transducer updated
imgCalfBoas = load("CalfData_Boas_01-Dec-2023.mat",varname1).(varname1);
pCalfBoas = load("CalfData_Boas_01-Dec-2023.mat",varname2).(varname2);

% Load Calf data using GE9LD (our xducer)
imgCalf = load("CalfData_GE9LD_27-Nov-2023.mat",varname1).(varname1);
pCalf = load("CalfData_GE9LD_27-Nov-2023.mat",varname2).(varname2);

% Load Resolution Target Data
imgRes = load("ResolutionTargetsData.mat",varname1).(varname1);
pRes = load("ResolutionTargetsData.mat",varname2).(varname2);

%% Your code here
% Images are stored in 1x10 structures where each structure includes the
% name of the method and the raw unscaled data for the image. The
% corresponding 'p' structure contains the associated parameters, i.e. the
% x coordinates and z coordinates, etc.

imB=imgAbd(1).data;
imJ=imgAbd(10).data;
imB=imB(435:955,345:680);
imJ=imJ(435:955,345:680);

% imB=imgCalfBoas(1).data;
% imJ=imgCalfBoas(10).data;
% imB=imB(365:870,290:480);
% imJ=imJ(365:870,290:480);
 
% imB=imgCalf(1).data;
% imJ=imgCalf(10).data;

imB=abs(imB);
imJ=abs(imJ);

imB=imB/mean(imB(:));
imJ=imJ/mean(imJ(:));

imB=imB.^0.25;

mB=mean(imB(:));
kB=var(imB(:))/mB^2;

% x0=0.05;
% F= @(x) [mB-mean(imJ(:).^x)];
% x=fsolve(F,x0);

x0=0.1;
F= @(x) [kB-var(imJ(:).^x)/mean(imJ(:).^x)^2];
x=fsolve(F,x0);

imJ=imJ.^x;

figure
subplot(1,2,1);imagesc(imB);
subplot(1,2,2);imagesc(imJ);
colormap(gray)

%% Comparison of all modalities
clear im

% abdomen
% for k=1:10
%    im(:,:,k)=abs(imgAbd(k).data(435:955,345:680));
%    im(:,:,k)=im(:,:,k)/mean(im(:,:,k),'all');
% end

% calf
% for k=1:10  
%    im(:,:,k)=abs(imgCalfBoas(k).data(340:820,220:485));
%    im(:,:,k)=im(:,:,k)/mean(im(:,:,k),'all');
% end

% calf zoom
for k=1:10
   im(:,:,k)=abs(imgCalfBoas(k).data(520:660,250:340));
%    im(:,:,k)=im(:,:,k)/mean(im(:,:,k),'all');
end

im(:,:,1)=im(:,:,1).^0.5;
mB=mean(im(:,:,1),'all');
kB=var(im(:,:,1),[],'all')/mB^2;
g(1)=1;

% match means
% for k=2:10
%     g0=0.05;
%     F= @(x) [mB-mean(im(:,:,k).^x,'all')];
%     g(k)=fsolve(F,g0);
% end

% match contrast
for k=2:10
    g0=0.05;
    F= @(x) [kB-var(im(:,:,k).^x,[],'all')/mean(im(:,:,k).^x,'all')^2];
    g(k)=fsolve(F,g0);
end

figure
for k=1:10
    subplot(2,5,k);imagesc(im(:,:,k).^g(k));
end
colormap(gray)


%% Example Plotting code
% plotterFunc1(pAbd,imgAbd)

% If you want to plot a different subsection of the image (i.e. only half
% of the image):
img2 = imgCalfBoas;
for i = 1:length(imgCalfBoas)
    img2(i).data = imgCalfBoas(i).data(365:870,290:480,:);
end
pNew = computeNewGrid(pCalfBoas,[290,480],[365,870]);

plotterFunc1(pNew,img2);

%% Helper functions

% This function scales the images correctly so there is no whitespace
% between them
function plotterFunc1(pLarge,images)

% Define the dimensions of each image in terms of your x and z coordinates.
img_width = range(pLarge.xCoord*1e3);  % Replace with actual width
img_height = range(pLarge.zCoord*1e3);  % Replace with actual height

% Define the number of rows and columns
num_rows = 2;
num_columns = ceil(length(images)/2);

% Calculate the total figure width and height based on the number of rows and columns
scale = 5;
fig_width = img_width * num_columns * scale;
fig_height = img_height * num_rows * scale;

% Create the figure with the calculated dimensions
fig = figure('Position', [50, 50, fig_width, fig_height]);  % [left, bottom, width, height]

t = tiledlayout(num_rows,num_columns,'TileSpacing','none','Padding','none');

for i = 1:length(images)
    nexttile(i,[1,1])
    if (images(i).name == "SLSC") || (images(i).name == "SLAC")
        chk = images(i).data(:,:,1);
%         chk = chk + abs(min(chk(:)));
        colormap gray
        imagesc(pLarge.xCoord*1e3,pLarge.zCoord*1e3,chk)
    elseif (images(i).name == "JCF") || (images(i).name == "UCF")
        % plotLogScaleImage(pLarge.xCoord*1e3,pLarge.zCoord*1e3,images(i).data(:,:,1))
        plotGammaScaleImage(pLarge.xCoord*1e3,pLarge.zCoord*1e3,images(i).data(:,:,1),0.1)
    else
        plotGammaScaleImage(pLarge.xCoord*1e3,pLarge.zCoord*1e3,images(i).data(:,:,1),0.25)
    end
        
    axis image
    % title(images(i).name)
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);

    % Add the subfigure label
    label_char = char(96 + i);  % 'a' is 97 in ASCII
    label_str = ['(' label_char ')'];
    text(min(pLarge.xCoord)*1e3, min(pLarge.zCoord)*1e3, label_str, 'FontSize', 12, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Color', 'w')
end

end