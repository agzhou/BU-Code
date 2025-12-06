uiopen;
wvlength = P.vSound/P.frequency/1e3; %mm

cd('C:\Users\BOAS-US\Documents\k-wave-toolbox-version-1.3\k-Wave');
%% figure for PAT
% Img = abs(img);
zaxis = (P.StartDepth:P.EndDepth)*wvlength;
xaxis = (-P.nCh/2:P.nCh/2-1)*P.pitch;
figure; imagesc(xaxis,zaxis,Img); colormap(hot(256));axis image;
xlabel('x-position (mm)'); ylabel('z-position (mm)');

figure; imagesc(xaxis,zaxis,sqrt(Img)); colormap(hot(256));axis image;
xlabel('x-position (mm)'); ylabel('z-position (mm)');

% figure; imagesc(xaxis,zaxis(150:230),Img(150:230,:)); colormap(hot(256));axis image;
% xlabel('x-position (mm)'); ylabel('z-position (mm)');colorbar;
% 
% figure; imagesc(xaxis,zaxis(150:230),sqrt(Img(150:230,:))); colormap(hot(256));axis image;
% xlabel('x-position (mm)'); ylabel('z-position (mm)'); colorbar;


%% calculate the frequency spectrum
Fs = 31.25e6; % sampling freq
func = double(RFRAW(1:P.actZsamples,1:128));
[f, func_as, func_ps] = spect(func, Fs);
func_asbar = sum(func_as,2)/P.nCh;
figure; plot(f/1e6,func_asbar);xlabel('Frequency/MHz');
figure; plot(f/1e6,10*log(func_asbar/max(func_asbar)));xlabel('Frequency/MHz');


%% figure for fus
zaxis = (P.startDepth:0.25:P.endDepth+1-0.25)*wvlength;
xaxis = (-P.nCh/2:0.25:P.nCh/2-0.25)*P.pitch;
figure; imagesc(xaxis,zaxis,abs(IQ(:,:,1))); colormap(gray(256));axis image;
xlabel('x-position (mm)'); ylabel('z-position (mm)');
