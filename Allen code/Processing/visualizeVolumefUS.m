% function visualizeVolume(volumeData)
%% Plot a bubble density map and make a rotating 3D view gif
% h = volshow(PDI_allBlocks_avg_rs, 'RenderingStyle', 'MaximumIntensityProjection')
% h = volshow(PDI_allBlocks_avg_rs, 'RenderingStyle', 'GradientOpacity')
% h = volshow(PDI_allBlocks_avg_rs, objectConfig, 'Alphamap', 'linear', 'BackgroundColor')

filename = 'H:\My Drive\Presentations\Stalls summit 2026\fUS volume from 02-10-2026\rotating_view.gif';

v = viewer3d(BackgroundColor="white", ...
    BackgroundGradient="off", ...
    OrientationAxes="on");
hFig = v.Parent;
% startPosition = v.CameraPosition;
% v.CameraPosition = startPosition;
% v.CameraUpVector = [0 -1 -0.5];
% v.CameraUpVector = [1 -1 -1];
v.CameraUpVector = [0 0 -1];
% v.CameraUpVector = [1, 1, -1];
% % startPosition = v.CameraPosition;
v.CameraZoom = 1.5;
% sx = 55;
% sy = sx;
% sz = 28.3;
% sn = sqrt(sx^2 + sy ^2 + sz ^2);
% sx = sx / sn; sy = sy / sn; sz = sz / sn;
% A = [sx 0 0 0; 0 sy 0 0; 0 0 sz 0; 0 0 0 1];
% viewertform = affinetform3d(A);

% d = volshow(abs(PDI_allBlocks_avg_rs) .^ 0.5, 'RenderingStyle', 'MaximumIntensityProjection', Parent=v, Transformation=viewertform);
d = volshow(abs(PDI_allBlocks_avg_rs) .^ 0.7, 'RenderingStyle', 'MaximumIntensityProjection', Parent=v);
% cmap = colormap_ULM;
% d = volshow(day7SM, 'RenderingStyle', 'VolumeRendering', 'Colormap', cmap, Parent=v, Transformation=viewertform);
%
%
% Make a rotating camera view
numberOfFrames = 12;
vec = linspace(0, 2*pi, numberOfFrames)';
sz = size(PDI_allBlocks_avg_rs);
% dist = sqrt(sz(1)^2 + sz(2)^2 + sz(3)^2);
dist = sqrt(sum((startPosition - v.CameraTarget) .^2));
myPosition = v.CameraTarget + ([cos(vec) sin(vec) ones(size(vec))]*dist);

for idx = 1:length(vec)
    % Update the current view
    v.CameraPosition = myPosition(idx, :);
    % v.CameraPosition = startPosition - myPosition(idx, :);
    % Capture the image using the getframe function
    I = getframe(hFig);
    [indI,cm] = rgb2ind(I.cdata,256);
    % Write the frame to the GIF file
    if idx==1
        % Do nothing. The first frame displays only the viewer, not the
        % volume.
    elseif idx == 2
        imwrite(indI,cm,filename,"gif",Loopcount=inf,DelayTime=0)
    else
        imwrite(indI,cm,filename,"gif",WriteMode="append",DelayTime=0)
    end
end

%%
v1f = abs(IQf(:, :, :, 1));
vf = figure;
V = volshow(v1f);

V_old = V;
V.Alphamap(1:100) = 0;          % Change transparency
V.BackgroundColor = [1, 1, 1];  % Make background white
V.ScaleFactors(3) = size(v1f, 1) / size(v1f, 3) * 8.8 / 8; % scale with # pixels and region size
V.CameraPosition = V.CameraPosition ./ 2;

cv = getframe(vf);     % get the current volume
rgb = frame2im(cv);      % convert the frame to rgb data

%% Make 3D video of the filtered IQ data
savepath = 'D:\Allen\Data\01-29-2025 AZ001 ULM RC15gV\run 1 left eye\Processed Data\';

datapath = 'D:\Allen\Data\01-29-2025 AZ001 ULM RC15gV\run 1 left eye\';
numFiles = 96;

vo = VideoWriter([savepath, 'IQf_10_150_buffer1']);

vo.Quality = 100;
vo.FrameRate = 30;
open(vo);

vf = figure;
findfigs

for f = 1:size(IQf, 4)
%     for f = 1:100
    v1f = abs(IQf(:, :, :, f));
    
    V = volshow(v1f);
    
%     V_old = V;
    V.Alphamap(1:100) = 0;          % Change transparency
    V.BackgroundColor = [1, 1, 1];  % Make background white
    V.ScaleFactors(3) = size(v1f, 1) / size(v1f, 3) * 8.8 / 8; % scale with # pixels and region size
    V.CameraPosition = V.CameraPosition ./ 2;
    
    cv = getframe(vf);     % get the current volume
    rgb = frame2im(cv);      % convert the frame to rgb data

    writeVideo(vo, rgb);
end
close(vo);
    
%% Make 2D video of the filtered IQ data
savepath = 'D:\Allen\Data\01-29-2025 AZ001 ULM RC15gV\run 1 left eye\Processed Data\';

datapath = 'D:\Allen\Data\01-29-2025 AZ001 ULM RC15gV\run 1 left eye\';
numFiles = 96;

vo = VideoWriter([savepath, 'IQf_10_150_buffer1_2D']);

vo.Quality = 100;
vo.FrameRate = 30;
open(vo);

vf = figure;
findfigs

for f = 1:size(IQf, 4)
%     for f = 1:100
    v1f = abs(IQf(40, :, :, f));
    
    imagesc(squeeze(v1f)')
    
    cv = getframe(vf);     % get the current volume
    rgb = frame2im(cv);      % convert the frame to rgb data

    writeVideo(vo, rgb);
end
close(vo);

%% 2D video of all buffers with loop
addpath('C:\Users\BOAS-US\Documents\Allen\GitHub\BU-Code\Allen code\Processing')

% savepath = 'G:\Allen\Data\01-29-2025 AZ001 ULM\RC15gV\run 1 left eye\Processed Data\';
savepath = 'F:\Allen\Data\01-29-2025 AZ001 ULM\RC15gV\run 1 left eye\FMAS Processed Data\';
mkdir(savepath)

% datapath = 'G:\Allen\Data\01-29-2025 AZ001 ULM\RC15gV\run 1 left eye\IQ Data - Verasonics Recon\';
datapath = 'G:\Allen\Data\01-29-2025 AZ001 ULM\RC15gV\run 1 left eye\IQ Data - Verasonics Recon\';

numFiles = 96;

% add a param loading line
load('G:\Allen\Data\01-29-2025 AZ001 ULM\RC15gV\run 1 left eye\params.mat')

sv_threshold_lower = 20;
sv_threshold_upper = 150;

vo = VideoWriter([savepath, 'FMAS_IQf_20_150_2D_xz_allbuffers']);

vo.Quality = 100;
vo.FrameRate = 100;
open(vo);

vf = figure;
findfigs
colormap gray

% change aspect ratio of the figure
zEnd = 142; zStart = 1; znumpix = 142; % hard coding for now
hwRatio = (P.endDepth - P.startDepth) / (P.Trans.spacing * P.numElements) * (zEnd-zStart + 1) / znumpix; % height to width ratio, idk if the adjustment works (1/12/25 change)

title(strcat("xz plane - ", num2str(P.na), " angles from -", num2str(P.maxAngle), " to ", num2str(P.maxAngle), " deg"))
xlabel('x pixels')
ylabel('z pixels')
colorbar
vf.Position(4) = ceil(vf.Position(3) * hwRatio);

filename_structure = ['IQ-', num2str(P.maxAngle), '-', num2str(P.na), '-', num2str(P.frameRate), '-', num2str(P.numFramesPerBuffer), '-1-'];

% for filenum = 1:numFiles
for filenum = 40:40
    load([datapath, filename_structure, num2str(filenum)])
    disp(strcat("File ", num2str(filenum), " loaded."))

    IQ = squeeze(IData + 1i .* QData);   % Combine I and Q, which are saved separately. It's easier to save the big reconstructed data with savefast, which doesn't support complex values.
    clear IData QData
    
    [PP, EVs, V_sort] = getSVs2D(IQ);

    [IQf] = applySVs2D(IQ, PP, EVs, V_sort, sv_threshold_lower, sv_threshold_upper);
    
    for f = 1:size(IQf, 4)
    %     for f = 1:100
        
        v1f = abs(IQf(40, :, :, f));
    
        imagesc(squeeze(v1f)')
        
        cv = getframe(vf);     % get the current volume
        rgb = frame2im(cv);      % convert the frame to rgb data
    
        writeVideo(vo, rgb);
    end

%     clear IQ IQf PP EVs V_sort
end
close(vo);
















%% with loop

savepath = 'D:\Allen\Data\01-29-2025 AZ001 ULM RC15gV\run 1 left eye\Processed Data\';
vo = VideoWriter([savepath, 'IQf_10_150']);

vo.Quality = 100;
vo.FrameRate = 10;
open(vo);

vf = figure;
findfigs

filename_structure = ['IQ-', num2str(P.maxAngle), '-', num2str(P.na), '-', num2str(P.frameRate), '-', num2str(P.numFramesPerBuffer), '-1-'];


for filenum = 1:numFiles
load([datapath, filename_structure, num2str(filenum)])
disp(strcat('File ', num2str(filenum), ' loaded.'))

% ADD SVD STEP HERE %%%%%%%%%%

for f = 1:size(IQf, 4)
%     for f = 1:100
    v1f = abs(IQf(:, :, :, f));
    
    V = volshow(v1f);
    
%     V_old = V;
    V.Alphamap(1:100) = 0;          % Change transparency
    V.BackgroundColor = [1, 1, 1];  % Make background white
    V.ScaleFactors(3) = size(v1f, 1) / size(v1f, 3) * 8.8 / 8; % scale with # pixels and region size
    V.CameraPosition = V.CameraPosition ./ 2;
    
    cv = getframe(vf);     % get the current volume
    rgb = frame2im(cv);      % convert the frame to rgb data

    writeVideo(vo, rgb);
end
end
close(vo);