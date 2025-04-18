



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