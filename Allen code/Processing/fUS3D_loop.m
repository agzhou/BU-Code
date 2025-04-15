%% load params and stuff
% IQpath = 'D:\Allen\Data\04-11-2025 AZ02 fUS RC15gV\run 1 all frames stacked\IQ Data - Verasonics recon\';
IQpath = uigetdir('G:\Allen\Data\', 'Select the IQ data path');
IQpath = [IQpath, '\'];

% Load parameters
% if ~exist('P', 'var')
%     load([IQpath, '..\params.mat'])
% end
% Load acquisition parameters: params.mat
if ~exist('P', 'var')
    % Choose and load the params.mat file (from the acquisition)
    [params_filename, params_pathname, ~] = uigetfile('*.mat', 'Select the params file', [IQpath, '..\params.mat']);
    load([params_pathname, params_filename])
end

% Load Verasonics reconstruction parameters: datapath\PData.mat
if ~exist('PData', 'var')
    load([IQpath, 'PData.mat'])
end

IQfilenameStructure = ['IQ-', num2str(P.maxAngle), '-', num2str(P.na), '-', num2str(P.frameRate), '-', num2str(P.numFramesPerBuffer), '-1-'];

savepath = uigetdir('G:\Allen\Data\', 'Select the save path');
savepath = [savepath, '\'];

%% Main loop
sv_threshold_lower = 40;
sv_threshold_upper = 400;

startFile = 1;
endFile = 148;

taustep = 1/P.frameRate;
% tau = taustep:taustep:(P.numFramesPerBuffer * taustep);
tau = 0:taustep:((P.numFramesPerBuffer - 1) * taustep);
tau_ms = tau .* 1000; % Assuming even time spacing between frames

tau1_index_CBF = 2;
tau2_index_CBF = 6;
tau1_index_CBV = 2;

% for filenum = startFile:endFile
% for filenum = 2:endFile
for filenum = 1
    tic
    load([IQpath, IQfilenameStructure, num2str(filenum)])
    
    IQ = squeeze(IData + 1i .* QData);
    clearvars IData QData
    
    % SVD decluttering
    [xp, yp, zp, nf] = size(IQ);
    
    [PP, EVs, V_sort] = getSVs2D(IQ);
    disp('SVs decomposed')
    [IQf] = applySVs2D(IQ, PP, EVs, V_sort, sv_threshold_lower, sv_threshold_upper);
    disp('SVD filtered images put together')

    clearvars IQ  
    
    g1 = g1test(IQf);
    
    [CBFi, CBVi] = g1_to_CBi(g1, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)

%     savefast([savepath, 'fUSdata-', num2str(filenum), '.mat'], g1, CBFi, CBVi);
    save([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'g1', 'CBFi', 'CBVi', '-v7.3', '-nocompression');
    disp("fUS result for file " + num2str(filenum) + " saved" )
    toc
    
end
savefast([savepath, 'fUS_proc_params.mat'], 'sv_threshold_lower', 'sv_threshold_upper', 'tau', 'tau_ms', 'tau1_index_CBF', 'tau2_index_CBF', 'tau1_index_CBV');

%%
volumeViewer(abs(IQf(:, :, :, 1)))
%%
figure; imagesc(abs(squeeze(max(IQf(:, :, :, 1), [], 1)))')
%%


%% Plot the magnitude of g1 at some point
figure; plot(tau_ms, abs(squeeze(g1(40, 45, 61, :))), '-o')
title('|g1| at 40, 45, 61')
xlabel('Tau [ms]')
ylabel('|g1|')

%%
[CBF, CBV] = g1_to_CBi(g1, tau_ms, 2, 3, 2); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
%%
figure; imagesc(squeeze(CBF(40, :, :))')
title('CBFi - xz plane')
xlabel('x pixels')
ylabel('z pixels')
figure; imagesc(squeeze(CBF(:, 40, :))')
title('CBFi - yz plane')
xlabel('x pixels')
ylabel('z pixels')

%% CBF MIP over the whole dimension
figure; imagesc(squeeze(max(CBF, [], 1))' .^ 1); colormap hot
title('CBFi - xz MIP')
xlabel('x pixels')
ylabel('z pixels')
figure; imagesc(squeeze(max(CBF, [], 2))' .^ 1); colormap hot
title('CBFi - yz MIP')
xlabel('x pixels')
ylabel('z pixels')

% figure; imagesc(squeeze(sum(CBF, 1))' .^ 1); colormap hot
% title('CBFi - sum over y')
% xlabel('x pixels')
% ylabel('z pixels')
%%
figure; imagesc(squeeze(CBV(40, :, :))'); colormap hot
title('CBVi - xz plane')
xlabel('x pixels')
ylabel('z pixels')
figure; imagesc(squeeze(CBV(:, 40, :))'); colormap hot
title('CBVi - yz plane')
xlabel('x pixels')
ylabel('z pixels')
%% CBV MIP over the whole dimension
gcp = 1; % gamma compression power
figure; imagesc(squeeze(max(CBV, [], 1))' .^ gcp); colormap hot; colorbar
title('CBVi - xz MIP')
xlabel('y pixels')
ylabel('z pixels')
figure; imagesc(squeeze(max(CBV, [], 2))' .^ gcp); colormap hot; colorbar
title('CBVi - yz MIP')
xlabel('x pixels')
ylabel('z pixels')

%% 
volumeViewer(CBV .^ gcp)