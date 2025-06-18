%% Look at the SVD and the SV magnitudes across superframes

%% load params and stuff
IQpath = uigetdir('D:\Allen\Data\', 'Select the IQ data path');
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

savepath = uigetdir('D:\Allen\Data\', 'Select the save path');
savepath = [savepath, '\'];

%% Define some parameters (add this to a prompt later)
% 
% parameterPrompt = {'Start file number', 'End file number', 'SVD lower bound', 'SVD upper bound', 'Tau 1 index for CBFspeed', 'Tau 2 index for CBFspeed', 'Tau 1 index for CBV'};
% parameterDefaults = {'1', '', '15', '180', '2', '6', '2'};
% parameterUserInput = inputdlg(parameterPrompt, 'Input Parameters', 1, parameterDefaults);
% 
% % define # of files manually for now
% % str2double(parameterUserInput{});
% startFile = str2double(parameterUserInput{1});
% endFile = str2double(parameterUserInput{2});
% numFiles = endFile - startFile + 1;
% sv_threshold_lower = str2double(parameterUserInput{3});
% sv_threshold_upper = str2double(parameterUserInput{4});

%% Main loop
for filenum = startFile:endFile
% for filenum = 11:endFile
    tic
    load([IQpath, IQfilenameStructure, num2str(filenum)])
    
    IQ = squeeze(IData + 1i .* QData);
    clearvars IData QData
    
    % SVD  
    [PP, EVs, V_sort] = getSVs1D(IQ);
    disp('SVs decomposed')
    % [IQf] = applySVs2D(IQ, PP, EVs, V_sort, sv_threshold_lower, sv_threshold_upper);
    % disp('SVD filtered images put together')

    save([savepath, 'SVD_data-', num2str(filenum), '.mat'], 'EVs', 'V_sort', '-v7.3', '-nocompression');
    disp("SVD result for file " + num2str(filenum) + " saved" )
    toc

end

%% Look at each superframe's singular value distribution

figure; hold on
for filenum = startFile:endFile
% for filenum = 1:10
    load([savepath, 'SVD_data-', num2str(filenum)]) % Load the saved SVD mat files

    % figure; histogram(EVs)
    plot(log10(EVs), 'o-')
%     plot(EVs, 'o-')

%     save([savepath, 'tlfUSdata-', num2str(filenum), '.mat'], 'CBFsi', 'CBVi', '-v7.3', '-nocompression');
% %     save([savepath, 'tlfUSdatatest-', num2str(filenum), '.mat'], 'CBFsi', 'CBVi', '-v7.3', '-nocompression');
%     disp("tl-fUS result for file " + num2str(filenum) + " saved" )

end
hold off
xlabel('Singular value index')
ylabel('log10(Magnitude)')

%% Get the IQf for each superframe
sv_threshold_lower = 20;
sv_threshold_upper = 180;

% for filenum = startFile:endFile
for filenum = 3:5
    % load([savepath, 'SVD_data-', num2str(filenum)]) % Load the saved SVD mat files

    tic
    load([IQpath, IQfilenameStructure, num2str(filenum)])
    
    IQ = squeeze(IData + 1i .* QData);
    clearvars IData QData
    
    % SVD decluttering
    [xp, yp, zp, nf] = size(IQ);
    [PP, EVs, V_sort] = getSVs2D(IQ);
    disp('SVs decomposed')
    [IQf] = applySVs2D(IQ, PP, EVs, V_sort, sv_threshold_lower, sv_threshold_upper);
    disp('IQ filtered')

    save([savepath, 'IQf-', num2str(filenum), '.mat'], 'IQf', '-v7.3', '-nocompression');
    toc
end

%% Inspect the dynamic range of the non-contrast-matched IQf across superframes

% figure; hold on
% xmax = -Inf;

% for filenum = startFile:endFile
for filenum = 1:5
    load([savepath, 'IQf-', num2str(filenum)]) % Load the saved IQf mat files

    temp_absIQf = abs(IQf);
    % temp_absIQf = abs(IQf) ./ max(abs(IQf), [], 'all');

    % figure; histogram(temp_absIQf(:), 'Normalization', 'pdf')
    figure; histogram(temp_absIQf(:), 'Normalization', 'count')
    xlabel('|IQf| magnitude')
    % figure; histogram(log10(temp_absIQf(:)))
    % xlabel('log10(|IQf| magnitude)')
    ylabel('Counts')
    title("Superframe " + num2str(filenum))

    % xmax = max(xmax, max(temp_absIQf(:)));

    % plot(log10(EVs), 'o-')
    % plot(EVs, 'o-')

%     save([savepath, 'tlfUSdata-', num2str(filenum), '.mat'], 'CBFsi', 'CBVi', '-v7.3', '-nocompression');
% %     save([savepath, 'tlfUSdatatest-', num2str(filenum), '.mat'], 'CBFsi', 'CBVi', '-v7.3', '-nocompression');
%     disp("tl-fUS result for file " + num2str(filenum) + " saved" )

end
% hold off
%%
ymax = 7e5;
fh_offset = 9;
for filenum = 1:1
    figure(fh_offset + filenum)
    xlim([0, xmax])
    % ylim([0, ymax])

end