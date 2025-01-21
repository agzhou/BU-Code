
%%

% IQ_coherent_sum = squeeze(sum(IQ, 3));
IQ_coherent_sum = squeeze(sum(IQstack, 3));
I_coherent_sum = abs(IQ_coherent_sum); % intensity
figure; imagesc(I_coherent_sum(:, :, 100))

%% SVD proc part 1
[PP, EVs, V_sort] = getSVs1D(IQ_coherent_sum);
%% SVD proc part 2
sv_threshold_lower = 11;
sv_threshold_upper = 100;
[IQ_f] = applySVs1D(IQ_coherent_sum, PP, EVs, V_sort, sv_threshold_lower, sv_threshold_upper);

%% Plot filtered data
abs_IQ_f = abs(IQ_f);
figure; imagesc(abs_IQ_f(:, :, 1))
figure; imagesc(abs_IQ_f(:, :, 100))

%%
figure; plot(EVs, '-o')
title('Singular values')
xlabel('SV #')
ylabel('SV magnitude')

figure; plot(log10(EVs), '-o')
title('Singular values log plot')
xlabel('SV #')
ylabel('log10(SV magnitude)')

%%
taustep = P.SeqControl(1).argument ./ 1e6 .* 1e3 .* P_new.na .* (0:P_new.numSubFrames - 1); % in ms

g1 = g1T_1D(IQ_f);

% pt = [40, 44, 63]; % y, x, z
% pt = [44, 40, 63];
pt = [97, 66];

figure; plot(taustep, abs(squeeze(g1(pt(1), pt(2), :))), '-o')
xlabel('tau (ms)')
ylabel('|g1|')
title(strcat("|g1| at ", num2str(pt), " pixel")) % (z, x)

mag_g1 = abs(g1);

%% Filtered Power Doppler for CBV comparison
IQ_cs_f_sq = IQ_f.^2;
I_f_PowerDoppler = abs(squeeze(sum(IQ_cs_f_sq, 3)));
figure; imagesc(I_f_PowerDoppler)
title('Filtered Power Doppler - xz plane')
xlabel('x pixels')
ylabel('z pixels')
%% Jianbo Power Doppler
[PDI]=sIQ2PDI(IQ_f);

figure; imagesc(PDI(:, :, 1))
title('Power Doppler - positive frequencies')
xlabel('x pixels')
ylabel('z pixels')
figure; imagesc(PDI(:, :, 2))
title('Power Doppler - negative frequencies')
xlabel('x pixels')
ylabel('z pixels')
figure; imagesc(PDI(:, :, 3))
title('Power Doppler - all frequencies')
xlabel('x pixels')
ylabel('z pixels')

%% same thing but rescaled
dr_scaling = 0.3; % dynamic range scaling factor

% figure; imagesc(PDI(:, :, 1).^dr_scaling)
% title('Power Doppler - positive frequencies')
% xlabel('x pixels')
% ylabel('z pixels')
%
figure; imagesc(PDI(:, :, 2).^dr_scaling)
title('Power Doppler - negative frequencies')
xlabel('x pixels')
ylabel('z pixels')
%
figure; imagesc(PDI(:, :, 3).^dr_scaling)
title('Power Doppler - all frequencies')
xlabel('x pixels')
ylabel('z pixels')
%%
[CBF, CBV] = g1_to_CBi(g1, taustep, 2, 15, 4); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
%% CBFi
figure; imagesc(CBF)
title('CBFi - xz plane')
xlabel('x pixels')
ylabel('z pixels')

%% CBVi
figure; imagesc(CBV)
title('CBVi - xz plane')
xlabel('x pixels')
ylabel('z pixels')
%%
%%
dr_scaling = 0.3; % dynamic range scaling factor

figure; imagesc(CBV.^dr_scaling)

title(strcat("CBVi - xz plane with exponential scaling factor = ", num2str(dr_scaling)))
xlabel('x pixels')
ylabel('z pixels')
colormap gray