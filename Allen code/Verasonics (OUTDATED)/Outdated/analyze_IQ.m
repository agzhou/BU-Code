% clear
% P_path = 'D:\Allen\Simulation Results\RC15gV\datatest\params.mat';
% RcvData_path = 'D:\Allen\Simulation Results\RC15gV\datatest\RC15gV_RcvData_31_Oct_2024_21_43_47_681';
% 
% load(P_path);
% load(RcvData_path)
%%
% figure; imagesc(squeeze(RcvData(:, :, 1)))
% figure; imagesc(squeeze(RcvData(:, :, 2)))

% numSupFramesToUse = P.numSupFrames;
% [RcvDataStacked, P_new] = stackSuperFrames(RcvData, P, numSupFramesToUse);
% r = reshapeRcvData(RcvDataStacked, P_new);
% r_old = r;
% r = r(:, :, 1:2);
% P_new.numSubFrames = 2;
%%
% IQ = RcvData2IQ3D(RcvData, P);

% IQ = RCA_DAS(RcvData, P, P.wl/2);
IQ = RCA_DAS(RcvData, P, P.wl/P.Receive(1).samplesPerWave/2);
% IQ = RcvData2IQ3D(r, P_new);


% IQ_coherent_sum = squeeze(sum(IQ, 4));
% I_coherent_sum = abs(IQ_coherent_sum); % intensity
% figure; imagesc(squeeze(I_coherent_sum(40, :, :))')
% title('xy')
% xlabel('x pixels')
% ylabel('z pixels')
% figure; imagesc(squeeze(I_coherent_sum(:, 40, :))')
% title('yz')
% xlabel('x pixels')
% ylabel('z pixels')

%%
IQ_coherent_sum = squeeze(sum(IQ, 4));
I_coherent_sum = abs(IQ_coherent_sum); % intensity
figure; imagesc(squeeze(I_coherent_sum(40, :, :, 1))')
title('xz')
xlabel('x pixels')
ylabel('z pixels')

figure; imagesc(squeeze(I_coherent_sum(:, 40, :, 1))')
title('yz')
xlabel('x pixels')
ylabel('z pixels')
%% RCA test
% scatter3(P.Media.MP(:, 1), P.Media.MP(:, 2), P.Media.MP(:, 3))
[r, P_new] = stackSuperFrames(RcvData, P, 1);
% figure; imagesc(RcvData(:, :, 1));
% r = squeeze(RcvData(:, :, 1));
% P.numSubFrames = P.numSubFrames * P.numSupFrames;
reshapedRcvData = reshapeRcvData(r, P_new);
clear r; clear RcvData;
IQ = RcvData2IQ3D(reshapedRcvData, P_new);

%% RCA superframe stacking with separate recon
r = squeeze(RcvData(:, :, 1));
P_new = P;
P_new.numSubFrames = P.numSubFrames * P.numSupFrames;
P_new.numSupFrames = 1;
reshapedRcvData = reshapeRcvData(r, P);
clear r;
IQ_sf1 = RcvData2IQ3D(reshapedRcvData, P);
clear reshapedRcvData;
temp_size = size(IQ_sf1);
IQ = zeros(temp_size(1), temp_size(2), temp_size(3), temp_size(4), temp_size(5) * P.numSupFrames);
IQ(:, :, :, :, 1:P.numSubFrames) = IQ_sf1;
clear IQ_sf1;

for nsupf = 2:P.numSupFrames
    IQ(:, :, :, :, (nsupf - 1) * P.numSubFrames + 1 : nsupf * P.numSubFrames) = RcvData2IQ3D(reshapeRcvData(squeeze(RcvData(:, :, nsupf)), P), P);

end
%% 12/3/24 single superframe saving test
P_new = P;
P_new.numSupFrames = 1;
reshapedRcvData = reshapeRcvData(RcvData, P_new);
%%
reshapedRcvData = reshapeRcvData(r, P);

% reshapedRcvData = reshapeRcvData(RcvData, P);
IQ = RcvData2IQ3D(reshapedRcvData, P);

% rRD2f = reshapedRcvData(:, :, 1:2);
% P2f = P; P2f.numSubFrames = 2;
% IQ = RcvData2IQ3D(rRD2f, P2f);


% IQ = RcvData2IQ3D(RcvData, P);

% for fn = 1:P.numSubFrames
%     plotRecon(IQ, P, fn)
% end
%% Linear array recon
[r, P_new] = stackSuperFrames(RcvData, P, P.numSupFrames);
% figure; imagesc(RcvData(:, :, 1));
% r = squeeze(RcvData(:, :, 1));
% P.numSubFrames = P.numSubFrames * P.numSupFrames;
reshapedRcvData = reshapeRcvData(r, P_new);

reshapedRcvData_full = reshapedRcvData;
reshapedRcvData = reshapedRcvData(:, :, 1:2);
P_new.numSubFrames = 2;

% clear r; clear RcvData;
IQ = RcvData2IQ2D(reshapedRcvData, P_new);

IQ_coherent_sum = squeeze(sum(IQ, 3));
I_coherent_sum = abs(IQ_coherent_sum); % intensity
figure; imagesc(I_coherent_sum(:, :, 1))
%% for use with Verasonics
addpath 'J:\My Drive\Verasonics files\Vantage-4.9.2-2308102000\Allen code'
IQ_coherent_sum = squeeze(IQ);
% IQ_coherent_sum = squeeze(sum(IQ, 4));
I_coherent_sum = abs(IQ_coherent_sum); % intensity

volumeViewer(squeeze(I_coherent_sum(:, :, :, 1)), scaleFactors = [1, 1, 1])

% xz plane
figure; imagesc(squeeze(I_coherent_sum(40, :, :, 1))')
title('Intensity - xz plane')
xlabel('x pixels')
ylabel('z pixels')
% 
% figure; imagesc(squeeze(sum(I_coherent_sum(:, :, :, 1), 2) ./ P.numElements)')
% title('Intensity - yz plane averaged over x')
% xlabel('y pixels')
% ylabel('z pixels')
%%
sv_threshold_lower = 8;
sv_threshold_upper = 100;
% sv_threshold_lower = 1;
% sv_threshold_upper = 48;
% [IQ_f, SVs] = svd_declutter_2D_acrossframes(IQ_coherent_sum, sv_threshold_lower, sv_threshold_upper);
[IQ_f, SVs, V_sort] = svd_declutter_2D_acrossframes_test_test(IQ_coherent_sum, sv_threshold_lower, sv_threshold_upper);
%%
figure; plot(SVs, '-o')
title('Singular values')
xlabel('SV #')
ylabel('SV magnitude')

figure; plot(log10(SVs), '-o')
title('Singular values log plot')
xlabel('SV #')
ylabel('log10(SV magnitude)')

volumeViewer(squeeze(abs(IQ_f(:, :, :, 1))), scaleFactors = [1, 1, 1])

%%
% I_f = abs(IQ_f);
%
% for fn = 1:P.numSubFrames
% for fn = 10
%     figure; imagesc(squeeze(I_f(40, :, :, fn)))
%     figure; imagesc(squeeze(I_f(:, 40, :, fn)))
% end

%  volumeViewer(squeeze(I_f(:, :, :, 1)), scaleFactors = [1, 1, 1])
%
% taustep = P.SeqControl(1).argument ./ 1e6 .* 1e3 .* P.na .* 2 .* (0:P.numSubFrames - 1); % in ms
taustep = P.SeqControl(1).argument ./ 1e6 .* 1e3 .* P_new.na .* 2 .* (0:P_new.numSubFrames - 1); % in ms

g1 = g1test(IQ_f);

% pt = [40, 44, 63]; % y, x, z
% pt = [44, 40, 63];
pt = [40, 47, 54];

figure; plot(taustep, abs(squeeze(g1(pt(1), pt(2), pt(3), :))), '-o')
xlabel('tau (ms)')
ylabel('|g1|')
title(strcat("|g1| at ", num2str(pt), " pixel")) % (x, y, z)

mag_g1 = abs(g1);
figure; imagesc(squeeze(mag_g1(40, :, :, 3))')
% test = g1(:, pt(2), pt(3), :); test = squeeze(sum(test, 1)) ./ P.numSubFrames;
% figure; plot(taustep, abs(test), '-o')
% xlabel('tau (ms)')
% ylabel('|g1|')
% title(strcat("avg |g1| across y at ", num2str(pt), " pixel")) % (x, y, z)

% figure; plot(taustep, abs(squeeze(I_coherent_sum(pt, :))), '-o')
% xlabel('tau (ms)')
% ylabel('I')
% title('I at (40, 41, 47) pixel') % (x, y, z)

%% Unfiltered Power Doppler for CBV comparison
IQ_cs_sq = IQ_coherent_sum.^2;
I_PowerDoppler = abs(squeeze(sum(IQ_cs_sq, 4)));
figure; imagesc(squeeze(sum(I_PowerDoppler, 1) ./ P.numElements)')
title('Unfiltered Power Doppler - xz plane, average over y')
xlabel('x pixels')
ylabel('z pixels')

figure; imagesc(squeeze(sum(I_PowerDoppler, 2) ./ P.numElements)')
title('Unfiltered Power Doppler - yz plane, average over x')
xlabel('x pixels')
ylabel('z pixels')

volumeViewer(I_PowerDoppler, scaleFactors = [1, 1, 1])

%% Filtered Power Doppler for CBV comparison
IQ_cs_f_sq = IQ_f.^2;
I_f_PowerDoppler = abs(squeeze(sum(IQ_cs_f_sq, 4)));
figure; imagesc(squeeze(I_f_PowerDoppler(40, :, :))')
title('Filtered Power Doppler - xz plane')
xlabel('x pixels')
ylabel('z pixels')

figure; imagesc(squeeze(I_f_PowerDoppler(:, 40, :))')
title('Filtered Power Doppler - yz plane')
xlabel('x pixels')
ylabel('z pixels')

% volumeViewer(I_f_PowerDoppler, scaleFactors = [1, 1, 1])

%%
[CBF, CBV] = g1_to_CBi(g1, taustep, 2, 37, 5); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
%%
figure; imagesc(squeeze(CBF(40, :, :))')
title('CBFi - xz plane')
xlabel('x pixels')
ylabel('z pixels')
figure; imagesc(squeeze(CBF(:, 40, :))')
title('CBFi - yz plane')
xlabel('x pixels')
ylabel('z pixels')

%%
figure; imagesc(squeeze(CBV(40, :, :))')
title('CBVi - xz plane')
xlabel('x pixels')
ylabel('z pixels')
figure; imagesc(squeeze(CBV(:, 40, :))')
title('CBVi - yz plane')
xlabel('x pixels')
ylabel('z pixels')

% figure; imagesc(squeeze(CBV(40, :, 30:end))')
% title('CBVi - xz plane')
% xlabel('x pixels')
% ylabel('z pixels')
% figure; imagesc(squeeze(CBV(:, 40, 30:end))')
% title('CBVi - yz plane')
% xlabel('x pixels')
% ylabel('z pixels')

% Ics_frame_avg = squeeze(sum(I_coherent_sum, 4)) ./ P.numSubFrames;
% figure; imagesc(squeeze(Ics_frame_avg(40, :, :))')

%%
volumeViewer(CBF, scaleFactors = [1, 1, 1])
volumeViewer(CBV, scaleFactors = [1, 1, 1])
%%
I_coherent_sum_frame_avg = sum(I_coherent_sum, 4) ./ P.numSubFrames;
figure; imagesc(squeeze(I_coherent_sum_frame_avg(40, :, :) ./ P.numElements)')
figure; imagesc(squeeze(I_coherent_sum_frame_avg(:, 40, :) ./ P.numElements)')
%% linear array test
reshapedRcvData = reshapeRcvData(RcvData, P);
IQ = RcvData2IQ2D(reshapedRcvData, P);

fn = 1;
plotRecon(IQ, P, fn)

%% RCA no pair angle test
IQ = RcvData2IQ3D_nopair(RcvData, P);
%
% fn = 2;
for fn = 1:P.numSubFrames
    plotRecon(IQ, P, fn)
end

% g1 = svd_proc_