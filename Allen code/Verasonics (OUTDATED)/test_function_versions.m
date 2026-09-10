% clear
% P_path = 'D:\Allen\Simulation Results\RC15gV\datatest\params.mat';
% RcvData_path = 'D:\Allen\Simulation Results\RC15gV\datatest\RC15gV_RcvData_31_Oct_2024_21_43_47_681';
% 
% load(P_path);
% load(RcvData_path)
% %%
% figure; imagesc(squeeze(RcvData(:, :, 1)))
% figure; imagesc(squeeze(RcvData(:, :, 2)))

%% RCA test
% IQ = RcvData2IQ3D(reshapeRcvData(RcvData, P), P);
% reshapedRcvData = reshapeRcvData(RcvData, P);
reshapedRcvData = reshapeRcvData(r, P);
IQ = RcvData2IQ3D_test(reshapedRcvData, P);
fn = 2;
plotRecon(IQ, P, fn)

IQ_coherent_sum = squeeze(sum(IQ, 4));
I_coherent_sum = abs(IQ_coherent_sum); 

% g1 = svd_proc_

% %% linear array test
% reshapedRcvData = reshapeRcvData(RcvData, P);
% IQ = RcvData2IQ2D(reshapedRcvData, P);
% 
% fn = 1;
% plotRecon(IQ, P, fn)
% 
% %% RCA test
% IQ = RcvData2IQ3D_nopair(RcvData, P);
% %%
% % fn = 2;
% for fn = 1:P.numSubFrames
%     plotRecon(IQ, P, fn)
% end
% 
% % g1 = svd_proc_