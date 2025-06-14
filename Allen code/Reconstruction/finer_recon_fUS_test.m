
%% Load a test file

if ~exist('RcvData', 'var')
    load('E:\Allen BME-BOAS-27 Data Backup\AZ01 fUS RCA\06-05-2025 awake manual right whisker stim\run 1 5 trials wooden stick right whiskers 11 angles -5 to 5 deg 2500 Hz\RF-5-11-2500-180-1-2.mat')
end
%% Unstack the RcvData and get the updated parameter structure for unstacked data
[P_unstacked] = updateParams_unstackedFrames_selfrecon(P);

r = unstackFrames(RcvData, P);
% clearvars RcvData;

%%
IQ_hr = RCA_DAS(r, P_unstacked, P_unstacked.wl/2, P_unstacked.wl/2, P_unstacked.wl/2);