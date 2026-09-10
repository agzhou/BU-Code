%% linear array
% r = reshapeRcvData(RcvData, P);
% temp
P.numSubFrames = P.numFramesPerBuffer;
IQ = LA_DAS(RcvData, P);



I = squeeze(abs(sum(IQ, 3)));
figure; imagesc(I(:, :, 1))
figure; imagesc(I(:, :, 2))