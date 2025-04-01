
%%
figure; imagesc(squeeze(RcvData(:, :, 1)))
figure; imagesc(squeeze(RcvData(:, :, 2)))

%%
test = reshapeRcvData(RcvData, P);
size(test)

figure; imagesc(squeeze(test(:, :, 1, 1)))
figure; imagesc(squeeze(test(:, :, 1, 2)))

figure; imagesc(squeeze(test(:, :, 4, 1)))
figure; imagesc(squeeze(test(:, :, 4, 2)))
%%
reshapedRcvData = reshapeRcvData(RcvData, P);
for nsupf = 1
    
    
    IQ = RcvData2IQ3D(squeeze(reshapedRcvData(:, :, :, 1)), P);
    nsubf = 1;
    plotRecon(IQ, P, nsubf)
end