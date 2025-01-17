
%%
filename1 = 'sf';
tic
savefast([savepath, filename1], 'RcvData')
toc
%% This one is too slow
filename2 = 'v7.3noc';
tic
save([savepath, filename2], "RcvData", "-v7.3", "-nocompression")
toc
%%
filename3 = 'v6';
tic
save([savepath, filename3], "RcvData", "-v6")
toc

%%
savepath = 'G:\Allen\Data\01-15-2025 test\L22-14v\run 5\';
filename4 = 'v6t';
a = randi(4, 2);
tic
% save([savepath, filename4], "a", "-v6")

save([savepath, filename4], "a", "-v6", "-append")
toc