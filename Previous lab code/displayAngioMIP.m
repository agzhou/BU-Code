clear
close all
clc

%% Jack version %%
dims = {'XY','YZ','ZX'};
ex_params = inputdlg({'Upsample factor:'},'Export Parameters',[1 35],{'2.5'});
[aFile,aFolder] = uigetfile('*.mat');
load(fullfile(aFolder,aFile))
bdm_interpolated = imresize3(bdm_interpolated,2.5);
[sX,sY,sZ] = size(bdm_interpolated);

%normalize
outStack = bdm_interpolated./prctile(bdm_interpolated(:),99);
%outStack = bdm_interpolated;

for i = 1:3
    outStack_1 = permute(outStack,[i mod(i,3)+1 mod(mod(i,3)+1,3)+1]); %cycle permutations
    for j = 1:size(outStack_1,3)
            imwrite((outStack_1(:,:,j)),fullfile(aFolder,[aFile(1:(end-4)),'_',dims{i},'raw2.tif']),'WriteMode','append')
    end
end
