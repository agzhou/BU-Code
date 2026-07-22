AVG_all_coronal2 = [];
savepath = 'F:\OShea Lab\BU1398\slice2'; %saves TIFF here

%load mat files
for k = 1:60
    fname = ['F:\OShea Lab\BU1398\slice2\IQData-15-63-10-1-10-',num2str(k),'.mat'];
    load(fname,'IData')
    AVG_all_coronal2 = cat(5,AVG_all_coronal2,IData);
end

%average files over time and show sample image
test = abs(squeeze(mean(AVG_all_coronal2,5)));
test = imresize(test,[201 258]); %make isotropic
figure(1)
imshow(test,[0 1E8])

%convert to uint8 and save tiff
test = test*(255/max(test(:)));
imwrite(test,fullfile(savepath,'sagittal_test.tiff'))

%%%% images are 12.8mm wide + user defined depth
%%%% before resizing, pixel is 100um in X and 50um in Z
%%%% after resizing, each pixel is 50um in X and Z