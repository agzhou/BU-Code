AVG_all_coronal2 = [];
for k = 1:70
    fname = ['F:\OShea Lab\BU1175\coronal 5\IQData-18-90-10-1-10-',num2str(k),'.mat'];
    load(fname,'IData')
    AVG_all_coronal2 = cat(5,AVG_all_coronal2,IData);
end

test = abs(squeeze(mean(AVG_all_coronal2,5)));
test = imresize(test,[201 258]);
figure(1)
imshow(test,[0 1E8])

test = test*(255/max(test(:)));
imwrite(test,'sagittal_test.tiff')