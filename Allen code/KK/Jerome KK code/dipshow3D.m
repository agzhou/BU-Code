function dipshow3D(window,data)
% displays 2D or 3D matrix with dipshow (scroll with 'n' and 'p')
data_=newim(size(data,2),size(data,1),size(data,3));
data_=mat2im(data);
dipshow(window,data_,'lin')
diptruesize(window,'off')
end

