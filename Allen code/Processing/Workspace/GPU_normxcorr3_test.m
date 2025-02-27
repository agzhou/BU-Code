A = rand(500, 500, 500);
tic
test = normxcorr3(PSFs, A);
toc
%%
B = gpuArray(A);
PSFsg = gpuArray(PSFs);
tic
testgpu = normxcorr3(PSFsg, B);
toc