Asize = 1000;
A = rand(Asize, Asize, Asize);
PSFs = rand(20, 20, 20);
tic
test = normxcorr3(PSFs, A);
toc
%%
B = gpuArray(A);
PSFsg = gpuArray(PSFs);
tic
testgpu = normxcorr3(PSFsg, B);
toc