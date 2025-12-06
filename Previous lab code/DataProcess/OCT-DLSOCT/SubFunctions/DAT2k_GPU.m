% Convert raw spectrum data to Reflectivity RR, GPU-based
function DAT_k = DAT2k_GPU(DAT, intDk)
[nk,Nx,ny] = size(DAT); % Dat: spectrum data, [nk, Nx, ny]
DAT_k=single(zeros(nk,Nx,ny)); % initialize RR
if nk*Nx*ny>2048*60000*2
    if rem(ny,2)==0
        nyPchk=2;
    else
        nyPchk=1;
    end
    nyChk=ny/nyPchk;
else
    nyPchk=ny;
    nyChk=ny/nyPchk;
end
%% transform from lamda to k, lamda-k interpolation, and ifft
k = linspace(1-intDk/2, 1+intDk/2, nk);
lam = single(gpuArray(1./fliplr(k)));
X=fliplr(lam);
XQ=fliplr(linspace(min(lam),max(lam),length(lam)));
DatGPU=gpuArray(DAT);
DC=flip(single(mean(DatGPU(:,:),2)));
DatGPU=flip(DatGPU,1);
%% lamda to k space
for iChk=1:nyChk
    V = DatGPU(:,:,(iChk-1)*nyPchk+1:iChk*nyPchk) - DC; % substract the reference signal, Subtract mean
    %% lamda to k space
    data_k= flip(interp1(X, V, XQ, 'linear'),1);
    DAT_k(:,:,(iChk-1)*nyPchk+1:iChk*nyPchk)=gather(data_k);
end

