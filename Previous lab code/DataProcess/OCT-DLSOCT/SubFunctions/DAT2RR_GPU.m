% Convert raw spectrum data to Reflectivity RR, GPU-based
% input: 
    % DAT [Dim.nk,RptA_n_P*Dim.nx*Dim.nyRpt,Dim.ny], Nx=RptA_n_P*Dim.nx*Dim.nyRpt
    % intpDk: Lambda to k interpolation factor
% output:
    % RR: [nz,Nx,ny]
function RR = DAT2RR_GPU(DAT, intpDk)
[nk,Nx,ny] = size(DAT); % Dat: spectrum data, [nk, Nx, ny]
nz = round(nk/2);
RR=single(zeros(nz,Nx,ny)); % initialize RR
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
k = linspace(1-intpDk/2, 1+intpDk/2, nk);
lam = single(gpuArray(1./fliplr(k)));
X=fliplr(lam);
XQ=fliplr(linspace(min(lam),max(lam),length(lam)));
DatGPU=gpuArray(DAT);
DC=flip(single(mean(DatGPU(:,:),2)));
DatGPU=flip(DatGPU,1);
for iChk=1:nyChk
    V = DatGPU(:,:,(iChk-1)*nyPchk+1:iChk*nyPchk) - DC; % substract the reference signal, Subtract mean
    %% lamda to k space
    data_k = flip(interp1(X, V, XQ, 'linear'),1);
    clear V
    %% ifft
    RRy(:,:,:)= ifft((data_k(:,:,:)),[],1);
    clear data_k
    RR(:,:,(iChk-1)*nyPchk+1:iChk*nyPchk)=gather(RRy(1:nz,:,:));
end

