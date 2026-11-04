%% PROPOSED replacement for the "Wavelet transform --> scalogram" section of LFO3D.m (NOT applied to the original)
% Needs cwtPrepare.m and cwtBatchedCore.m on the path.
% Inputs expected from the rest of LFO3D.m: PDIallSF, vs, num_voxels, sfRate, nsf

fb = cwtfilterbank('SignalLength', nsf, 'SamplingFrequency', sfRate, 'Wavelet', 'amor', 'Boundary', 'reflection');
%   Optional, big speedup + memory saving if you only care about the LFO band:
%   ..., 'FrequencyLimits', [fLow fHigh])

cwtOutput = "abs";     % "complex" (phase kept) | "abs" | "abs2"; abs is ~2x faster and half the memory
useGPU    = gpuDeviceCount > 0;
chunkSize = 5000;      % voxels per batch; lower this if the GPU/RAM runs out

PDIstackedT = single(reshape(PDIallSF{3}, [num_voxels, nsf])).';   % [time x voxels], time down the columns
P = cwtPrepare(fb, 'single');                                      % filters pulled out once (self-checks against cwt)
faxis = P.faxis; num_freqs = numel(faxis);
if useGPU, P.psiT = gpuArray(P.psiT); end

voxelsToRun = 1:num_voxels;           % e.g. replace with find(brainMask(:)) to skip voxels outside the brain
if cwtOutput == "complex", outCls = complex(zeros(1, 'single')); else, outCls = zeros(1, 'single'); end
PDI_scalogram_T = zeros(nsf, numel(voxelsToRun), num_freqs, 'like', outCls);   % [time x voxel x freq]

for s = 1:chunkSize:numel(voxelsToRun)
    ii = s:min(s+chunkSize-1, numel(voxelsToRun));
    Xc = PDIstackedT(:, voxelsToRun(ii));
    if useGPU, Xc = gpuArray(Xc); end
    PDI_scalogram_T(:, ii, :) = gather(cwtBatchedCore(Xc, P, cwtOutput));
end
% If downstream code needs the old [voxel x freq x time] layout:
%   PDI_scalogram = permute(PDI_scalogram_T, [2 3 1]);   % makes a full copy -> mind memory
