%% Load CR and RC IQ volumes

%% SVD filter CR and RC IQ volume stacks separately
[xpo, ypo, zpo, nfo] = size(IQ_CR); % Original sizes
% figure; imagesc(squeeze(max(abs(IQ_CR(:, :, :, 1)), [], 1)).')

xrange = 1:xpo;
yrange = 1:ypo;
zrange = 30:130;
% zrange = 1:zpo;
% zrange = zstart:zend;
voxelRange = {xrange, yrange, zrange};
xp = length(xrange); yp = length(yrange); zp = length(zrange); nf = size(IQ_CR, 4);
    
% Mask the IQ to some region (according to the function input)
xrange = voxelRange{1};
yrange = voxelRange{2};
zrange = voxelRange{3};
IQ_CR = IQ_CR(xrange, yrange, zrange, :);
IQ_RC = IQ_RC(xrange, yrange, zrange, :);

% [brainMask, maskInfo] = createBrainVoxelMask(sum(abs(IQ), 4), 'MarginVoxels', 0);
% IQ = IQ .* brainMask;   % zero out non-brain voxels before the SVD sees them

% SVD decluttering
% [xp, yp, zp, nf] = size(IQ);
svt_low = 15; % Lower singular value threshold
IQf_RC = SVD_CF(IQ_RC, svt_low, nf);
IQf_CR = SVD_CF(IQ_CR, svt_low, nf);

IQf_OPW = IQf_RC + IQf_CR;
PDI_OPW = mean(abs(IQf_OPW).^2, 4);

XD = IQf_RC .* conj(IQf_CR);
% PDI_XD = mean(abs(XD).^2, 4);
PDI_XD = 2.*real(mean(XD, 4));
PDI_XD = abs(mean(IQf_RC .* conj(IQf_CR) + IQf_CR .* conj(IQf_RC), 4));
%% Plot
PDI_OPW_xz_MIP = squeeze(max(PDI_OPW, [], 1)).';
PDI_XD_xz_MIP = squeeze(max(PDI_XD, [], 1)).';
figure; imagesc(PDI_OPW_xz_MIP .^ 0.5); colormap hot
figure; imagesc(PDI_XD_xz_MIP .^ 0.5); colormap hot

figure;
yyaxis left
plot(squeeze(PDI_OPW_xz_MIP(30, :)))
yyaxis right
plot(squeeze(PDI_XD_xz_MIP(30, :)))
% hold off