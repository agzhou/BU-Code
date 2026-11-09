% No frequency separation

function [PDI, g1, noise] = IQ2g1T_3D_nofreqsep(IQ, voxelRange, sv_threshold_lower, sv_threshold_upper, HPF, nTau)

    % Mask the IQ to some region (according to the function input)
    xrange = voxelRange{1};
    yrange = voxelRange{2};
    zrange = voxelRange{3};
    IQ = IQ(xrange, yrange, zrange, :);

    % [brainMask, maskInfo] = createBrainVoxelMask(sum(abs(IQ), 4), 'MarginVoxels', 0);
    % IQ = IQ .* brainMask;   % zero out non-brain voxels before the SVD sees them
    
    % SVD decluttering
    % [xp, yp, zp, nf] = size(IQ);
    
    [CM, EVs, V] = getSVs2D(IQ);
    % disp('SVs decomposed')
    % [IQf, noise] = applySVs2D_accel(IQ, CM, EVs, V, sv_threshold_lower, sv_threshold_upper);
    [IQf] = applySVs2D_accel(IQ, CM, EVs, V, sv_threshold_lower, sv_threshold_upper);
    % disp('SVD filtered images put together')
    clearvars CM EVs V IQ

%     figure; imagesc(squeeze(max(abs(IQ(:, :, :, 1)), [], 1))' .^ 0.5)

    % High pass filter (apply to the post-SVD clutter filtered data)
    HPF.dim = length(size(IQf)); % Operate on the time dimension
    IQf_HPF = filter(HPF.b, HPF.a, IQf, [], HPF.dim);
    clearvars IQf
    [PDI] = sum(abs(IQf_HPF).^2, 4)./size(IQf_HPF, 4);
    g1 = g1T_fft(IQf_HPF, nTau); % Use the base filtered IQ


    % [PDI] = sum(abs(IQf).^2, 4)./size(IQf, 4);
    % g1 = g1T_fft(IQf, nTau); % Use the base filtered IQ
    
end