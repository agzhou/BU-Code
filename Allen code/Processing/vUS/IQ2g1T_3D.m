function [PDI, CDI, g1] = IQ2g1T_3D(IQ, P, voxelRange, sv_threshold_lower, sv_threshold_upper, HPF, nTau)

    % Mask the IQ to some region (according to the function input)
    xrange = voxelRange{1};
    yrange = voxelRange{2};
    zrange = voxelRange{3};
    IQ = IQ(xrange, yrange, zrange, :);
    
    % SVD decluttering
    % [xp, yp, zp, nf] = size(IQ);
    
    [CM, EVs, V] = getSVs2D(IQ);
    % disp('SVs decomposed')
    [IQf, noise] = applySVs2D(IQ, CM, EVs, V, sv_threshold_lower, sv_threshold_upper);
    % disp('SVD filtered images put together')
    clearvars CM EVs V IQ

%     figure; imagesc(squeeze(abs(IQf(:, :, 1))) .^ 0.5)

    % High pass filter (apply to the post-SVD clutter filtered data)
    HPF.dim = length(size(IQf)); % Operate on the time dimension
    IQf_HPF = filter(HPF.b, HPF.a, IQf, [], HPF.dim);
    clearvars IQf

    % Use the IQf with separated negative and positive frequency components
    [IQf_separated, IQf_FT_separated, nFTpts] = separatePosNegFreqs(IQf_HPF); % Outputs are cell arrays in the order of: negative, positive, all frequencies
    [PDI] = calcPowerDoppler(IQf_separated, noise);
    [CDI] = calcColorDoppler(IQf_FT_separated, P);

    % PDI = sum(abs(IQf) .^ 2, 3) ./ size(IQf, 3);
    % PDI = sum(abs(IQf) .^ 2, 3) ./ size(IQf, 3) ./ noise;
    % figure; imagesc(x_mm, z_mm, squeeze(PDI .^ 0.5)); colormap hot; colorbar; title('Power Doppler'); xlabel('x [mm]'); ylabel('z [mm]')
    % figure; imagesc(x_mm, z_mm, squeeze(abs(IQ(:, :, 1)))); colorbar; title('IQ'); xlabel('x [mm]'); ylabel('z [mm]')

    % Calculate g1
    % Store g1 for each frequency component in a cell array
    g1 = cell(size(IQf_separated));
    
    for j = ctp
        g1{j} = g1T_fft(IQf_separated{j}, nTau); % Use the base filtered IQ
        % g1{j} = g1T_fft(IQf_separated_masked{j}, nTau); % Use the filtered IQ with system noise removed
    end
end