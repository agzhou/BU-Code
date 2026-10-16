% SVD clutter filter
function [IQf] = SVD_CF(IQ, sv_threshold_lower, sv_threshold_upper)
    d = length(size(IQ)); % Dimensionality of the IQ data (3 --> 2D + time, 4 --> 3D + time)
    switch d
        case 3
            [CM, EVs, V] = getSVs1D(IQ);
            % disp('SVs decomposed')
            % [IQf, noise] = applySVs1D_accel(IQ, CM, EVs, V, sv_threshold_lower, sv_threshold_upper);
            [IQf] = applySVs1D_accel(IQ, CM, EVs, V, sv_threshold_lower, sv_threshold_upper);

        case 4
            [CM, EVs, V] = getSVs2D(IQ);
            % disp('SVs decomposed')
            % [IQf, noise] = applySVs2D_accel(IQ, CM, EVs, V, sv_threshold_lower, sv_threshold_upper);
            [IQf] = applySVs2D_accel(IQ, CM, EVs, V, sv_threshold_lower, sv_threshold_upper);

    end
end