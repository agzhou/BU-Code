% Input: IQf_separated (3x1 cell array with the negative, positive, and all
% frequency components of the filtered IQ (IQf)

% Output: CDI (same except Color Doppler index)
function [CDI] = calcColorDoppler(IQf_FT_separated, P)

    % Determine thresholds on the Fourier spectrum to get a more accurate
    % Doppler frequency estimation
    max_vz_expected = 0.05; % [m/s]
    abs_fD_max_expected = 2 * (P.Trans.frequency * 1e6) * max_vz_expected / P.Resource.Parameters.speedOfSound;

    CDI = cell(size(IQf_FT_separated));

    freqDim = length(size(IQf_FT_separated{1})); % Usually the frequency dimension is the last dimension. 3 for 2D data and 4 for 3D data.
    for i = 1:length(CDI)
%     for i = 1
        switch i
            case 1 % Negative frequencies
                Flim = [-P.frameRate/2, 0];
            case 2 % Positive
                Flim = [0, P.frameRate/2];
            case 3 % All
                Flim = [-P.frameRate/2, P.frameRate/2];
        end
        fi = linspace(Flim(1), Flim(2), size(IQf_FT_separated{i}, freqDim)); % frequencies
        
        % Note: this line only works for 3D right now
        fi = repmat(permute(fi, [1, 4, 3, 2]), size(IQf_FT_separated{i}, 1), size(IQf_FT_separated{i}, 2), size(IQf_FT_separated{i}, 3), 1);
%         fDi = sum( abs(IQf_FT_separated{i} .^ 2) .* fi, freqDim) ./ sum( abs(IQf_FT_separated{i} .^ 2), freqDim); % Estimate the Doppler frequency
%         fDi = sum( abs(IQf_FT_separated{i}) .^ 2 .* fi, freqDim) ./ sum( abs(IQf_FT_separated{i}) .^ 2, freqDim); % Estimate the Doppler frequency
        
        % Apply the threshold
%         switch i
%             case 1
%                 FTsegment = abs(IQf_FT_separated{i});
%                 FTsegment = FTsegment(fi < -abs_fD_max_expected);
%                 
%             case 2
%                 Flim = [0, P.frameRate/2];
%             case 3
%                 Flim = [-P.frameRate/2, P.frameRate/2];
%         end
        FTsegment = abs(IQf_FT_separated{i});
        mask = abs(fi) > abs_fD_max_expected;
        FTsegment(~mask) = 0;

        threshold = 0.5 * max(FTsegment, [], freqDim);
        
        IQf_FT_thresholded = abs(IQf_FT_separated{i});

        IQf_FT_thresholded(mask) = 0; % Threshold method 1 (throw out high |freq| noise)
%         IQf_FT_thresholded(IQf_FT_thresholded < threshold) = 0; % The
%         magnitude based method in the Improved Color Doppler paper

%         %%%%%%%%%%%%%%%%% currently not thresholding %%%%%%%%%%%%%%%%%%%%%
        fDi = sum( abs(IQf_FT_thresholded) .^ 2 .* fi, freqDim) ./ sum( abs(IQf_FT_thresholded) .^ 2, freqDim); % Estimate the Doppler frequency

        CDI{i} = -fDi * P.Resource.Parameters.speedOfSound / 2 / (P.Trans.frequency * 1e6);
    end

end