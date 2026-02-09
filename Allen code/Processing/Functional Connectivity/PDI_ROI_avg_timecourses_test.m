%% Plot each ROI's PDI timecourse
% pmax = max(PDI_ROI_GMS_timecourses_mat(:, 1));
% pmin = min(PDI_ROI_GMS_timecourses_mat(:, 1));

% pmax = prctile(PDI_ROI_GMS_timecourses_mat, 80, 'all');
% pmin = prctile(PDI_ROI_GMS_timecourses_mat, 40, 'all');

% pmax = mean(abs(prctile(PDI_ROI_GMS_timecourses_mat, 90, 1)), 'omitnan');
% pmin = mean(prctile(PDI_ROI_GMS_timecourses_mat, 10, 1), 'omitnan');
mr = mean(abs(prctile(PDI_ROI_GMS_timecourses_mat, 95, 1) - prctile(PDI_ROI_GMS_timecourses_mat, 5, 1)), 'omitnan'); % mean range..

% Plot multiple stackedplots to visualize ROI PDI timecourses
num_cols_per_sp = 16;
num_sps = ceil(roi.num_regions/num_cols_per_sp); % # of stackedplot to use since they only allow 25 columns max

for spi = 1:num_sps
    % - NOTE: stackedplot only allows for 25 columns max - % 
    % Normal
    figure
    temp_ind_spi = (spi - 1)*num_cols_per_sp + 1:spi*num_cols_per_sp;
    % % ROI_PDI_timecourse_sp = stackedplot(t, PDI_ROI_timecourses_mat, 'DisplayLabels', region_acronyms);
    % ROI_PDI_timecourse_sp = stackedplot(t, [PDI_ROI_timecourses_mat(:, temp_ind_spi), GVTD], 'DisplayLabels', [roi.acronyms(temp_ind_spi); {'GVTD'}]);
    ROI_PDI_timecourse_sp = stackedplot(t, [PDI_ROI_GMS_timecourses_mat(:, temp_ind_spi) - mean(PDI_ROI_GMS_timecourses_mat(:, temp_ind_spi), 1, 'omitnan'), GVTD], 'DisplayLabels', [roi.acronyms(temp_ind_spi); {'GVTD'}]);
    title("ROI average PDI timecourses")
    xlabel("Time [s]")
    % fontsize(14, 'points')
    % ylim([pmin, pmax])
    for ai = 1:size(ROI_PDI_timecourse_sp.AxesProperties, 1) - 1 % Go through the index for each axes object and equalize scaling (except GVTD)
        ROI_PDI_timecourse_sp.AxesProperties(ai).YLimits = [-mr, mr];
    end
end

for spi = 1:num_sps*2
    % - NOTE: stackedplot only allows for 25 columns max - % 
    
    % Hemisphere-separated
    figure
    temp_ind_spi = (spi - 1)*num_cols_per_sp + 1:spi*num_cols_per_sp;

    % ROI_PDI_timecourse_sp = stackedplot(t, [PDI_ROI_hemis_timecourses_mat(:, temp_ind_spi), GVTD], 'DisplayLabels', [roi.acronyms_hemis_interleaved(temp_ind_spi); {'GVTD'}]);
    ROI_PDI_timecourse_sp = stackedplot(t, [PDI_ROI_hemis_GMS_timecourses_mat(:, temp_ind_spi) - mean(PDI_ROI_hemis_GMS_timecourses_mat(:, temp_ind_spi), 1, 'omitnan'), GVTD], 'DisplayLabels', [roi.acronyms_hemis_interleaved(temp_ind_spi); {'GVTD'}]);
    title("ROI average (hemisphere-separated) PDI timecourses")
    xlabel("Time [s]")
    % fontsize(14, 'points')
    for ai = 1:size(ROI_PDI_timecourse_sp.AxesProperties, 1) - 1 % Go through the index for each axes object and equalize scaling (except GVTD)
        ROI_PDI_timecourse_sp.AxesProperties(ai).YLimits = [-mr, mr];
    end
end
