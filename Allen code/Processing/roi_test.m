%% Quick test on 12-17-2025 to look at 12-12-2025 RCA fUS data

% ROI where we expect the activation
roi_x = 20:40;
roi_y = 60:80;
roi_z = 1:20;

[x, y, z] = meshgrid(roi_x, roi_y, roi_z);

roi_indices = sub2ind(size(PDI), x(:), y(:), z(:));
roi_PDI = calc_ROI_avg(PDIallSF, roi_indices);
roi_CBVi = calc_ROI_avg(CBViallSF, roi_indices);
roi_CBFsi = calc_ROI_avg(CBFsiallSF, roi_indices);

% for ti = 1:P.numTrials
%     roi_PDI_usi{ti} = calc_ROI_avg(trial_PDI_usi{ti}, roi_indices);
%     roi_CBVi_usi{ti} = calc_ROI_avg(CBViallSF{ti}, roi_indices);
%     roi_CBFsi_usi{ti} = calc_ROI_avg(CBFsiallSF{ti}, roi_indices);
% end
%%
% figure; plot(roi_PDI); title("ROI PDI"); xlabel("Superframe index")
% figure; plot(roi_CBVi); title("ROI CBVi"); xlabel("Superframe index")
% figure; plot(roi_CBFsi); title("ROI CBFsi"); xlabel("Superframe index")

%% Plot ROI data with actual time
figure; plot(t, roi_PDI, 'LineWidth', 2); title("ROI PDI"); xlabel("Time [s]")
hold on
xline(stimOnsets, 'r', 'LineWidth', 2)
hold off
