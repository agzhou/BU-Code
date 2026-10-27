function [frameRate] = checkFrameRate(P)
    %% Check the frame rate - hardware could lag behind the nominal one
    % Load RF time tags: RFtimeTags.mat
    % if ~exist('RFtimeTags', 'var')
        % Choose and load the RFTimeTagsPerChannel.mat file (from the acquisition and running readRFTimeTags.m)
        [RFTT_filename, RFTT_pathname, ~] = uigetfile('*.mat', 'Select the RF timetags file', ['RFTimeTagsPerChannel.mat']);
        load([RFTT_pathname, RFTT_filename], 'RFTimeTags_allSFStacked')
    % end
    
    dt_all = diff(squeeze(mean(RFTimeTags_allSFStacked, 1)), 1); % Difference in frame timetags
    dt_subframe = dt_all(dt_all < 0.8*(1/P.sfRate - P.numFramesPerBuffer/P.frameRate)); % Remove the large inter-superframe jumps
    dt = mean(dt_subframe);
    frameRate_actual = 1/dt; % Actual frame rate
    
    % Set the downstream-used frameRate to the actual frame rate
    if abs(P.frameRate - frameRate_actual) > 25
        warning("Nominal frame rate (" + num2str(P.frameRate) + " Hz) was not supported by the hardware during data acquisition. Using the actual frame rate of " + num2str(frameRate_actual) + " Hz instead.")
        frameRate = frameRate_actual;
    else
        frameRate = P.frameRate;
    end
end