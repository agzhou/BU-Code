
% Adapting Brad Rauscher's code
% Input: one frame of pupil data (grayscale image)

function ROIs = defineROIs(frame)

    ff = figure; % frame figure
    imagesc(frame); % colormap gray
    title('Draw a ROI around the eye');
    eye_ROI = drawassisted('Color', 'r'); % Draw the eye ROI
    pause;
    
    ROIs.eye_mask = createMask(eye_ROI); % Create a mask from the eye ROI definition
    
    % Draw a point in the pupil to get the intensity there
    clf; imagesc(frame); % Reset the figure
    title('Draw point in the pupil');
    pupilpos = drawpoint('Color', 'r'); % Pick a point in the pupil and store its location
    pause;
    % point.pupil.pos = pupilpos.Position;
    pupil_intensity = frame(round(pupilpos.Position(2)), round(pupilpos.Position(1))); % Store the intensity at the picked pupil point
    
    clf; imagesc(frame);
    title('Draw point in the iris');
    irispos = drawpoint('Color', 'r'); % Pick a point in the iris and store its location
    pause;
    % point.iris.pos = irispos.Position;
    ival = imageData(round(irispos.Position(2)),round(irispos.Position(1)));
    
    ROIs.pupil_thrVal = (pval+ival)/2;
    
    clf; imshow(imageData);
    title('Draw line across eye (horizontal)');
    eyelength = drawline('Color','r');
    pause;
    ROIs.eye_length = ((eyelength.Position(1,1)-eyelength.Position(2,1))^2+(eyelength.Position(1,2)-eyelength.Position(2,2))^2)^0.5;
    
    clf; imshow(imageData);
    title('Select ROI around long whiskers');
    r = drawrectangle('Color','r');
    roi1 = r.Position;
    pause;
    
    clf; imshow(imageData);
    title('Select ROI around whisker pad');
    r2 = drawrectangle('Color','r');
    roi2 = r2.Position;
    pause;
    
    ROIs.whisker_rois = [roi1; roi2];
    
    close(bFig);

end