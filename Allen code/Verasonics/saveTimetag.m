% Note: savefast only works well with numerical arrays. see https://www.mathworks.com/matlabcentral/fileexchange/39721-save-mat-files-more-quickly


function saveTimetag()
%     tic
    savepath = evalin('base', 'savepath');
    
    filename = strcat('startTimeTag', '.mat');

    timetag = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS');
    savefast([savepath, filename], 'timetag')

end