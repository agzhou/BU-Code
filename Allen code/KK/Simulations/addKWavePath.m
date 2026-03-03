%% Add the k-Wave path automatically
% (using my GitHub structure: clone the k-wave repo)
function addKWavePath()
    tmp = matlab.desktop.editor.getActive;
    fp = fileparts(tmp.Filename);
    fp_cell = regexp(fp, filesep, 'split');
    
    addpath(fullfile(fp_cell{1:find(contains(fp_cell, "GitHub"), 1)}) + "\k-wave");
end