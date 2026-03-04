%% Add the MUST path automatically
% (using my GitHub structure: I put the MUST toolbox in the GitHub/BU-Code/ directoryrepo)
function addMUSTPath()
    tmp = matlab.desktop.editor.getActive;
    fp = fileparts(tmp.Filename);
    fp_cell = regexp(fp, filesep, 'split');
    
    addpath(genpath(fullfile(fp_cell{1:find(contains(fp_cell, "BU-Code"), 1)}) + "\MUST"));
end