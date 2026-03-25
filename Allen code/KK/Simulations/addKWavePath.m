%% Add the k-Wave path automatically
% (using my GitHub structure: clone the k-wave repo from https://github.com/ucl-bug/k-wave)
% Note: to use the C++ code, download from
% http://www.k-wave.org/download.php and put the executables into the
% k-Wave 'binaries' folder

function addKWavePath()
    tmp = matlab.desktop.editor.getActive;
    fp = fileparts(tmp.Filename);
    fp_cell = regexp(fp, filesep, 'split');

    addpath(genpath(fullfile(fp_cell{1:find(contains(fp_cell, "GitHub"), 1)}) + "\k-wave"));    
    % addpath(genpath(fullfile(fp_cell{1:find(contains(fp_cell, "GitHub"), 1)}) + "\k-Wave"));
    % addpath(genpath(fullfile(fp_cell{1:find(contains(fp_cell, "BU-Code"), 1)}) + "\k-Wave"));

end