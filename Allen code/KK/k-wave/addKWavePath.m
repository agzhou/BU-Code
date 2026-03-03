%% Add the k-Wave path
% tmp = matlab.desktop.editor.getActive;
fp = fileparts(tmp.Filename)
% cd();
test = regexp(fp, filesep, 'split');

% filename = mfilename;
fullfile(test{1:find(contains(test, "k-wave"), 1)});