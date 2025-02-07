% Try the simpletracker.m by Jean-Yves Tinevez

%%
addpath('C:\Users\BOAS-US\Documents\Allen\GitHub\BU-Code\Allen code\Processing\simpletracker\src')
addpath('C:\Users\BOAS-US\Documents\Allen\GitHub\BU-Code\Allen code\Processing\munkres')

load('G:\Allen\Data\01-29-2025 AZ001 ULM\RC15gV\run 1 left eye\Processed Data\allCentroids.mat')

max_linking_distance = Inf; % maximum spatial distance that we track a particle for between frames
max_gap_closing = 3;    % maximum # of frames to check for gap filling
debug = true;

points = allCentroids;

[ tracks, adjacency_tracks ] = simpletracker(points,...
    'MaxLinkingDistance', max_linking_distance, ...
    'MaxGapClosing', max_gap_closing, ...
    'Debug', debug);