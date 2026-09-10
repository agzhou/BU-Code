%
% File name: processFIOresults.m
%            Import data from text file script for importing data from the
%            fio test script generated file and plotting in Matlab.
%
% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage NXT Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%
% Copyright © 2013-2025 Verasonics, Inc.


% Copyright (C) 2001-2025, Verasonics, Inc.
% All worldwide rights and remedies under all intellectual
% property laws and industrial property laws are reserved.

close all
%filename = "/home/verasonics/git/bryan/vantage-sw-2/Example_Scripts/Community_Portal/FastSave/LINUX_RAID0_4TBSN850x_2TB_nvme_write_benchmark.csv"
%driveString = "RAID WD BLACK SN850"

%filename = "/home/verasonics/git/bryan/vantage-sw-2/Example_Scripts/Community_Portal/FastSave/LINUX_1TBWDSN8100_2GB_nvme_write_benchmark.csv"
%driveString = "WD BLACK SN8100"

filename = "./LINUX_WD850_XFS_nvme_write_benchmark.csv"
driveString = "WD BLACK SN850 - XFS"

osString = "Linux"

FileSizeString = "2GB"

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 3, "Encoding", "UTF-8");

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["time_seconds", "dataWritten_GB", "writespeed_MBps"];
opts.VariableTypes = ["double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
fioresults = readtable(filename, opts);


%% Clear temporary variables
clear opts

plot(fioresults.dataWritten_GB,fioresults.writespeed_MBps)
titleString = sprintf("os: %s, drive: %s, fileSize: %s", osString, driveString, FileSizeString)
title(titleString)
xlabel('Data Written (GB)')
ylabel('Write Speed (MB/s)')

filename = sprintf('FIO_%s_%s_%s.png',osString,driveString,FileSizeString);
print('-dpng',filename)
