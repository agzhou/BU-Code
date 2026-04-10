% demoRun
%
%  A simple script that can be used to demonstrate the fastSave tools.  It performs the following tasks:
%   1. Checks that the fastsave .mex file appropriate for the OS is on the path
%   2. Checks for optimal formatting parameters for the drives and verifies that they are mounted
%   3. Erases previously collected data on all of the drives.
%   4. Trim NVME drives
%   5. Run VSX
%   6. Run validation script

% Copyright (C) 2001-2025, Verasonics, Inc.
% All worldwide rights and remedies under all intellectual
% property laws and industrial property laws are reserved.

SetUpL22_14vXFlashAngles_fastSave
filename = 'SetUpL22_14vXFlashAngles_fastSave';

%SetUpVermon1024_8MHz_FlashAngles_Light_superframe
%filename = 'MatFiles/Vermon1024_8MHz-FlashAngles-Light_superframe.mat';


%% Check drives and wipe them for best performance
if isunix()
    fsDriveNames = {'/media/verasonics/WD1', '/media/verasonics/WD2', '/media/verasonics/WD3', '/media/verasonics/WD4'};
elseif ispc()
    fsDriveNames = {'E:\', 'F:\', 'G:\','H:\'};
end

mexfileName = ['fastsave.' mexext];
if ~exist(mexfileName, 'file')
    error('%s could not be found on the path.  Please verify the file location is on the path.', mexfileName)
end
% Test drive formatting
disp('[FastSave] Checking formatting type')
for driveI = 1:length(fsDriveNames)
    testDrive(fsDriveNames{driveI});
end

%  3. Delete previiously stored data on destination drive.
answer = questdlg('Delete all RF data from storage devices?', 'Question', 'Yes', 'No', 'No');
if strcmp(answer, 'Yes')
    disp('[FastSave] Clearing disk space of fastwrite drives')

    for driveI = 1:length(fsDriveNames)
        delete([fsDriveNames{driveI} '/*.rf'])
    end
end


% 4. Trim NVMe Drives
disp('[FastSave] TRIM fastwrite drives - (requires password for sudo permission)')
for driveI = 1:length(fsDriveNames)
    fprintf('TRIMing %s\n',fsDriveNames{driveI});
if isunix()
    trimString = sprintf("sudo fstrim %s -v", fsDriveNames{driveI});
elseif ispc()
    trimString = sprintf("defrag %s /L ", fsDriveNames{driveI});
end
    system(trimString);
end

%% Pop up a resource manager for either linux or windows
if isunix()
    system('export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH; gnome-terminal --geometry=139x23--26+87 -- bash -c "watch df; exec bash"');
elseif ispc()
    system('taskmgr &');
end

% Run script
VSX
% Run validation script
validateCollectedData

%% Test the destination path to make sure that it is mounted and in the best format
function testDrive(drivePath)
[driveLetter, ~, ~] = fileparts(drivePath);
driveLetter = driveLetter(1);
if isunix()
    [~, result] = system(['lsblk -f |grep ' driveLetter]);
    if isempty(result)
        error(['Drive ' driveLetter ' is not mounted'])
    elseif ~contains(result,'ext4')
        error(['Drive ' driveLetter ' is not formatted as ext4.  Disk write speed will be slower than expected'])
    end
elseif ispc()

    [~, result] = system('powershell -Command "Get-CimInstance -ClassName Win32_LogicalDisk | select-object DeviceID, FileSystem"');
    driveTable = splitlines(result);
    isMounted = any(contains(driveTable, driveLetter));
    if (isMounted)
        isNTFS = contains(driveTable(contains(driveTable, driveLetter)),'NTFS');
        if (~isNTFS)
            error('Drive %s is not mounted as NTFS', driveLetter)
        end
    else
        error('Drive %s is not mounted', driveLetter)
    end
end
end
