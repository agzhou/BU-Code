% Create a large random data array (e.g., 1 GB of double-precision numbers)
dataSizeMB = 1000; 
numElements = (dataSizeMB * 1024 * 1024) / 8;
largeData = rand(1, numElements);

filename = 'F:\test_speed.mat';

% Time the write operation with compression disabled
tic;
% save(filename, 'largeData', '-v7.3', '-nocompression');
savefast(filename, 'largeData');
elapsedTime = toc;

% Calculate speed in MB/s
speedMBps = dataSizeMB / elapsedTime;
fprintf('Write Speed: %.2f MB/s\n', speedMBps);

% Clean up the file
delete(filename);