% Description: load binary .rf files, e.g., those saved by Verasonics'
% FastSave.

% Inputs:
%   fileName: Complete path of a .rf file, including the .rf extension
%   P: parameter struct, which only needs to include the Resource.RcvBuffer info
function [RFData] = loadBinaryRF(fileName, P)
    fid = fopen(fileName, 'rb'); 
    RFData = fread(fid, [P.Resource.RcvBuffer(1).rowsPerFrame, P.Resource.RcvBuffer(1).colsPerFrame], '*int16');
    fclose(fid);
end