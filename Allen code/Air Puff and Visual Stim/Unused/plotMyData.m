% From Matlab documentation
% (https://www.mathworks.com/help/daq/acquire-data-in-the-background.html)

function plotMyData(obj,evt)
% obj is the DataAcquisition object passed in. evt is not used.
     data = read(obj,obj.ScansAvailableFcnCount,"OutputFormat","Matrix");
     plot(data)
 end