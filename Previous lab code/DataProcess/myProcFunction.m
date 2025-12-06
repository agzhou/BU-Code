function myProcFunction(RData)
persistent myHandle
% If ‘myPlotChnl’ exists, read it for the channel to plot.
if evalin('base','exist(''myPlotChnl'',''var'')')
    channel = evalin('base','myPlotChnl');
else
    channel = 32; % Channel no. to plot
end
% Create the figure if it doesn’t exist.
if isempty(myHandle)||~ishandle(myHandle)
    figure;
    myHandle = axes('XLim',[0,1500],'YLim',[-16384 16384], ...
        'NextPlot','replacechildren');
end
% Plot the RF data.
plot(myHandle,RData(:,channel));
drawnow
return