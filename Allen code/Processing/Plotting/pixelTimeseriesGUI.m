function pixelTimeseriesGUI(Data, dispImg, varargin)
% PIXELTIMESERIESGUI  Click a pixel on a 2D image to plot its timeseries.
%
%   pixelTimeseriesGUI(Data) opens a figure with an image on the left,
%   the pixel's timeseries at top right, and (for complex Data) a
%   Re(signal) vs Im(signal) trajectory plot at bottom right. Data is a
%   [nz x nx x nt] matrix (rows x cols x time), i.e. Data(row,col,:) is
%   the timeseries at pixel (row,col). Click (and drag) on the image to
%   inspect pixels. The displayed image defaults to the std-over-time
%   projection of Data.
%
%   pixelTimeseriesGUI(Data, dispImg) displays dispImg, a separate
%   [nz x nx] image, instead of the default projection. dispImg may
%   have a different pixel grid size than Data (e.g. a higher-resolution
%   B-mode background behind coarser flow/Doppler data) -- clicks on
%   dispImg are mapped proportionally onto Data's grid. Pass [] to use
%   the default projection while still supplying other options below.
%
%   pixelTimeseriesGUI(...,'Time',t) uses t as the x-axis for the
%   timeseries plot (default 1:nt).
%
%   pixelTimeseriesGUI(...,'ComplexMode',mode) controls how complex
%   Data is plotted: 'realimag' (default, plots Re and Im as separate
%   lines), 'abs', 'angle', 'real', or 'imag'. Ignored for real Data.
%
%   Additional Name/Value pairs: 'CLim', 'Colormap'.

if nargin<2
    dispImg = [];
end

p = inputParser;
addParameter(p,'Time',[]);
addParameter(p,'CLim',[]);
addParameter(p,'Colormap','gray');
addParameter(p,'ComplexMode','realimag');
parse(p,varargin{:});

[nz,nx,nt] = size(Data);

t = p.Results.Time;
if isempty(t)
    t = 1:nt;
end

if isempty(dispImg)
    dispImg = std(double(Data),0,3);
end
[nzDisp,nxDisp] = size(dispImg);

fig = figure('Name','Pixel Timeseries Viewer','NumberTitle','off');

imgAx = subplot(2,2,[1 3],'Parent',fig);
imH = imagesc(imgAx, dispImg);
axis(imgAx,'image');
colormap(imgAx, p.Results.Colormap);
colorbar(imgAx);
if ~isempty(p.Results.CLim)
    set(imgAx,'CLim',p.Results.CLim);
end
hold(imgAx,'on');
markerH = plot(imgAx, NaN, NaN, 'r+', 'MarkerSize',12,'LineWidth',2);
title(imgAx,'Click a pixel to view its timeseries');
xlabel(imgAx,'X (col)'); ylabel(imgAx,'Z (row)');

tsAx = subplot(2,2,2,'Parent',fig);
hold(tsAx,'on');
lineH1 = plot(tsAx, t, NaN(size(t)),'DisplayName','Signal');
lineH2 = plot(tsAx, t, NaN(size(t)),'DisplayName','Im','Visible','off');
xlabel(tsAx,'Time');
ylabel(tsAx,'Signal');
title(tsAx,'Click a pixel on the image');
grid(tsAx,'on');

riAx = subplot(2,2,4,'Parent',fig);
riLineH = plot(riAx, NaN, NaN, '.-');
axis(riAx,'equal');
xlim(riAx,[-1 1]);
ylim(riAx,[-1 1]);
xlabel(riAx,'Re(Signal)');
ylabel(riAx,'Im(Signal)');
title(riAx,'Re vs Im (complex data only)');
grid(riAx,'on');

vars = struct('Data',Data,'Time',t,'ImgAxes',imgAx,'TSAxes',tsAx,'RIAxes',riAx,...
    'Marker',markerH,'Line1',lineH1,'Line2',lineH2,'RILine',riLineH,'Fig',fig,...
    'DispSize',[nzDisp,nxDisp],'DataSize',[nz,nx],...
    'ComplexMode',p.Results.ComplexMode);

set(imH,'ButtonDownFcn',{@imageClick,vars});
set(fig,'WindowButtonUpFcn',@(~,~) set(fig,'WindowButtonMotionFcn',''));

end

% Callback subfunctions to support UI actions
function imageClick(~,~,vars)
    updatePlot(vars);
    set(vars.Fig,'WindowButtonMotionFcn',{@dragMotion,vars});
end

function dragMotion(~,~,vars)
    updatePlot(vars);
end

function updatePlot(vars)
    pt = get(vars.ImgAxes,'CurrentPoint');
    colDisp = round(pt(1,1));
    rowDisp = round(pt(1,2));
    nzDisp = vars.DispSize(1); nxDisp = vars.DispSize(2);
    if colDisp<1 || colDisp>nxDisp || rowDisp<1 || rowDisp>nzDisp
        return
    end

    % map click location on the (possibly higher/lower-res) display
    % image onto Data's own pixel grid
    nz = vars.DataSize(1); nx = vars.DataSize(2);
    col = min(nx, max(1, round(colDisp/nxDisp*nx)));
    row = min(nz, max(1, round(rowDisp/nzDisp*nz)));

    ts = squeeze(vars.Data(row,col,:));

    if isreal(ts)
        set(vars.Line1,'YData',ts,'DisplayName','Signal','Visible','on');
        set(vars.Line2,'Visible','off');
        ylabel(vars.TSAxes,'Signal');
        legend(vars.TSAxes,'off');
        set(vars.RILine,'XData',NaN,'YData',NaN);
    else
        set(vars.RILine,'XData',real(ts),'YData',imag(ts));
        switch vars.ComplexMode
            case 'abs'
                set(vars.Line1,'YData',abs(ts),'DisplayName','|Signal|','Visible','on');
                set(vars.Line2,'Visible','off');
                ylabel(vars.TSAxes,'|Signal|');
                legend(vars.TSAxes,'off');
            case 'angle'
                set(vars.Line1,'YData',angle(ts),'DisplayName','angle(Signal)','Visible','on');
                set(vars.Line2,'Visible','off');
                ylabel(vars.TSAxes,'Phase (rad)');
                legend(vars.TSAxes,'off');
            case 'real'
                set(vars.Line1,'YData',real(ts),'DisplayName','Re(Signal)','Visible','on');
                set(vars.Line2,'Visible','off');
                ylabel(vars.TSAxes,'Re(Signal)');
                legend(vars.TSAxes,'off');
            case 'imag'
                set(vars.Line1,'YData',imag(ts),'DisplayName','Im(Signal)','Visible','on');
                set(vars.Line2,'Visible','off');
                ylabel(vars.TSAxes,'Im(Signal)');
                legend(vars.TSAxes,'off');
            otherwise % 'realimag'
                set(vars.Line1,'YData',real(ts),'DisplayName','Re','Visible','on');
                set(vars.Line2,'YData',imag(ts),'DisplayName','Im','Visible','on');
                ylabel(vars.TSAxes,'Signal');
                legend(vars.TSAxes,'show');
        end
    end

    set(vars.Marker,'XData',colDisp,'YData',rowDisp);
    title(vars.TSAxes, sprintf('Timeseries at data pixel (row=%d, col=%d)',row,col));
    title(vars.RIAxes, sprintf('Re vs Im at data pixel (row=%d, col=%d)',row,col));
    ylim(vars.TSAxes,'auto');
    axis(vars.RIAxes,'equal');
    xlim(vars.RIAxes,[-1 1]);
    ylim(vars.RIAxes,[-1 1]);
end
