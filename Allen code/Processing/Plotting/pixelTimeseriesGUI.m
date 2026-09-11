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
%   pixelTimeseriesGUI({Data1, Data2, ...}, dispImg) overlays multiple
%   datasets sharing the same [nz x nx] pixel grid (e.g. raw data vs. a
%   fit) -- each is plotted as its own line/trajectory, color-coded and
%   labeled in a legend. All datasets are indexed at the same (row,col).
%
%   pixelTimeseriesGUI(Data, dispImg) displays dispImg, a separate
%   [nz x nx] image, instead of the default projection. dispImg may
%   have a different pixel grid size than Data (e.g. a higher-resolution
%   B-mode background behind coarser flow/Doppler data) -- clicks on
%   dispImg are mapped proportionally onto Data's grid. Pass [] to use
%   the default projection while still supplying other options below.
%
%   pixelTimeseriesGUI(...,'Time',t) sets the x-axis for the timeseries
%   plot (default 1:nt for each dataset). For multiple datasets with
%   different numbers of samples, pass a cell array of time vectors
%   matching Data, e.g. {t1, t2}; a single vector is reused for all
%   datasets.
%
%   pixelTimeseriesGUI(...,'DataNames',names) labels each dataset in the
%   legend, e.g. {'Data','Fit'} (default 'Data1','Data2',...).
%
%   pixelTimeseriesGUI(...,'ComplexMode',mode) controls how complex
%   Data is plotted: 'realimag' (default, plots Re and Im as separate
%   lines), 'abs', 'angle', 'real', or 'imag'. Ignored for real Data.
%
%   Additional Name/Value pairs: 'CLim', 'Colormap'.

if nargin<2
    dispImg = [];
end

if ~iscell(Data)
    Data = {Data};
end
nSets = numel(Data);

p = inputParser;
addParameter(p,'Time',[]);
addParameter(p,'DataNames',[]);
addParameter(p,'CLim',[]);
addParameter(p,'Colormap','gray');
addParameter(p,'ComplexMode','realimag');
parse(p,varargin{:});

[nz,nx,~] = size(Data{1});

t = p.Results.Time;
if isempty(t)
    t = arrayfun(@(i) 1:size(Data{i},3), 1:nSets, 'UniformOutput',false);
elseif ~iscell(t)
    t = repmat({t}, 1, nSets);
end

names = p.Results.DataNames;
if isempty(names)
    names = arrayfun(@(i) sprintf('Data%d',i), 1:nSets, 'UniformOutput',false);
end

if isempty(dispImg)
    dispImg = std(double(Data{1}),0,3);
end
[nzDisp,nxDisp] = size(dispImg);

colors = lines(nSets);

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
line1 = gobjects(1,nSets);
line2 = gobjects(1,nSets);
for i = 1:nSets
    line1(i) = plot(tsAx, t{i}, NaN(size(t{i})), 'Color',colors(i,:), ...
        'DisplayName',names{i});
    line2(i) = plot(tsAx, t{i}, NaN(size(t{i})), '--', 'Color',colors(i,:), ...
        'DisplayName',[names{i} ' (Im)'], 'Visible','off');
end
xlabel(tsAx,'Time');
ylabel(tsAx,'Signal');
title(tsAx,'Click a pixel on the image');
grid(tsAx,'on');

riAx = subplot(2,2,4,'Parent',fig);
hold(riAx,'on');
riLine = gobjects(1,nSets);
for i = 1:nSets
    riLine(i) = plot(riAx, NaN, NaN, '.-', 'Color',colors(i,:), 'DisplayName',names{i});
end
axis(riAx,'equal');
xlim(riAx,[-1 1]);
ylim(riAx,[-1 1]);
xlabel(riAx,'Re(Signal)');
ylabel(riAx,'Im(Signal)');
title(riAx,'Re vs Im (complex data only)');
grid(riAx,'on');

vars = struct('Data',{Data},'Time',{t},'Names',{names},'nSets',nSets,...
    'ImgAxes',imgAx,'TSAxes',tsAx,'RIAxes',riAx,...
    'Marker',markerH,'Line1',line1,'Line2',line2,'RILine',riLine,'Fig',fig,...
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

    anyComplex = false;
    for i = 1:vars.nSets
        ts = squeeze(vars.Data{i}(row,col,:));

        if isreal(ts)
            set(vars.Line1(i),'YData',ts,'DisplayName',vars.Names{i},'Visible','on');
            set(vars.Line2(i),'Visible','off');
            set(vars.RILine(i),'XData',NaN,'YData',NaN);
        else
            anyComplex = true;
            set(vars.RILine(i),'XData',real(ts),'YData',imag(ts));
            switch vars.ComplexMode
                case 'abs'
                    set(vars.Line1(i),'YData',abs(ts),...
                        'DisplayName',vars.Names{i},'Visible','on');
                    set(vars.Line2(i),'Visible','off');
                case 'angle'
                    set(vars.Line1(i),'YData',angle(ts),...
                        'DisplayName',vars.Names{i},'Visible','on');
                    set(vars.Line2(i),'Visible','off');
                case 'real'
                    set(vars.Line1(i),'YData',real(ts),...
                        'DisplayName',vars.Names{i},'Visible','on');
                    set(vars.Line2(i),'Visible','off');
                case 'imag'
                    set(vars.Line1(i),'YData',imag(ts),...
                        'DisplayName',vars.Names{i},'Visible','on');
                    set(vars.Line2(i),'Visible','off');
                otherwise % 'realimag'
                    set(vars.Line1(i),'YData',real(ts),...
                        'DisplayName',[vars.Names{i} ' (Re)'],'Visible','on');
                    set(vars.Line2(i),'YData',imag(ts),...
                        'DisplayName',[vars.Names{i} ' (Im)'],'Visible','on');
            end
        end
    end

    switch vars.ComplexMode
        case 'abs',   ylabel(vars.TSAxes,'|Signal|');
        case 'angle', ylabel(vars.TSAxes,'Phase (rad)');
        case 'real',  ylabel(vars.TSAxes,'Re(Signal)');
        case 'imag',  ylabel(vars.TSAxes,'Im(Signal)');
        otherwise,    ylabel(vars.TSAxes,'Signal');
    end

    if vars.nSets>1 || (anyComplex && strcmp(vars.ComplexMode,'realimag'))
        legend(vars.TSAxes,'show');
    else
        legend(vars.TSAxes,'off');
    end
    if vars.nSets>1
        legend(vars.RIAxes,'show');
    end

    set(vars.Marker,'XData',colDisp,'YData',rowDisp);
    title(vars.TSAxes, sprintf('Timeseries at data pixel (row=%d, col=%d)',row,col));
    title(vars.RIAxes, sprintf('Re vs Im at data pixel (row=%d, col=%d)',row,col));
    ylim(vars.TSAxes,'auto');
    axis(vars.RIAxes,'equal');
    xlim(vars.RIAxes,[-1 1]);
    ylim(vars.RIAxes,[-1 1]);
end
