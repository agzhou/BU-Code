function voxelTimeseriesGUI(Data, dispVol, varargin)
% VOXELTIMESERIESGUI  Pick a voxel from orthogonal slice views to plot its timeseries.
%
%   Volumetric counterpart to pixelTimeseriesGUI. Data is a
%   [nz x nx x ny x nt] array (or a cell array {Data1, Data2, ...} of
%   such arrays sharing the same [nz x nx x ny] grid, e.g. raw data vs.
%   a fit), i.e. Data(z,x,y,:) is the timeseries at voxel (z,x,y).
%
%   voxelTimeseriesGUI(Data) opens a figure with three orthogonal slice
%   views through the current voxel (Z-X at fixed Y, Z-Y at fixed X,
%   X-Y at fixed Z), the voxel's timeseries, and (for complex Data) a
%   Re(signal) vs Im(signal) trajectory plot. Click (and drag) on any
%   of the three slice views to move the voxel -- the other two views
%   update to slice through the new location. The displayed volume
%   defaults to the std-over-time projection of Data.
%
%   voxelTimeseriesGUI(Data, dispVol) displays dispVol, a separate
%   [nz x nx x ny] volume, instead of the default projection (must
%   match Data's spatial size). Pass [] to use the default projection
%   while still supplying other options below.
%
%   voxelTimeseriesGUI(...,'Time',t) sets the x-axis for the timeseries
%   plot (default 1:nt for each dataset). For multiple datasets with
%   different numbers of samples, pass a cell array of time vectors
%   matching Data; a single vector is reused for all datasets.
%
%   voxelTimeseriesGUI(...,'DataNames',names) labels each dataset in
%   the legend, e.g. {'Data','Fit'} (default 'Data1','Data2',...).
%
%   voxelTimeseriesGUI(...,'ComplexMode',mode) controls how complex
%   Data is plotted: 'realimag' (default, plots Re and Im as separate
%   lines), 'abs', 'angle', 'real', or 'imag'. Ignored for real Data.
%
%   Additional Name/Value pairs: 'CLim', 'Colormap'.

if nargin<2
    dispVol = [];
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

[nz,nx,ny,~] = size(Data{1});

t = p.Results.Time;
if isempty(t)
    t = arrayfun(@(i) 1:size(Data{i},4), 1:nSets, 'UniformOutput',false);
elseif ~iscell(t)
    t = repmat({t}, 1, nSets);
end

names = p.Results.DataNames;
if isempty(names)
    names = arrayfun(@(i) sprintf('Data%d',i), 1:nSets, 'UniformOutput',false);
end

if isempty(dispVol)
    dispVol = std(double(Data{1}),0,4);
end
if ~isequal(size(dispVol),[nz,nx,ny])
    error('voxelTimeseriesGUI:sizeMismatch',...
        'dispVol is [%d x %d x %d] but must match Data''s spatial size [%d x %d x %d].',...
        size(dispVol,1),size(dispVol,2),size(dispVol,3),nz,nx,ny);
end

colors = lines(nSets);

% initial voxel = volume center
z0 = max(1,round(nz/2));
x0 = max(1,round(nx/2));
y0 = max(1,round(ny/2));

fig = figure('Name','Voxel Timeseries Viewer','NumberTitle','off');

% View 1: Z-X slice at fixed Y
ax1 = subplot(2,3,1,'Parent',fig);
im1 = imagesc(ax1, squeeze(dispVol(:,:,y0)));
axis(ax1,'image'); colormap(ax1,p.Results.Colormap); colorbar(ax1);
if ~isempty(p.Results.CLim), set(ax1,'CLim',p.Results.CLim); end
hold(ax1,'on');
mk1 = plot(ax1, NaN, NaN, 'r+', 'MarkerSize',12,'LineWidth',2);
xlabel(ax1,'X (col)'); ylabel(ax1,'Z (row)');

% View 2: Z-Y slice at fixed X
ax2 = subplot(2,3,2,'Parent',fig);
im2 = imagesc(ax2, squeeze(dispVol(:,x0,:)));
axis(ax2,'image'); colormap(ax2,p.Results.Colormap); colorbar(ax2);
if ~isempty(p.Results.CLim), set(ax2,'CLim',p.Results.CLim); end
hold(ax2,'on');
mk2 = plot(ax2, NaN, NaN, 'r+', 'MarkerSize',12,'LineWidth',2);
xlabel(ax2,'Y (col)'); ylabel(ax2,'Z (row)');

% View 3: X-Y slice at fixed Z
ax3 = subplot(2,3,3,'Parent',fig);
im3 = imagesc(ax3, squeeze(dispVol(z0,:,:)));
axis(ax3,'image'); colormap(ax3,p.Results.Colormap); colorbar(ax3);
if ~isempty(p.Results.CLim), set(ax3,'CLim',p.Results.CLim); end
hold(ax3,'on');
mk3 = plot(ax3, NaN, NaN, 'r+', 'MarkerSize',12,'LineWidth',2);
xlabel(ax3,'Y (col)'); ylabel(ax3,'X (row)');

tsAx = subplot(2,3,[4 5],'Parent',fig);
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
title(tsAx,'Click a voxel on any slice view');
grid(tsAx,'on');

riAx = subplot(2,3,6,'Parent',fig);
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
    'DispVol',dispVol,'VolSize',[nz,nx,ny],...
    'Ax1',ax1,'Ax2',ax2,'Ax3',ax3,'Im1',im1,'Im2',im2,'Im3',im3,...
    'Mk1',mk1,'Mk2',mk2,'Mk3',mk3,...
    'TSAxes',tsAx,'RIAxes',riAx,'Line1',line1,'Line2',line2,'RILine',riLine,...
    'Fig',fig,'ComplexMode',p.Results.ComplexMode,...
    'Voxel',[z0,x0,y0]);

set(im1,'ButtonDownFcn',{@sliceClick,1});
set(im2,'ButtonDownFcn',{@sliceClick,2});
set(im3,'ButtonDownFcn',{@sliceClick,3});
set(fig,'WindowButtonUpFcn',@(~,~) set(fig,'WindowButtonMotionFcn',''));

guidata(fig, vars);
updateVoxel(fig);

end

% Callback subfunctions to support UI actions
function sliceClick(src,~,viewID)
    fig = ancestor(src,'figure');
    vars = guidata(fig);
    vars = pickVoxelFromView(vars, viewID);
    guidata(fig, vars);
    updateVoxel(fig);
    set(fig,'WindowButtonMotionFcn',{@dragMotion,viewID});
end

function dragMotion(fig,~,viewID)
    vars = guidata(fig);
    vars = pickVoxelFromView(vars, viewID);
    guidata(fig, vars);
    updateVoxel(fig);
end

function vars = pickVoxelFromView(vars, viewID)
    % Reads the clicked point on the given slice axes and updates the
    % two voxel coordinates that axes' image spans (row = y-axis of
    % the image, col = x-axis).
    nz = vars.VolSize(1); nx = vars.VolSize(2); ny = vars.VolSize(3);
    voxel = vars.Voxel;

    switch viewID
        case 1
            ax = vars.Ax1;
        case 2
            ax = vars.Ax2;
        case 3
            ax = vars.Ax3;
    end
    pt = get(ax,'CurrentPoint');
    c = round(pt(1,1));
    r = round(pt(1,2));

    switch viewID
        case 1 % Z-X slice: col=X, row=Z
            if c>=1 && c<=nx && r>=1 && r<=nz
                voxel(2) = c; voxel(1) = r;
            end
        case 2 % Z-Y slice: col=Y, row=Z
            if c>=1 && c<=ny && r>=1 && r<=nz
                voxel(3) = c; voxel(1) = r;
            end
        case 3 % X-Y slice: col=Y, row=X
            if c>=1 && c<=ny && r>=1 && r<=nx
                voxel(3) = c; voxel(2) = r;
            end
    end
    vars.Voxel = voxel;
end

function updateVoxel(fig)
    vars = guidata(fig);
    z0 = vars.Voxel(1); x0 = vars.Voxel(2); y0 = vars.Voxel(3);

    set(vars.Im1,'CData', squeeze(vars.DispVol(:,:,y0)));
    set(vars.Im2,'CData', squeeze(vars.DispVol(:,x0,:)));
    set(vars.Im3,'CData', squeeze(vars.DispVol(z0,:,:)));

    set(vars.Mk1,'XData',x0,'YData',z0);
    set(vars.Mk2,'XData',y0,'YData',z0);
    set(vars.Mk3,'XData',y0,'YData',x0);

    title(vars.Ax1, sprintf('Z-X slice (Y=%d)',y0));
    title(vars.Ax2, sprintf('Z-Y slice (X=%d)',x0));
    title(vars.Ax3, sprintf('X-Y slice (Z=%d)',z0));

    anyComplex = false;
    for i = 1:vars.nSets
        ts = squeeze(vars.Data{i}(z0,x0,y0,:));

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

    title(vars.TSAxes, sprintf('Timeseries at voxel (Z=%d, X=%d, Y=%d)',z0,x0,y0));
    title(vars.RIAxes, sprintf('Re vs Im at voxel (Z=%d, X=%d, Y=%d)',z0,x0,y0));
    ylim(vars.TSAxes,'auto');
    axis(vars.RIAxes,'equal');
    xlim(vars.RIAxes,[-1 1]);
    ylim(vars.RIAxes,[-1 1]);
end
