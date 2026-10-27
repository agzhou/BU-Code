function voxelSpectrogramGUI(S, dispVol, varargin)
% VOXELSPECTROGRAMGUI  Pick a voxel from orthogonal slice views to plot its spectrogram.
%
%   Spectrogram counterpart to voxelTimeseriesGUI. S is a
%   [nz x nx x ny x nWin x nFreq] array of per-window spectra (e.g. the
%   sliding-window FFT built in LFO3D.m), i.e. squeeze(S(z,x,y,:,:)) is
%   the [nWin x nFreq] spectrogram at voxel (z,x,y). S may be complex
%   (raw FFT output) or real (already a magnitude/power).
%
%   voxelSpectrogramGUI(S) opens a figure with three orthogonal slice
%   views through the current voxel (Z-X at fixed Y, Z-Y at fixed X,
%   X-Y at fixed Z), the voxel's spectrogram (time vs frequency), and
%   its window-averaged spectrum. Click (and drag) on any of the three
%   slice views to move the voxel -- the other views update. The
%   displayed volume defaults to the mean spectral magnitude of S.
%
%   voxelSpectrogramGUI(S, dispVol) displays dispVol, a separate
%   [nz x nx x ny] volume (e.g. the power Doppler map), instead of the
%   default (must match S's spatial size). Pass [] to use the default
%   while still supplying other options below.
%
%   voxelSpectrogramGUI(...,'Time',t) sets the time axis of the
%   spectrogram, one value per window (default 1:nWin), e.g. the start
%   or center time of each window [s].
%
%   voxelSpectrogramGUI(...,'Freq',f) sets the frequency axis, one value
%   per spectral bin (default 1:nFreq), e.g. the fftshifted axis [Hz].
%
%   voxelSpectrogramGUI(...,'Mode',mode) sets the initial quantity
%   plotted: 'abs' (default), 'power', 'db' or 'angle'. Can be changed
%   from the Display menu in the figure.
%
%   voxelSpectrogramGUI(...,'FreqRange',[fmin fmax]) limits the
%   displayed frequency range (default: the full axis).
%
%   Controls in the figure: Display (quantity shown), Hide DC (blank the
%   zero-frequency bin, which usually dominates PDI spectra), f >= 0 only
%   (show positive frequencies only), Lock CLim (freeze the spectrogram
%   color limits while moving between voxels, to compare voxels on the
%   same scale).
%
%   Additional Name/Value pairs: 'CLim', 'Colormap' (slice views),
%   'SpecCLim', 'SpecColormap' (spectrogram).

if nargin<2
    dispVol = [];
end

p = inputParser;
addParameter(p,'Time',[]);
addParameter(p,'Freq',[]);
addParameter(p,'Mode','abs');
addParameter(p,'FreqRange',[]);
addParameter(p,'CLim',[]);
addParameter(p,'Colormap','gray');
addParameter(p,'SpecCLim',[]);
addParameter(p,'SpecColormap','turbo');
parse(p,varargin{:});

if ndims(S)~=5 && ndims(S)~=4
    error('voxelSpectrogramGUI:badSize',...
        'S must be [nz x nx x ny x nWin x nFreq]; got %s.',mat2str(size(S)));
end
[nz,nx,ny,nw,nf] = size(S);

t = p.Results.Time;
if isempty(t)
    t = 1:nw;
    tLabel = 'Window';
else
    tLabel = 'Time [s]';
end
t = t(:).';
if numel(t)~=nw
    error('voxelSpectrogramGUI:sizeMismatch',...
        '''Time'' has %d elements but S has %d windows.',numel(t),nw);
end

f = p.Results.Freq;
if isempty(f)
    f = 1:nf;
    fLabel = 'Frequency bin';
else
    fLabel = 'Frequency [Hz]';
end
f = f(:).';
if numel(f)~=nf
    error('voxelSpectrogramGUI:sizeMismatch',...
        '''Freq'' has %d elements but S has %d frequency bins.',numel(f),nf);
end

if isempty(dispVol)
    % mean spectral magnitude, accumulated one window at a time so a large
    % complex S is never duplicated
    dispVol = zeros(nz,nx,ny);
    for wi = 1:nw
        dispVol = dispVol + mean(abs(double(S(:,:,:,wi,:))),5);
    end
    dispVol = dispVol/nw;
end
if ~isequal(size(dispVol),[nz,nx,ny])
    error('voxelSpectrogramGUI:sizeMismatch',...
        'dispVol is [%d x %d x %d] but must match S''s spatial size [%d x %d x %d].',...
        size(dispVol,1),size(dispVol,2),size(dispVol,3),nz,nx,ny);
end

modes = {'abs','power','db','angle'};
modeNames = {'|S|','Power |S|^2','dB','Phase'};
modeIdx = find(strcmpi(p.Results.Mode,modes),1);
if isempty(modeIdx)
    error('voxelSpectrogramGUI:badMode','''Mode'' must be one of: %s.',strjoin(modes,', '));
end

[~,dcIdx] = min(abs(f)); % bin closest to zero frequency

% initial voxel = volume center
z0 = max(1,round(nz/2));
x0 = max(1,round(nx/2));
y0 = max(1,round(ny/2));

fig = figure('Name','Voxel Spectrogram Viewer','NumberTitle','off');

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

% Spectrogram (time on x, frequency on y). imagesc only uses the first and
% last time/frequency values, so a single window needs a nonzero extent.
tExt = [t(1), t(end)];
if tExt(1)==tExt(2), tExt = tExt + [-0.5, 0.5]; end
fExt = [f(1), f(end)];
if fExt(1)==fExt(2), fExt = fExt + [-0.5, 0.5]; end
specAx = axes('Parent',fig,'Position',[0.07 0.17 0.50 0.27]);
specIm = imagesc(specAx, tExt, fExt, NaN(nf,nw));
set(specAx,'YDir','normal','Color',[0.85 0.85 0.85]);
colormap(specAx,p.Results.SpecColormap);
specCB = colorbar(specAx);
xlabel(specAx,tLabel);
ylabel(specAx,fLabel);

% Window-averaged spectrum
spAx = axes('Parent',fig,'Position',[0.69 0.17 0.27 0.27]);
hold(spAx,'on');
spLine = plot(spAx, f, NaN(size(f)), 'k');
xlabel(spAx,fLabel);
grid(spAx,'on');

% Controls
ui.ModeLbl = uicontrol(fig,'Style','text','Units','normalized',...
    'Position',[0.02 0.02 0.07 0.04],'String','Display:','HorizontalAlignment','right');
ui.Mode = uicontrol(fig,'Style','popupmenu','Units','normalized',...
    'Position',[0.095 0.025 0.12 0.04],'String',modeNames,'Value',modeIdx);
ui.HideDC = uicontrol(fig,'Style','checkbox','Units','normalized',...
    'Position',[0.25 0.02 0.12 0.04],'String','Hide DC','Value',1);
ui.PosOnly = uicontrol(fig,'Style','checkbox','Units','normalized',...
    'Position',[0.38 0.02 0.14 0.04],'String','f >= 0 only','Value',0);
ui.Lock = uicontrol(fig,'Style','checkbox','Units','normalized',...
    'Position',[0.53 0.02 0.14 0.04],'String','Lock CLim','Value',0);

if ~isempty(p.Results.SpecCLim)
    set(specAx,'CLim',p.Results.SpecCLim);
    set(ui.Lock,'Value',1);
end

vars = struct('S',{S},'Time',t,'Freq',f,'DCIdx',dcIdx,'FreqRange',p.Results.FreqRange,...
    'Modes',{modes},'ModeNames',{modeNames},...
    'DispVol',dispVol,'VolSize',[nz,nx,ny],'nWin',nw,'nFreq',nf,...
    'Ax1',ax1,'Ax2',ax2,'Ax3',ax3,'Im1',im1,'Im2',im2,'Im3',im3,...
    'Mk1',mk1,'Mk2',mk2,'Mk3',mk3,...
    'SpecAxes',specAx,'SpecIm',specIm,'SpecCB',specCB,...
    'SpAxes',spAx,'SpLine',spLine,'UI',ui,...
    'Fig',fig,'Voxel',[z0,x0,y0]);

set(im1,'ButtonDownFcn',{@sliceClick,1});
set(im2,'ButtonDownFcn',{@sliceClick,2});
set(im3,'ButtonDownFcn',{@sliceClick,3});
set(fig,'WindowButtonUpFcn',@(~,~) set(fig,'WindowButtonMotionFcn',''));

set(ui.Mode,'Callback',@(~,~) modeChanged(fig));
set([ui.HideDC, ui.PosOnly, ui.Lock],'Callback',@(~,~) updateVoxel(fig));

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

function modeChanged(fig)
    % A color scale locked for one quantity is meaningless for another
    vars = guidata(fig);
    set(vars.UI.Lock,'Value',0);
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
    nw = vars.nWin; nf = vars.nFreq;
    ui = vars.UI;

    set(vars.Im1,'CData', squeeze(vars.DispVol(:,:,y0)));
    set(vars.Im2,'CData', squeeze(vars.DispVol(:,x0,:)));
    set(vars.Im3,'CData', squeeze(vars.DispVol(z0,:,:)));

    set(vars.Mk1,'XData',x0,'YData',z0);
    set(vars.Mk2,'XData',y0,'YData',z0);
    set(vars.Mk3,'XData',y0,'YData',x0);

    title(vars.Ax1, sprintf('Z-X slice (Y=%d)',y0));
    title(vars.Ax2, sprintf('Z-Y slice (X=%d)',x0));
    title(vars.Ax3, sprintf('X-Y slice (Z=%d)',z0));

    % [nWin x nFreq] spectrogram of the selected voxel
    sv = double(reshape(vars.S(z0,x0,y0,:,:), nw, nf));

    mode = vars.Modes{get(ui.Mode,'Value')};
    switch mode
        case 'abs'
            Z = abs(sv);
            sp = mean(Z,1);
            cbLabel = '|S|';
            spLabel = 'Mean |S|';
        case 'power'
            Z = abs(sv).^2;
            sp = mean(Z,1);
            cbLabel = '|S|^2';
            spLabel = 'Mean |S|^2';
        case 'db'
            Z = 20*log10(abs(sv) + eps);
            sp = 10*log10(mean(abs(sv).^2,1) + eps);
            cbLabel = '|S| [dB]';
            spLabel = 'Mean power [dB]';
        case 'angle'
            Z = angle(sv);
            sp = mean(abs(sv),1); % a mean phase is not meaningful; show magnitude
            cbLabel = 'Phase [rad]';
            spLabel = 'Mean |S|';
    end
    Z = Z.';  % [nFreq x nWin]: frequency on y, time on x

    if get(ui.HideDC,'Value')
        Z(vars.DCIdx,:) = NaN;
        sp(vars.DCIdx) = NaN;
    end

    % NaN (hidden DC) cells are drawn transparent over the gray axes color
    set(vars.SpecIm,'CData',Z,'AlphaData',~isnan(Z));
    set(vars.SpLine,'YData',sp);

    if get(ui.Lock,'Value')
        set(vars.SpecAxes,'CLimMode','manual');
    else
        set(vars.SpecAxes,'CLimMode','auto');
    end
    vars.SpecCB.Label.String = cbLabel;
    ylabel(vars.SpAxes,spLabel);

    % frequency limits (shared by the spectrogram y-axis and spectrum x-axis)
    f = vars.Freq;
    if ~isempty(vars.FreqRange)
        flim = vars.FreqRange;
    else
        flim = [min(f), max(f)];
    end
    if get(ui.PosOnly,'Value')
        flim(1) = max(flim(1),0);
    end
    if flim(2)>flim(1)
        set(vars.SpecAxes,'YLim',flim);
        set(vars.SpAxes,'XLim',flim);
    end
    ylim(vars.SpAxes,'auto');

    title(vars.SpecAxes, sprintf('Spectrogram at voxel (Z=%d, X=%d, Y=%d)',z0,x0,y0));
    title(vars.SpAxes, 'Window-averaged spectrum');
end
