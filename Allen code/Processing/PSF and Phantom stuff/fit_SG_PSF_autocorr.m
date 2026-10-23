function [sigma, info] = fit_SG_PSF_autocorr(varargin)
%% Description:
%   Single-Gaussian (SG) analog of fit_RCA_PSF_autocorr.m: fit the two sigmas of the single-Gaussian g1,s model
%   (sigma_t = sigma_x = sigma_y, and sigma_z) to the MEASURED AUTOCORRELATION of a point-scatterer volume. The
%   output plugs directly into the sigma of vUS_3D_quad_vec.m / vUS_3D_quad_complex_Jac.m / vUS_3D_num_wrapper.m
%   (sigma = [sigma_t, sigma_t, sigma_z]).
%
%   The measured quantity is the same as in fit_RCA_PSF_autocorr.m: for scatterers displaced by d, the
%   single-scatterer field autocorrelation is the normalized autocorrelation of the complex PSF h,
%
%       g1,s(d) = A(d)/A(0),     A(d) = int conj(h(r)) h(r + d) dr,     A(0) = int |h|^2 dr,
%
%   and only |g1,s| is fit (the carrier exp(i 2 k0 v_z tau f) is carried explicitly by the vUS model).
%
%   Model (the single-scatterer term inside the Tangelder-profile integral of the vUS model, for a displacement
%   d = (dx, dy, dz) = f v tau):
%
%       |g1,s| = Clat(r) * Cax(dz),    r = sqrt(dx^2 + dy^2)
%       Clat   = exp(-r^2 / (4 sigma_t^2))
%       Cax    = exp(-dz^2 / (4 sigma_z^2))
%
%   (equivalent to exp(-dx^2/4sx^2 - dy^2/4sy^2 - dz^2/4sz^2) with sx = sy = sigma_t.)
%
%   THE LATERAL PSF OF AN RCA IS NOT ISOTROPIC. Its autocorrelation is a cross: along the x and y axes it has a
%   narrow core plus a long shoulder, and along the diagonals it decays much faster. A single isotropic Gaussian
%   cannot match all directions, so the answer depends on WHICH lateral lags are fit. 'LateralFit' chooses:
%       'map'      (default) all lateral lags in the window, equally weighted: a compromise over directions.
%                  Because the quadrant has more off-axis lags than on-axis ones, this leans toward the faster
%                  (diagonal) decay.
%       'axes'     only lags along the x and y axes (flow along an array axis, where the shoulder matters).
%       'diagonal' only lags along the xy diagonal (flow at 45 deg).
%   info.sigmaByMode gives sigma_t for all three, and info.rmsByDirection the residual of the chosen sigma along
%   each direction, so the spread (the model error for flow direction) is visible. If the spread is large, the
%   single-Gaussian model has an irreducible direction-dependent error and fit_RCA_PSF_autocorr.m is the better
%   description.
%
%   Procedure (same as fit_RCA_PSF_autocorr.m, with one lateral parameter):
%     1. Autocorrelation by FFT with zero padding: the lateral map at dz = 0 is the sum over z of each slice's 2D
%        autocorrelation, the axial curve at dx = dy = 0 the sum over (x, y) of each column's 1D autocorrelation.
%        Exact at integer-voxel lags; no interpolation.
%     2. Normalize by zero lag. 3. Check separability directly at a few combined lags.
%     4. Window in two stages: a pilot fit in a wide window (6x lateral, 4x axial the measured 1/e lag), then a
%        final window of WindowSigmas x the pilot sigma (default 3.5).
%     5. 1D least-squares fits for sigma_t and sigma_z.
%
% Usage:
%   fit_SG_PSF_autocorr                          % interactive: pick the folder with PSF.mat and PData.mat
%   sigma = fit_SG_PSF_autocorr(datapath)         % folder with PSF.mat, PData.mat; params.mat there or one up
%   sigma = fit_SG_PSF_autocorr(IQ, voxelSize)    % IQ [x, y, z] complex, voxelSize [dx, dy, dz] in meters
%   [sigma, info] = fit_SG_PSF_autocorr(..., 'LateralFit', 'axes', 'CompareSigma', [368.1 368.1 126.5]*1e-6)
%
% Inputs:
%   IQ: [nx, ny, nz] complex image of a single point scatterer; z is axial.
%   voxelSize: [dx, dy, dz] [m] or a scalar; PData.PDelta(1:3) .* P.wl for Verasonics data.
%   datapath: folder with PSF.mat (IQ), PData.mat (PData), and params.mat (P) next to it or one folder up.
%
% Name-value options:
%   'LateralFit': 'map' (default), 'axes' or 'diagonal'; see above.
%   'Plot': draw data vs model (default true only when called with no arguments)
%   'CompareSigma': [sigma_t, sigma_z], or [sigma_x, sigma_y, sigma_z] (sigma_x is used laterally), e.g. the
%       sigma currently in use. Reported and plotted against the same measured autocorrelation.
%   'WindowSigmas': final fit window half-width in units of sigma (default 3.5)
%   'MaxLateralLag', 'MaxAxialLag': override the final window [m]
%   'CropLateral', 'CropAxial': keep only +-this distance [m] around the brightest voxel before the
%       autocorrelation. The result depends on how much of the cross-artifact arms the volume contains (see
%       fit_RCA_PSF_autocorr.m); use this to see the spread.
%   'FitAmplitude': also fit an overall amplitude A0 <= 1 (default false), for NOISY measured point scatterers;
%       white noise only adds energy at zero lag and scales every other lag by signal/(signal + noise).
%   'CheckSeparability': direct check at a few combined lags (default true)
%   'Verbose': print a summary (default true)
%
% Outputs:
%   sigma: [sigma_t, sigma_t, sigma_z] [m], FIELD (amplitude) sigmas in the sigma convention of the vUS models
%   info: struct with the window, residual statistics, per-direction results, curves, separability, energy
%       fraction, and CompareSigma statistics.

    %% Parse the inputs
    interactive = false;
    if nargin == 0
        datapath = uigetdir('G:\', 'Select the PSF data path');
        if isequal(datapath, 0), sigma = []; info = []; return, end
        [IQ, vox] = loadPSFFolder(datapath);
        rest = {}; interactive = true;
    elseif ischar(varargin{1}) || isstring(varargin{1})
        [IQ, vox] = loadPSFFolder(char(varargin{1}));
        rest = varargin(2:end);
    else
        if nargin < 2, error('fit_SG_PSF_autocorr:voxelSize', 'Provide the voxel size as the second input.'), end
        IQ = varargin{1}; vox = varargin{2}; rest = varargin(3:end);
    end
    vox = double(vox(:)).';
    if isscalar(vox), vox = vox*[1 1 1]; end
    if ndims(IQ) ~= 3 || numel(vox) ~= 3 || any(vox <= 0)
        error('fit_SG_PSF_autocorr:input', 'IQ must be 3D [x, y, z] and voxelSize [dx, dy, dz] positive.')
    end
    p = inputParser;
    p.addParameter('LateralFit', 'map', @(x) any(strcmpi(x, {'map', 'axes', 'diagonal'})));
    p.addParameter('Plot', interactive, @(x) isscalar(x));
    p.addParameter('CompareSigma', [], @(x) isempty(x) || (isnumeric(x) && any(numel(x) == [2 3])));
    p.addParameter('WindowSigmas', 3.5, @(x) isscalar(x) && x > 0);
    p.addParameter('MaxLateralLag', [], @(x) isempty(x) || (isscalar(x) && x > 0));
    p.addParameter('MaxAxialLag', [], @(x) isempty(x) || (isscalar(x) && x > 0));
    p.addParameter('CropLateral', [], @(x) isempty(x) || (isscalar(x) && x > 0));
    p.addParameter('CropAxial', [], @(x) isempty(x) || (isscalar(x) && x > 0));
    p.addParameter('FitAmplitude', false, @(x) isscalar(x));
    p.addParameter('CheckSeparability', true, @(x) isscalar(x));
    p.addParameter('Verbose', true, @(x) isscalar(x));
    p.parse(rest{:});
    opt = p.Results; mode = lower(opt.LateralFit);
    if ~isempty(opt.CropLateral) || ~isempty(opt.CropAxial)
        [~, i0] = max(abs(IQ(:)));
        [px0, py0, pz0] = ind2sub(size(IQ), i0);
        if ~isempty(opt.CropLateral)
            wx0 = round(opt.CropLateral/vox(1)); wy0 = round(opt.CropLateral/vox(2));
            IQ = IQ(max(1, px0-wx0):min(end, px0+wx0), max(1, py0-wy0):min(end, py0+wy0), :);
        end
        if ~isempty(opt.CropAxial)
            wz0 = round(opt.CropAxial/vox(3));
            IQ = IQ(:, :, max(1, pz0-wz0):min(end, pz0+wz0));
        end
    end
    [nx, ny, nz] = size(IQ);

    %% 1. Autocorrelation at integer lags
    [Rlat, Rax, E0] = psfAutocorr(IQ);
    Cax = abs(Rax(1:nz)).' / E0;                        % [1, nz], a row like zlag
    zlag = (0:nz-1) * vox(3);

    %% 2. Pilot window from the measured 1/e decay lags
    Rlat_all = min((nx - 1)*vox(1), (ny - 1)*vox(2));
    [DX, DY, C] = lateralMap(Rlat, vox, Rlat_all, E0);
    rl = hypot(DX, DY);
    de_lat = min(rl(C < exp(-1) & rl > 0));              % shortest 1/e lag over all lateral directions
    if isempty(de_lat), de_lat = Rlat_all; end
    ia = find(Cax < exp(-1), 1);
    if isempty(ia), de_ax = zlag(end); else, de_ax = zlag(ia); end
    Rlat1 = min(6*de_lat, Rlat_all); Rax1 = min(4*de_ax, zlag(end));
    st1 = fitLateral(DX, DY, C, Rlat1, de_lat, vox, opt.FitAmplitude, mode);
    sa1 = fitAxial(zlag, Cax, Rax1, de_ax, vox(3), opt.FitAmplitude);

    %% 3. Final window and fit
    if isempty(opt.MaxLateralLag), Rlat2 = opt.WindowSigmas * st1; else, Rlat2 = opt.MaxLateralLag; end
    if isempty(opt.MaxAxialLag),   Rax2  = opt.WindowSigmas * sa1; else, Rax2  = opt.MaxAxialLag;   end
    Rlat2 = min(max(Rlat2, 6*max(vox(1:2))), Rlat_all);
    Rax2  = min(max(Rax2,  6*vox(3)), zlag(end));
    [st, A0lat] = fitLateral(DX, DY, C, Rlat2, de_lat, vox, opt.FitAmplitude, mode);
    [sz, A0ax]  = fitAxial(zlag, Cax, Rax2, de_ax, vox(3), opt.FitAmplitude);
    sigma = [st, st, sz];

    %% 4. Diagnostics
    mLat = rl <= Rlat2 & rl > 0;
    resLat = A0lat*exp(-rl.^2/(4*st^2)) - C;
    mAx = zlag <= Rax2 & zlag > 0;
    resAx = A0ax*exp(-zlag.^2/(4*sz^2)) - Cax;
    dirMask = directionMasks(DX, DY, vox, mLat);
    info.sigma = sigma;
    info.sigma_um = sigma * 1e6;
    info.voxelSize = vox;
    info.lateralFit = mode;
    info.window = struct('lateral', Rlat2, 'axial', Rax2, 'pilotLateral', Rlat1, 'pilotAxial', Rax1, ...
        'pilotSigma', [st1, sa1], 'nLateralLags', nnz(mLat), 'nAxialLags', nnz(mAx));
    info.rmsLateral = sqrt(mean(resLat(mLat).^2));
    info.maxLateral = max(abs(resLat(mLat)));
    info.rmsAxial = sqrt(mean(resAx(mAx).^2));
    info.maxAxial = max(abs(resAx(mAx)));
    info.rmsByDirection = struct('axes', sqrt(mean(resLat(dirMask.axes).^2)), 'diagonal', sqrt(mean(resLat(dirMask.diagonal).^2)), ...
        'other', sqrt(mean(resLat(dirMask.other).^2)));
    info.amplitude = [A0lat, A0ax];
    info.measured = struct('DX', DX, 'DY', DY, 'lateral', C, 'zlag', zlag, 'axial', Cax);

    % sigma_t that each lateral-lag choice would give, in the same window
    modes = {'map', 'axes', 'diagonal'}; sbm = zeros(1, 3);
    for k = 1:3, sbm(k) = fitLateral(DX, DY, C, Rlat2, de_lat, vox, opt.FitAmplitude, modes{k}); end
    info.sigmaByMode = cell2struct(num2cell(sbm), modes, 2);

    % Energy of the PSF inside the fit window, around the brightest voxel
    [~, imax] = max(abs(IQ(:)));
    [px, py, pz] = ind2sub([nx, ny, nz], imax);
    wx = ceil(Rlat2/vox(1)); wy = ceil(Rlat2/vox(2)); wz = ceil(Rax2/vox(3));
    Ew = sum(abs(IQ(max(1, px-wx):min(nx, px+wx), max(1, py-wy):min(ny, py+wy), max(1, pz-wz):min(nz, pz+wz))).^2, 'all');
    info.energyFraction = Ew / E0;

    info.undersampled = [vox(1) > st, vox(2) > st, vox(3) > sz];
    if any(info.undersampled)
        warning('fit_SG_PSF_autocorr:undersampled', ...
            ['Voxel size [%.1f %.1f %.1f] um is larger than the fitted sigma_t / sigma_z (%.1f / %.1f um). ' ...
             'The integer-lag autocorrelation is not trustworthy; reconstruct the PSF on a finer grid.'], ...
             vox(1)*1e6, vox(2)*1e6, vox(3)*1e6, st*1e6, sz*1e6)
    end

    % Separability, directly at a few combined lags
    info.separability = [];
    if opt.CheckSeparability
        l1 = max(1, round(st/vox(1))); l2 = max(l1 + 1, round(2*st/vox(1)));
        lz = max(1, round(sz/vox(3)) + [0 1]);
        lagsLat = [l1 0; l2 0; l1 l1]; rows = zeros(0, 6);
        for a = 1:size(lagsLat, 1)
            for zi = 1:numel(lz)
                lx = lagsLat(a, 1); ly = lagsLat(a, 2);
                if lx >= nx || ly >= ny || lz(zi) >= nz, continue, end
                direct = abs(directLag(IQ, lx, ly, lz(zi))) / E0;
                prod_ = abs(Rlat(1 + lx, 1 + ly)) / E0 * Cax(1 + lz(zi));
                rows(end+1, :) = [lx*vox(1), ly*vox(2), lz(zi)*vox(3), direct, prod_, direct - prod_]; %#ok<AGROW>
            end
        end
        info.separability = array2table(rows, 'VariableNames', {'dx_m', 'dy_m', 'dz_m', 'direct', 'product', 'diff'});
        info.maxSeparabilityError = max(abs(rows(:, 6)));
    end

    % Compare against another sigma set over the same window
    info.compare = [];
    if ~isempty(opt.CompareSigma)
        cs = double(opt.CompareSigma(:)).';
        if numel(cs) == 3
            if abs(cs(1) - cs(2)) > 1e-9*cs(1), warning('fit_SG_PSF_autocorr:compare', 'CompareSigma has sigma_x ~= sigma_y; using sigma_x laterally.'), end
            cs = [cs(1), cs(3)];
        end
        cl = A0lat*exp(-rl.^2/(4*cs(1)^2)) - C; ca = A0ax*exp(-zlag.^2/(4*cs(2)^2)) - Cax;
        info.compare = struct('sigma', [cs(1) cs(1) cs(2)], 'rmsLateral', sqrt(mean(cl(mLat).^2)), 'maxLateral', max(abs(cl(mLat))), ...
            'rmsAxial', sqrt(mean(ca(mAx).^2)), 'maxAxial', max(abs(ca(mAx))), ...
            'rmsByDirection', struct('axes', sqrt(mean(cl(dirMask.axes).^2)), 'diagonal', sqrt(mean(cl(dirMask.diagonal).^2))));
    end

    %% 5. Report
    if opt.Verbose
        fprintf('Single-Gaussian autocorrelation fit (voxel %.1f x %.1f x %.1f um, volume %d x %d x %d), lateral lags: %s\n', ...
            vox(1)*1e6, vox(2)*1e6, vox(3)*1e6, nx, ny, nz, mode);
        fprintf('  window: +-%.0f um lateral (%d lags) x +-%.0f um axial (%d lags); %.0f%% of the PSF energy is inside it\n', ...
            1e6*Rlat2, nnz(mLat), 1e6*Rax2, nnz(mAx), 100*info.energyFraction);
        fprintf('Final: sigma_t (= sigma_x = sigma_y) / sigma_z = %.1f / %.1f um (field)\n', info.sigma_um(1), info.sigma_um(3));
        fprintf('  sigma_t if fit to: all lags %.1f | axes only %.1f | diagonal only %.1f um\n', 1e6*sbm);
        fprintf('  lateral RMS %.4f (max %.3f): along axes %.3f, along diagonal %.3f, elsewhere %.3f | axial RMS %.4f (max %.3f)\n', ...
            info.rmsLateral, info.maxLateral, info.rmsByDirection.axes, info.rmsByDirection.diagonal, info.rmsByDirection.other, info.rmsAxial, info.maxAxial);
        if opt.FitAmplitude, fprintf('  fitted amplitude (1 = no noise floor): lateral %.3f, axial %.3f\n', A0lat, A0ax); end
        if ~isempty(info.separability)
            fprintf('  separability: max |direct 3D - lateral*axial| = %.4f over %d lags\n', info.maxSeparabilityError, height(info.separability));
        end
        if ~isempty(info.compare)
            fprintf('  for comparison, sigma = [%.1f %.1f %.1f] um: lateral RMS %.4f (axes %.3f, diagonal %.3f), axial RMS %.4f\n', ...
                1e6*info.compare.sigma, info.compare.rmsLateral, info.compare.rmsByDirection.axes, info.compare.rmsByDirection.diagonal, info.compare.rmsAxial);
        end
    end

    %% 6. Plot
    if opt.Plot
        plotFit(info, sigma, info.compare, Rlat2, Rax2, [A0lat, A0ax], dirMask)
    end
end

%% ----------------------------------------------------------------------------------------------------
function [IQ, vox] = loadPSFFolder(datapath)
    if ~endsWith(datapath, filesep), datapath = [datapath filesep]; end
    S = load([datapath 'PSF.mat'], 'IQ');                    IQ = S.IQ;
    D = load([datapath 'PData.mat'], 'PData');
    if isfile([datapath 'params.mat'])                       % params.mat is next to PSF.mat or one folder up
        Pf = load([datapath 'params.mat'], 'P');
    elseif isfile([datapath '..' filesep 'params.mat'])
        Pf = load([datapath '..' filesep 'params.mat'], 'P');
    else
        error('fit_SG_PSF_autocorr:params', 'params.mat not found in %s or one folder up.', datapath)
    end
    vox = D.PData.PDelta(1:3) * Pf.P.wl;                     % [m]
end

function [Rlat, Rax, E0] = psfAutocorr(IQ)
    % Lateral autocorrelation at dz = 0 and axial autocorrelation at dx = dy = 0, by zero-padded FFT, in chunks.
    [nx, ny, nz] = size(IQ);
    Px = 2*nx; Py = 2*ny; Pz = 2*nz;
    P2 = zeros(Px, Py);
    nzc = max(1, floor(2e7 / (Px*Py)));
    for z0 = 1:nzc:nz
        F = fft2(IQ(:, :, z0:min(nz, z0 + nzc - 1)), Px, Py);
        P2 = P2 + sum(abs(F).^2, 3);
    end
    Rlat = ifft2(P2);
    P1 = zeros(Pz, 1);
    nxc = max(1, floor(2e7 / (ny*Pz)));
    for x0 = 1:nxc:nx
        F = fft(IQ(x0:min(nx, x0 + nxc - 1), :, :), Pz, 3);
        P1 = P1 + squeeze(sum(sum(abs(F).^2, 1), 2));
    end
    Rax = ifft(P1);
    E0 = real(Rlat(1, 1));
end

function [DX, DY, C] = lateralMap(Rlat, vox, Rmax, E0)
    [Px, Py] = size(Rlat);
    ix = 0:min(floor(Rmax/vox(1)), Px/2 - 1);
    iy = 0:min(floor(Rmax/vox(2)), Py/2 - 1);
    [I, J] = ndgrid(ix, iy);
    C1 = abs(Rlat(sub2ind([Px, Py], I + 1, J + 1)));
    C2 = abs(Rlat(sub2ind([Px, Py], I + 1, mod(-J, Py) + 1)));   % lag (i, -j)
    C = (C1 + C2) / (2*E0);
    DX = I*vox(1); DY = J*vox(2);
end

function dm = directionMasks(DX, DY, vox, mLat)
    tol = 0.5*min(vox(1:2));
    dm.axes = mLat & (DX == 0 | DY == 0);
    dm.diagonal = mLat & abs(DX - DY) < tol;
    dm.other = mLat & ~dm.axes & ~dm.diagonal;
end

function [st, A0] = fitLateral(DX, DY, C, R, de, vox, fitAmp, mode)
    r = hypot(DX, DY);
    m = r <= R & r > 0;
    switch mode
        case 'axes',     m = m & (DX == 0 | DY == 0);
        case 'diagonal', m = m & abs(DX - DY) < 0.5*min(vox(1:2));
    end
    if nnz(m) < 4, error('fit_SG_PSF_autocorr:window', 'Only %d lateral lags (%s) in the fit window; use a finer grid or a larger window.', nnz(m), mode), end
    rr = r(m); c = C(m);
    o = optimoptions('lsqnonlin', 'Display', 'off', 'FunctionTolerance', 1e-14, 'StepTolerance', 1e-12);
    lb = [min(vox(1:2))/10, 1]; ub = [max(R, min(vox(1:2))), 1];
    if fitAmp, lb(2) = 0.3; end
    q = lsqnonlin(@(q) q(2)*exp(-rr.^2/(4*q(1)^2)) - c, [max(de/2, min(vox(1:2))/5), 1], lb, ub, o);
    st = q(1); A0 = q(2);
end

function [sa, A0] = fitAxial(zlag, Cax, R, de, dz, fitAmp)
    m = zlag <= R & zlag > 0;
    if nnz(m) < 4, error('fit_SG_PSF_autocorr:window', 'Only %d axial lags in the fit window; use a finer grid or a larger window.', nnz(m)), end
    d = zlag(m); c = Cax(m);
    o = optimoptions('lsqnonlin', 'Display', 'off', 'FunctionTolerance', 1e-14, 'StepTolerance', 1e-12);
    lb = [dz/10, 1]; ub = [max(R, dz), 1];
    if fitAmp, lb(2) = 0.3; end
    q = lsqnonlin(@(q) q(2)*exp(-d.^2/(4*q(1)^2)) - c, [max(de/2, dz/5), 1], lb, ub, o);
    sa = q(1); A0 = q(2);
end

function a = directLag(IQ, lx, ly, lz)
    % sum conj(IQ(r)) * IQ(r + [lx, ly, lz]) over r, in z-chunks to keep memory use modest
    [nx, ny, nz] = size(IQ);
    a = 0;
    for z0 = 1:32:nz - lz
        z1 = min(nz - lz, z0 + 31);
        A = IQ(1:nx-lx, 1:ny-ly, z0:z1);
        B = IQ(1+lx:nx, 1+ly:ny, z0+lz:z1+lz);
        a = a + sum(conj(A) .* B, 'all');
    end
end

function plotFit(info, sigma, cmp, Rlat, Rax, A0, dm)
    st = sigma(1); sz = sigma(3); M = info.measured; um = 1e6; vox = info.voxelSize;
    figure('Name', 'Single-Gaussian PSF autocorrelation fit'); tiledlayout(2, 2, 'TileSpacing', 'compact');
    col = lines(3);

    nexttile; ax = M.DY == 0 & M.DX <= 1.6*Rlat; d = M.DX(ax); dd = linspace(0, max(d), 300);
    h1 = plot(d*um, M.lateral(ax), 'o', 'Color', col(1,:), 'MarkerFaceColor', col(1,:)); hold on
    dg = abs(M.DX - M.DY) < 0.5*min(vox(1:2)) & M.DX*sqrt(2) <= 1.6*Rlat; dgd = M.DX(dg)*sqrt(2);
    h2 = plot(dgd*um, M.lateral(dg), 's', 'Color', col(2,:), 'MarkerFaceColor', col(2,:));
    h3 = plot(dd*um, A0(1)*exp(-dd.^2/(4*st^2)), 'k-', 'LineWidth', 1.8);
    labs = {'measured, along x', 'measured, diagonal', sprintf('single Gaussian, \\sigma_t = %.0f \\mum', st*um)}; hs = [h1 h2 h3];
    if ~isempty(cmp), h4 = plot(dd*um, A0(1)*exp(-dd.^2/(4*cmp.sigma(1)^2)), 'k--', 'LineWidth', 1.4); hs(end+1) = h4; labs{end+1} = sprintf('comparison, \\sigma_t = %.0f \\mum', cmp.sigma(1)*um); end %#ok<AGROW>
    xline(Rlat*um, ':', 'HandleVisibility', 'off'); hold off; grid on; xlabel('lateral lag [\mum]'); ylabel('|g_{1,s}|'); ylim([0 1.02]); title('lateral cuts'); legend(hs, labs, 'Location', 'northeast', 'FontSize', 8)

    nexttile; m = hypot(M.DX, M.DY) <= 1.6*Rlat & ~(M.DX == 0 & M.DY == 0); r = hypot(M.DX(m), M.DY(m)); dd = linspace(0, max(r), 300);
    scatter(r*um, M.lateral(m), 10, 'filled', 'MarkerFaceAlpha', 0.5); hold on
    plot(dd*um, A0(1)*exp(-dd.^2/(4*st^2)), 'k-', 'LineWidth', 1.8); xline(Rlat*um, ':'); hold off; grid on
    xlabel('lateral lag radius [\mum]'); ylabel('|g_{1,s}|'); ylim([0 1.02]); title('all lateral lags (spread = anisotropy a single Gaussian cannot follow)')

    nexttile; z = M.zlag(M.zlag <= 1.6*Rax); c = M.axial(1:numel(z)); dd = linspace(0, max(z), 300);
    plot(z*um, c, 'ko', 'MarkerFaceColor', 'k'); hold on
    plot(dd*um, A0(2)*exp(-dd.^2/(4*sz^2)), 'LineWidth', 1.8)
    if ~isempty(cmp), plot(dd*um, A0(2)*exp(-dd.^2/(4*cmp.sigma(3)^2)), '--', 'LineWidth', 1.4); end
    xline(Rax*um, ':'); hold off; grid on; xlabel('axial lag [\mum]'); ylim([0 1.02]); title(sprintf('axial, \\sigma_z = %.0f \\mum', sz*um))

    nexttile; m = hypot(M.DX, M.DY) <= Rlat; R = A0(1)*exp(-(M.DX.^2 + M.DY.^2)/(4*st^2)) - M.lateral; R(~m) = NaN;
    imagesc(M.DX(:, 1)*um, M.DY(1, :)*um, R.'); axis image; set(gca, 'YDir', 'normal'); colorbar; clim(max(abs(R(:)), [], 'omitnan')*[-1 1]);
    xlim([0 Rlat]*um); ylim([0 Rlat]*um); xlabel('dx [\mum]'); ylabel('dy [\mum]'); title('lateral residual (model - measured)')
end
