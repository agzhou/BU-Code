function [sigma, info] = fit_2D_PSF_autocorr(varargin)
%% Description:
%   Fit the two sigmas of the 2D (linear array) single-Gaussian g1,s model, sigma_x and sigma_z, to the MEASURED
%   AUTOCORRELATION of the complex PSF of a point scatterer. 2D analog of fit_SG_PSF_autocorr.m, with the same
%   robustness measures, and a function version of fit_2D_PSF_autocorr_test.m. The output plugs directly into the
%   sigma = [sigma_x, sigma_z] of vUS_2D_quad_batchModel.m / vUS_2D_num_wrapper.m / vUS_2D_newmodel.m.
%
%   For scatterers displaced by d = (dx, dz), the single-scatterer field autocorrelation is the normalized
%   autocorrelation of the complex PSF h,
%
%       g1,s(d) = A(d)/A(0),     A(d) = int conj(h(r)) h(r + d) dr,     A(0) = int |h|^2 dr,
%
%   and only |g1,s| is fit (the carrier exp(i 2 k0 v_z tau f) is carried explicitly by the vUS model).
%
%   Model:   |g1,s| = Cx(dx) * Cz(dz),    Cx = exp(-dx^2/(4 sigma_x^2)),    Cz = exp(-dz^2/(4 sigma_z^2))
%
%   Differences from fit_2D_PSF_autocorr_test.m (the script this replaces), and why:
%     - The whole image enters the autocorrelation (zero-padded FFT, exact at integer-pixel lags) and is
%       normalized by the TOTAL energy A(0). The script uses a fixed +-20 pixel window around the scatterer and
%       normalizes by that window's energy, so the answer depends on the window size and the energy outside it
%       (lateral sidelobes) is ignored. info.energyFraction reports how much energy is inside the fit window, and
%       'CropLateral'/'CropAxial' reproduce a restricted-volume result on purpose.
%     - The fit window is derived from the data in two stages (pilot window from the measured 1/e lags, then
%       WindowSigmas x the pilot sigma), not fixed at 10 pixels. The Gaussian is not fit to far lags where the
%       real PSF has sidelobes and tails the model can't represent.
%     - The scatterer position is not needed (shift invariance), so the "abs(PData.Origin)" hack is gone.
%     - Several start points per fit; their agreement and any sigma at a bound are reported.
%     - Separability |g1,s(dx,dz)| = Cx Cz is checked on the full 2D autocorrelation inside the window, not assumed.
%     - Undersampling is detected: if a pixel is larger than the fitted sigma, the integer-lag autocorrelation is
%       not trustworthy (a warning is issued).
%     - 'FitAmplitude' handles noisy measured scatterers: white noise only adds energy at zero lag, so every lag
%       > 0 is scaled by signal/(signal + noise) after normalization, and a fit without an amplitude shrinks sigma.
%
% Usage:
%   fit_2D_PSF_autocorr                         % interactive: pick the folder with PSF.mat and PData.mat
%   sigma = fit_2D_PSF_autocorr(datapath)        % folder with PSF.mat, PData.mat; params.mat there or one up
%   sigma = fit_2D_PSF_autocorr(IQ, voxelSize)   % IQ [nz, nx] complex, voxelSize [dx, dz] in meters
%   [sigma, info] = fit_2D_PSF_autocorr(..., 'FitMode', 'map', 'CompareSigma', [60.55 73.95]*1e-6, 'Plot', true)
%
% Inputs:
%   IQ: [nz, nx] complex image of a single point scatterer (dimension 1 is axial z, dimension 2 lateral x, as in
%       the Verasonics 2D reconstruction and the 2D PSF.mat files).
%   voxelSize: [dx, dz] pixel size [m] (note the order, matching sigma), or a scalar.
%       [PData.PDelta(1), PData.PDelta(3)] .* P.wl for Verasonics data.
%   datapath: folder with PSF.mat (IQ), PData.mat (PData), and params.mat (P) next to it or one folder up.
%
% Name-value options:
%   'FitMode': 'curves' (default) fits sigma_x to |A(dx, 0)| and sigma_z to |A(0, dz)|, as the script does.
%       'map' fits both to the whole 2D map |A(dx, dz)| inside the elliptical window (product model).
%   'Plot': draw data vs model (default true only when called with no arguments)
%   'CompareSigma': [sigma_x, sigma_z] [m], e.g. the sigma currently in use; reported and plotted
%   'WindowSigmas': final fit window half-width in units of sigma (default 3.5)
%   'MaxLateralLag', 'MaxAxialLag': override the final window [m]
%   'CropLateral', 'CropAxial': keep only +-this distance [m] around the brightest pixel first
%   'FitAmplitude': also fit an overall amplitude A0 <= 1 (default false); use for NOISY measured scatterers
%   'Verbose': print a summary (default true)
%
% Outputs:
%   sigma: [sigma_x, sigma_z] [m], FIELD (amplitude) sigmas in the sigma convention of the vUS_2D models
%   info: struct with the window, residual statistics, measured and fitted curves, energy fraction, separability
%       (info.maxSeparabilityError), sampling check, and CompareSigma statistics.
%
% Caveats: one depth and an ideal medium (the PSF of a linear array changes with depth: repeat for the depths you
%   image); the Gaussian form is an approximation (check info.rmsX, info.rmsZ and the plot); the result depends on
%   how much sidelobe energy the image contains (use 'CropLateral' to see the spread).

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
        if nargin < 2, error('fit_2D_PSF_autocorr:voxelSize', 'Provide the pixel size [dx, dz] as the second input.'), end
        IQ = varargin{1}; vox = varargin{2}; rest = varargin(3:end);
    end
    vox = double(vox(:)).';
    if isscalar(vox), vox = vox*[1 1]; end
    if ~ismatrix(IQ) || min(size(IQ)) < 8 || numel(vox) ~= 2 || any(vox <= 0)
        error('fit_2D_PSF_autocorr:input', 'IQ must be a 2D image [nz, nx] and voxelSize [dx, dz] positive.')
    end
    p = inputParser;
    p.addParameter('FitMode', 'curves', @(x) any(strcmpi(x, {'curves', 'map'})));
    p.addParameter('Plot', interactive, @(x) isscalar(x));
    p.addParameter('CompareSigma', [], @(x) isempty(x) || (isnumeric(x) && numel(x) == 2));
    p.addParameter('WindowSigmas', 3.5, @(x) isscalar(x) && x > 0);
    p.addParameter('MaxLateralLag', [], @(x) isempty(x) || (isscalar(x) && x > 0));
    p.addParameter('MaxAxialLag', [], @(x) isempty(x) || (isscalar(x) && x > 0));
    p.addParameter('CropLateral', [], @(x) isempty(x) || (isscalar(x) && x > 0));
    p.addParameter('CropAxial', [], @(x) isempty(x) || (isscalar(x) && x > 0));
    p.addParameter('FitAmplitude', false, @(x) isscalar(x));
    p.addParameter('Verbose', true, @(x) isscalar(x));
    p.parse(rest{:});
    opt = p.Results; mapMode = strcmpi(opt.FitMode, 'map');
    dx = vox(1); dz = vox(2);
    if ~isempty(opt.CropLateral) || ~isempty(opt.CropAxial)
        [~, i0] = max(abs(IQ(:)));
        [pz0, px0] = ind2sub(size(IQ), i0);
        if ~isempty(opt.CropLateral)
            wx0 = round(opt.CropLateral/dx); IQ = IQ(:, max(1, px0-wx0):min(end, px0+wx0));
        end
        if ~isempty(opt.CropAxial)
            wz0 = round(opt.CropAxial/dz); IQ = IQ(max(1, pz0-wz0):min(end, pz0+wz0), :);
        end
    end
    [nz, nx] = size(IQ);

    %% 1. Autocorrelation at integer lags (exact; zero-padded FFT of the whole image)
    R = ifft2(abs(fft2(IQ, 2*nz, 2*nx)).^2);
    E0 = real(R(1, 1));
    zlag = (0:nz-1)*dz; xlag = (0:nx-1)*dx;
    Cz = abs(R(1:nz, 1)).' / E0;                         % |A(0, dz)|, a row like zlag
    Cx = abs(R(1, 1:nx)) / E0;                           % |A(dx, 0)|
    [ZL, XL] = ndgrid(zlag, xlag);
    jneg = mod(-(0:nx-1), 2*nx) + 1;                     % lag (dz, -dx): average with (dz, +dx) (mirror symmetry)
    C2 = (abs(R(1:nz, 1:nx)) + abs(R(1:nz, jneg))) / (2*E0);

    %% 2. Pilot window from the measured 1/e decay lags
    Rx_all = xlag(end); Rz_all = zlag(end);
    ix = find(Cx < exp(-1), 1); if isempty(ix), de_x = Rx_all; else, de_x = xlag(ix); end
    iz = find(Cz < exp(-1), 1); if isempty(iz), de_z = Rz_all; else, de_z = zlag(iz); end
    Rx1 = min(6*de_x, Rx_all); Rz1 = min(4*de_z, Rz_all);
    sx1 = fitCurve(xlag, Cx, Rx1, de_x, dx, opt.FitAmplitude);
    sz1 = fitCurve(zlag, Cz, Rz1, de_z, dz, opt.FitAmplitude);

    %% 3. Final window and fit
    if isempty(opt.MaxLateralLag), Rx2 = opt.WindowSigmas * sx1; else, Rx2 = opt.MaxLateralLag; end
    if isempty(opt.MaxAxialLag),   Rz2 = opt.WindowSigmas * sz1; else, Rz2 = opt.MaxAxialLag;   end
    Rx2 = min(max(Rx2, 6*dx), Rx_all);                   % at least 6 pixels, at most what the image holds
    Rz2 = min(max(Rz2, 6*dz), Rz_all);
    if mapMode
        [sx, sz, A0x, nAgree, atBound] = fitMap(XL, ZL, C2, Rx2, Rz2, [de_x, de_z], vox, opt.FitAmplitude);
        A0z = A0x;
    else
        [sx, A0x, nAgreeX, bX] = fitCurve(xlag, Cx, Rx2, de_x, dx, opt.FitAmplitude);
        [sz, A0z, nAgreeZ, bZ] = fitCurve(zlag, Cz, Rz2, de_z, dz, opt.FitAmplitude);
        nAgree = min(nAgreeX, nAgreeZ); atBound = [bX, bZ];
    end
    sigma = [sx, sz];

    %% 4. Diagnostics
    mX = xlag <= Rx2 & xlag > 0; mZ = zlag <= Rz2 & zlag > 0;
    resX = A0x*exp(-xlag.^2/(4*sx^2)) - Cx;
    resZ = A0z*exp(-zlag.^2/(4*sz^2)) - Cz;
    win = (XL/Rx2).^2 + (ZL/Rz2).^2 <= 1 & ~(XL == 0 & ZL == 0);
    Mmodel = A0x*exp(-XL.^2/(4*sx^2)) .* exp(-ZL.^2/(4*sz^2));
    resMap = Mmodel - C2;
    info.sigma = sigma; info.sigma_um = sigma*1e6; info.voxelSize = vox; info.fitMode = lower(opt.FitMode);
    info.window = struct('lateral', Rx2, 'axial', Rz2, 'pilotLateral', Rx1, 'pilotAxial', Rz1, 'pilotSigma', [sx1, sz1], ...
        'nLateralLags', nnz(mX), 'nAxialLags', nnz(mZ), 'nMapLags', nnz(win));
    info.rmsX = sqrt(mean(resX(mX).^2));  info.maxX = max(abs(resX(mX)));
    info.rmsZ = sqrt(mean(resZ(mZ).^2));  info.maxZ = max(abs(resZ(mZ)));
    info.rmsMap = sqrt(mean(resMap(win).^2)); info.maxMap = max(abs(resMap(win)));
    info.amplitude = [A0x, A0z]; info.startsAgreeing = nAgree; info.atBound = atBound;
    info.measured = struct('xlag', xlag, 'zlag', zlag, 'lateral', Cx, 'axial', Cz, 'map', C2);

    % Separability of the DATA: |A(dx, dz)| vs |A(dx, 0)| * |A(0, dz)| over the window (independent of the fit)
    sepErr = C2 - Cz(:) * Cx;
    info.maxSeparabilityError = max(abs(sepErr(win)));
    info.rmsSeparabilityError = sqrt(mean(sepErr(win).^2));

    % Energy of the PSF inside the fit window, around the brightest pixel
    [~, imax] = max(abs(IQ(:))); [pz, px] = ind2sub([nz, nx], imax);
    wx = ceil(Rx2/dx); wz = ceil(Rz2/dz);
    info.energyFraction = sum(abs(IQ(max(1, pz-wz):min(nz, pz+wz), max(1, px-wx):min(nx, px+wx))).^2, 'all') / E0;

    % Sampling and bounds
    info.undersampled = [dx > sx, dz > sz];
    if any(info.undersampled)
        warning('fit_2D_PSF_autocorr:undersampled', ...
            ['Pixel size [dx dz] = [%.1f %.1f] um is larger than the fitted sigma_x / sigma_z (%.1f / %.1f um). ' ...
             'The integer-lag autocorrelation is not trustworthy; reconstruct the PSF on a finer grid.'], dx*1e6, dz*1e6, sx*1e6, sz*1e6)
    end
    if any(atBound)
        warning('fit_2D_PSF_autocorr:atBound', 'A sigma hit its bound; the fit is not well determined (window or sampling).')
    end

    % Compare against another sigma set over the same window
    info.compare = [];
    if ~isempty(opt.CompareSigma)
        cs = double(opt.CompareSigma(:)).';
        cx = A0x*exp(-xlag.^2/(4*cs(1)^2)) - Cx; cz = A0z*exp(-zlag.^2/(4*cs(2)^2)) - Cz;
        cm = A0x*exp(-XL.^2/(4*cs(1)^2)) .* exp(-ZL.^2/(4*cs(2)^2)) - C2;
        info.compare = struct('sigma', cs, 'rmsX', sqrt(mean(cx(mX).^2)), 'maxX', max(abs(cx(mX))), ...
            'rmsZ', sqrt(mean(cz(mZ).^2)), 'maxZ', max(abs(cz(mZ))), 'rmsMap', sqrt(mean(cm(win).^2)), 'maxMap', max(abs(cm(win))));
    end

    %% 5. Report
    if opt.Verbose
        fprintf('2D single-Gaussian autocorrelation fit (pixel dx = %.1f, dz = %.1f um, image %d x %d), fit mode: %s\n', dx*1e6, dz*1e6, nz, nx, info.fitMode);
        fprintf('  pilot sigma_x / sigma_z = %.1f / %.1f um\n', 1e6*sx1, 1e6*sz1);
        fprintf('  window: +-%.0f um lateral (%d lags) x +-%.0f um axial (%d lags); %.0f%% of the PSF energy is inside it\n', ...
            1e6*Rx2, nnz(mX), 1e6*Rz2, nnz(mZ), 100*info.energyFraction);
        fprintf('Final: sigma_x / sigma_z = %.1f / %.1f um (field)\n', info.sigma_um);
        fprintf('  RMS: lateral %.4f (max %.3f), axial %.4f (max %.3f), 2D map %.4f (max %.3f); %d start points agree\n', ...
            info.rmsX, info.maxX, info.rmsZ, info.maxZ, info.rmsMap, info.maxMap, nAgree);
        fprintf('  separability of the data: max |A(dx,dz) - A(dx,0)A(0,dz)| = %.4f in the window\n', info.maxSeparabilityError);
        if opt.FitAmplitude, fprintf('  fitted amplitude (1 = no noise floor): lateral %.3f, axial %.3f\n', A0x, A0z); end
        if ~isempty(info.compare)
            fprintf('  for comparison, sigma = [%.1f %.1f] um: RMS lateral %.4f, axial %.4f, 2D map %.4f\n', ...
                1e6*info.compare.sigma, info.compare.rmsX, info.compare.rmsZ, info.compare.rmsMap);
        end
    end

    %% 6. Plot
    if opt.Plot
        plotFit(IQ, info, sigma, info.compare, Rx2, Rz2, [A0x, A0z], win, resMap, [pz, px])
    end
end

%% ----------------------------------------------------------------------------------------------------
function [IQ, vox] = loadPSFFolder(datapath)
    if ~endsWith(datapath, filesep), datapath = [datapath filesep]; end
    S = load([datapath 'PSF.mat'], 'IQ');                    IQ = squeeze(S.IQ);
    D = load([datapath 'PData.mat'], 'PData');
    if isfile([datapath 'params.mat'])                       % params.mat is next to PSF.mat or one folder up
        Pf = load([datapath 'params.mat'], 'P');
    elseif isfile([datapath '..' filesep 'params.mat'])
        Pf = load([datapath '..' filesep 'params.mat'], 'P');
    else
        error('fit_2D_PSF_autocorr:params', 'params.mat not found in %s or one folder up.', datapath)
    end
    vox = [D.PData.PDelta(1), D.PData.PDelta(3)] * Pf.P.wl;  % [dx, dz] [m]
end

function [s, A0, nAgree, atBound] = fitCurve(lag, C, R, de, dl, fitAmp)
    m = lag <= R & lag > 0;
    if nnz(m) < 4, error('fit_2D_PSF_autocorr:window', 'Only %d lags in a fit window; use a finer grid or a larger window.', nnz(m)), end
    d = lag(m); c = C(m);
    o = optimoptions('lsqnonlin', 'Display', 'off', 'FunctionTolerance', 1e-14, 'StepTolerance', 1e-12);
    lb = [dl/10, 1]; ub = [max(R, dl), 1];
    if fitAmp, lb(2) = 0.3; end
    s0 = max(de/2, dl/5) * [0.5 1 2]; qs = zeros(numel(s0), 2); costs = zeros(numel(s0), 1);
    for k = 1:numel(s0)
        [qs(k, :), costs(k)] = lsqnonlin(@(q) q(2)*exp(-d.^2/(4*q(1)^2)) - c, [min(max(s0(k), lb(1)*1.01), ub(1)*0.99), 1], lb, ub, o);
    end
    [cost, best] = min(costs);
    nAgree = sum(abs(costs - cost) < 1e-6*cost + 1e-10);
    s = qs(best, 1); A0 = qs(best, 2);
    atBound = abs(s - lb(1)) < 1e-3*lb(1) || abs(s - ub(1)) < 1e-3*ub(1);
end

function [sx, sz, A0, nAgree, atBound] = fitMap(XL, ZL, C2, Rx, Rz, de, vox, fitAmp)
    m = (XL/Rx).^2 + (ZL/Rz).^2 <= 1 & ~(XL == 0 & ZL == 0);
    if nnz(m) < 10, error('fit_2D_PSF_autocorr:window', 'Only %d map lags in the fit window; use a finer grid or a larger window.', nnz(m)), end
    xl = XL(m); zl = ZL(m); c = C2(m);
    o = optimoptions('lsqnonlin', 'Display', 'off', 'FunctionTolerance', 1e-14, 'StepTolerance', 1e-12, 'MaxIterations', 400);
    lb = [vox(1)/10, vox(2)/10, 1]; ub = [max(Rx, vox(1)), max(Rz, vox(2)), 1];
    if fitAmp, lb(3) = 0.3; end
    f = [0.5 1 2]; P0 = [max(de(1)/2, vox(1)/5)*f(:), max(de(2)/2, vox(2)/5)*f(:), ones(3, 1)];
    qs = zeros(3, 3); costs = zeros(3, 1);
    for k = 1:3
        P0(k, 1:2) = min(max(P0(k, 1:2), lb(1:2)*1.01), ub(1:2)*0.99);
        [qs(k, :), costs(k)] = lsqnonlin(@(q) q(3)*exp(-xl.^2/(4*q(1)^2)).*exp(-zl.^2/(4*q(2)^2)) - c, P0(k, :), lb, ub, o);
    end
    [cost, best] = min(costs);
    nAgree = sum(abs(costs - cost) < 1e-6*cost + 1e-10);
    sx = qs(best, 1); sz = qs(best, 2); A0 = qs(best, 3);
    atBound = abs(sx - lb(1)) < 1e-3*lb(1) || abs(sx - ub(1)) < 1e-3*ub(1) || abs(sz - lb(2)) < 1e-3*lb(2) || abs(sz - ub(2)) < 1e-3*ub(2);
end

function plotFit(IQ, info, sigma, cmp, Rx, Rz, A0, win, resMap, pk)
    sx = sigma(1); sz = sigma(2); M = info.measured; um = 1e6; vox = info.voxelSize; [nz, nx] = size(IQ);
    figure('Name', '2D PSF autocorrelation fit'); tiledlayout(2, 2, 'TileSpacing', 'compact');

    nexttile; hx = min(nx, ceil(4*Rx/vox(1))); hz = min(nz, ceil(4*Rz/vox(2)));
    xi = max(1, pk(2)-hx):min(nx, pk(2)+hx); zi = max(1, pk(1)-hz):min(nz, pk(1)+hz);
    E = abs(IQ(zi, xi)); imagesc((xi - pk(2))*vox(1)*um, (zi - pk(1))*vox(2)*um, 20*log10(E/max(E(:)))); axis image; clim([-40 0]); colormap(gca, 'gray'); colorbar
    xlabel('x [\mum]'); ylabel('z [\mum]'); title('|IQ| of the point scatterer [dB]')

    nexttile; x = M.xlag(M.xlag <= 1.6*Rx); dd = linspace(0, max(x), 300);
    plot(x*um, M.lateral(1:numel(x)), 'ko', 'MarkerFaceColor', 'k'); hold on
    plot(dd*um, A0(1)*exp(-dd.^2/(4*sx^2)), 'LineWidth', 1.8); labs = {'measured', sprintf('fit, \\sigma_x = %.0f \\mum', sx*um)};
    if ~isempty(cmp), plot(dd*um, A0(1)*exp(-dd.^2/(4*cmp.sigma(1)^2)), '--', 'LineWidth', 1.4); labs{end+1} = sprintf('comparison, \\sigma_x = %.0f \\mum', cmp.sigma(1)*um); end
    xline(Rx*um, ':', 'HandleVisibility', 'off'); hold off; grid on; xlabel('lateral lag [\mum]'); ylabel('|g_{1,s}|'); ylim([0 1.02]); title('lateral'); legend(labs, 'FontSize', 8)

    nexttile; z = M.zlag(M.zlag <= 1.6*Rz); dd = linspace(0, max(z), 300);
    plot(z*um, M.axial(1:numel(z)), 'ko', 'MarkerFaceColor', 'k'); hold on
    plot(dd*um, A0(2)*exp(-dd.^2/(4*sz^2)), 'LineWidth', 1.8); labs = {'measured', sprintf('fit, \\sigma_z = %.0f \\mum', sz*um)};
    if ~isempty(cmp), plot(dd*um, A0(2)*exp(-dd.^2/(4*cmp.sigma(2)^2)), '--', 'LineWidth', 1.4); labs{end+1} = sprintf('comparison, \\sigma_z = %.0f \\mum', cmp.sigma(2)*um); end
    xline(Rz*um, ':', 'HandleVisibility', 'off'); hold off; grid on; xlabel('axial lag [\mum]'); ylim([0 1.02]); title('axial'); legend(labs, 'FontSize', 8)

    nexttile; R = resMap; R(~win) = NaN;
    imagesc(M.xlag*um, M.zlag*um, R); axis image; colorbar; clim(max(abs(R(:)), [], 'omitnan')*[-1 1]);
    xlim([0 Rx]*um); ylim([0 Rz]*um); xlabel('dx [\mum]'); ylabel('dz [\mum]'); title('2D residual (product model - measured)')
end
