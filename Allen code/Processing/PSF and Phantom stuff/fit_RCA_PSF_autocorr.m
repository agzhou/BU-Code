function [sigma, info] = fit_RCA_PSF_autocorr(varargin)
%% Description:
%   Fit the three sigmas of the RCA g1,s model (sigma_narrow, sigma_wide, sigma_axial) to the MEASURED
%   AUTOCORRELATION of a simulated (or measured) point-scatterer volume, instead of to its envelope abs(IQ)
%   as fit_RCA_PSF_v2.m does. The output plugs directly into sigma_RCA of main_vUS_3D.m /
%   vUS_3D_quad_RCA_xy_vec.m.
%
%   Why the autocorrelation: a scatterer that has moved by d relative to its position at the previous time
%   point changes the voxel signal from h(r) to h(r - d), where h is the complex point-scatterer image. For
%   scatterers distributed uniformly in space, the single-scatterer field autocorrelation is therefore
%
%       g1,s(d) = int conj(h(r)) h(r + d) dr / int |h(r)|^2 dr,
%
%   the normalized autocorrelation of the COMPLEX PSF. Only |g1,s(d)| is fit here: the carrier phase,
%   exp(i 2 k0 v_z tau f), is carried explicitly by the vUS model. The envelope abs(IQ) is not the same
%   thing. If the complex field has any phase structure beyond the carrier (as the compounded RCA PSF does),
%   it decorrelates faster than its envelope's autocorrelation, so a sigma fit to abs(IQ) makes the model
%   decay too slowly. On the 5-angle RC15gV simulation this is a factor of about 1.4 to 1.7 in sigma.
%
%   Model (identical to the single-scatterer term of vUS_3D_quad_RCA_xy_vec.m, plug flow), for a lateral
%   displacement (dx, dy) and axial displacement dz:
%
%       |g1,s| = Clat(dx, dy) * Cax(dz)
%       Cax    = exp(-dz^2 / (4 sigma_axial^2))
%       Clat   = [ (E1 + E4)/2 + c E23 ] / (1 + c),    S = sigma_narrow^2 + sigma_wide^2,   c = 2 sigma_narrow sigma_wide / S
%       E1  = exp(-dx^2/(4 sigma_narrow^2) - dy^2/(4 sigma_wide^2))
%       E4  = exp(-dx^2/(4 sigma_wide^2)   - dy^2/(4 sigma_narrow^2))
%       E23 = exp(-(dx^2 + dy^2) / (2 S))
%
%   (E1 and E4 are the autocorrelations of the two "arms" of the cross, E23 is their cross-correlation. The
%   model has no free amplitude: |g1,s(0)| = 1 by construction.)
%
%   Procedure:
%     1. Autocorrelation by FFT with zero padding (no wrap-around): the lateral map at dz = 0 is the sum over
%        z of each slice's 2D autocorrelation; the axial curve at dx = dy = 0 is the sum over (x, y) of each
%        column's 1D autocorrelation. Both are exact at integer-voxel lags. No Fourier upsampling and no
%        spline interpolation is used: both ring when the PSF is only a voxel or two wide. If the PSF is
%        undersampled, reconstruct on a finer grid (a quarter wavelength is enough for RC15gV).
%     2. Normalize by the autocorrelation at zero lag (the total energy), so |g1,s(0)| = 1. This includes the
%        energy of the low-amplitude cross-artifact arms, which is what uniformly distributed scatterers see.
%        info.energyFraction reports how much of it lies inside the final fit window.
%     3. Separability: the model is a product of a lateral and an axial factor. This is checked directly at a
%        few combined lags (info.separability), not assumed.
%     4. Lateral fit, [sigma_narrow, sigma_wide], to the first-quadrant lag map. |A(dx, dy)| and |A(dx, -dy)|
%        are averaged. Fit parameters are [sigma_narrow, sigma_wide/sigma_narrow] with the ratio >= 1, as in
%        fit_RCA_PSF_v2.m, since the model is symmetric under narrow <-> wide. Nine start points; the best is
%        kept and the number that agree is reported.
%     5. Axial fit, sigma_axial, to |A(0, 0, dz)|.
%     6. Window: like fit_RCA_PSF_v2.m, in two stages. A pilot fit uses a wide window from the measured 1/e
%        decay lags (6x laterally, 4x axially). The final fit uses a window of WindowSigmas x the pilot sigma
%        (sigma_wide laterally, sigma_axial axially; default 3.5, where the model has fallen to exp(-3.5^2/4) =
%        4.7%). Far beyond that the real PSF has arms the model can't represent, which would only bias
%        sigma_wide upward. It is deliberately not iterated to self-consistency, for the same reason.
%
% Usage:
%   fit_RCA_PSF_autocorr                          % interactive: pick the folder with PSF.mat and PData.mat
%   sigma = fit_RCA_PSF_autocorr(datapath)         % folder with PSF.mat, PData.mat; params.mat there or one up
%   sigma = fit_RCA_PSF_autocorr(IQ, voxelSize)    % IQ [x, y, z] complex, voxelSize [dx, dy, dz] in meters
%   [sigma, info] = fit_RCA_PSF_autocorr(..., 'CompareSigma', [57.6 286.5 50.1]*1e-6, 'Plot', true)
%
% Inputs:
%   IQ: [nx, ny, nz] complex image of a single point scatterer (e.g., PSF.mat from the PSF simulations).
%       Dimensions are [x, y, z] as in fit_RCA_PSF_v2.m; z is axial.
%   voxelSize: [dx, dy, dz] voxel size [m], or a scalar. PData.PDelta(1:3) .* P.wl for Verasonics data.
%   datapath: folder containing PSF.mat (IQ), PData.mat (PData) and params.mat (P) next to it or one folder up.
%
% Name-value options:
%   'Plot': draw data vs model (default true only when called with no arguments)
%   'CompareSigma': [narrow, wide, axial] [m], e.g. the sigma_RCA you are using now. Reported and plotted
%       against the same measured autocorrelation, in the same window.
%   'WindowSigmas': final fit window half-width in units of sigma (default 3.5)
%   'MaxLateralLag', 'MaxAxialLag': override the final fit window [m] (default: derived as above)
%   'CropLateral', 'CropAxial': keep only +-this distance [m] around the brightest voxel before the
%       autocorrelation (default: use the whole volume). See the caveat on the arms below.
%   'FitAmplitude': also fit an overall amplitude A0 <= 1 (default false). Use it for a MEASURED point
%       scatterer (noisy): white noise adds energy only at zero lag, so every lag > 0 is scaled by
%       signal/(signal + noise) after normalization, and a fit without A0 shrinks sigma_narrow to compensate
%       (about 11% for 0.3% noise in the synthetic test). Leave it off for noise-free simulations, where A0
%       would only trade off against the model's tail. info.amplitude reports the fitted [lateral, axial] A0.
%   'CheckSeparability': direct check at a few combined lags (default true; ~1 s)
%   'Verbose': print a summary (default true)
%
% Outputs:
%   sigma: [sigma_narrow, sigma_wide, sigma_axial] [m], FIELD (amplitude) sigmas in the convention of
%       sigma_RCA (g1,s has exp(-d^2/(4 sigma^2)))
%   info: struct with the fit window, residual statistics, the measured and fitted curves, the energy fraction
%       in the window, the separability check, a sampling check, and CompareSigma statistics.
%
% Caveats:
%   - The result depends on how much of the cross-artifact arms the volume contains, because g1,s is
%     normalized by the TOTAL energy and the arms carry most of it (about 70% outside a +-0.6 mm x +-0.17 mm
%     window in the 5-angle RC15gV simulation). Cropping the 5-angle simulation laterally from +-4.4 mm to
%     +-1 mm moved [sigma_narrow, sigma_wide, sigma_axial] from [32, 171, 36] to [39, 202, 41] um. Treat the
%     spread as the uncertainty (use 'CropLateral' to see it), and prefer the extent that matches the
%     field of view you image. info.energyFraction reports how much energy is in the fit window.
%   - One depth, one transmit sequence, an ideal medium. The RCA PSF changes with depth; repeat for the depths
%     you image (the simulated point scatterer's z position sets it).
%   - Sampling: the narrow arm is about 30 to 60 um. The function warns if a voxel is larger than the fitted
%     sigma_narrow or sigma_axial, because the autocorrelation at integer lags is then not trustworthy.
%   - The Gaussian-cross form is an approximation. The residual is largest on the lateral diagonal, where the
%     measured autocorrelation drops to ~0 by about 3 sigma_narrow and the model keeps a small tail. Check
%     info.rmsLateral and the plot before trusting the fit.

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
        if nargin < 2, error('fit_RCA_PSF_autocorr:voxelSize', 'Provide the voxel size as the second input.'), end
        IQ = varargin{1}; vox = varargin{2}; rest = varargin(3:end);
    end
    vox = double(vox(:)).';
    if isscalar(vox), vox = vox*[1 1 1]; end
    if ndims(IQ) ~= 3 || numel(vox) ~= 3 || any(vox <= 0)
        error('fit_RCA_PSF_autocorr:input', 'IQ must be 3D [x, y, z] and voxelSize [dx, dy, dz] positive.')
    end
    p = inputParser;
    p.addParameter('Plot', interactive, @(x) isscalar(x));
    p.addParameter('CompareSigma', [], @(x) isempty(x) || (isnumeric(x) && numel(x) == 3));
    p.addParameter('WindowSigmas', 3.5, @(x) isscalar(x) && x > 0);
    p.addParameter('MaxLateralLag', [], @(x) isempty(x) || (isscalar(x) && x > 0));
    p.addParameter('MaxAxialLag', [], @(x) isempty(x) || (isscalar(x) && x > 0));
    p.addParameter('CropLateral', [], @(x) isempty(x) || (isscalar(x) && x > 0));
    p.addParameter('CropAxial', [], @(x) isempty(x) || (isscalar(x) && x > 0));
    p.addParameter('FitAmplitude', false, @(x) isscalar(x));
    p.addParameter('CheckSeparability', true, @(x) isscalar(x));
    p.addParameter('Verbose', true, @(x) isscalar(x));
    p.parse(rest{:});
    opt = p.Results;
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

    %% 1. Autocorrelation at integer lags (exact; no interpolation)
    [Rlat, Rax, E0] = psfAutocorr(IQ);
    Cax = abs(Rax(1:nz)).' / E0;                        % [1, nz] axial lags 0 ... nz-1 (a row, like zlag)
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
    [sn1, sw1] = fitLateral(DX, DY, C, Rlat1, de_lat, vox, opt.FitAmplitude);
    sa1 = fitAxial(zlag, Cax, Rax1, de_ax, vox(3), opt.FitAmplitude);

    %% 3. Final window and fit
    if isempty(opt.MaxLateralLag), Rlat2 = opt.WindowSigmas * sw1; else, Rlat2 = opt.MaxLateralLag; end
    if isempty(opt.MaxAxialLag),   Rax2  = opt.WindowSigmas * sa1; else, Rax2  = opt.MaxAxialLag;   end
    Rlat2 = min(max(Rlat2, 6*max(vox(1:2))), Rlat_all);   % at least 6 voxels, at most what the volume holds
    Rax2  = min(max(Rax2,  6*vox(3)), zlag(end));
    [sn, sw, A0lat, nAgree, atBound] = fitLateral(DX, DY, C, Rlat2, de_lat, vox, opt.FitAmplitude);
    [sa, A0ax] = fitAxial(zlag, Cax, Rax2, de_ax, vox(3), opt.FitAmplitude);
    sigma = [sn, sw, sa];

    %% 4. Diagnostics
    mLat = rl <= Rlat2 & rl > 0;
    resLat = A0lat*latModel(sn, sw, DX, DY) - C;
    mAx = zlag <= Rax2 & zlag > 0;
    resAx = A0ax*exp(-zlag.^2/(4*sa^2)) - Cax;
    info.sigma = sigma;
    info.sigma_um = sigma * 1e6;
    info.voxelSize = vox;
    info.window = struct('lateral', Rlat2, 'axial', Rax2, 'pilotLateral', Rlat1, 'pilotAxial', Rax1, ...
        'pilotSigma', [sn1, sw1, sa1], 'nLateralLags', nnz(mLat), 'nAxialLags', nnz(mAx));
    info.rmsLateral = sqrt(mean(resLat(mLat).^2));
    info.maxLateral = max(abs(resLat(mLat)));
    info.rmsAxial = sqrt(mean(resAx(mAx).^2));
    info.maxAxial = max(abs(resAx(mAx)));
    info.amplitude = [A0lat, A0ax];                      % 1 unless 'FitAmplitude'
    info.startsAgreeing = nAgree;
    info.atBound = atBound;
    info.measured = struct('DX', DX, 'DY', DY, 'lateral', C, 'zlag', zlag, 'axial', Cax);

    % Energy of the PSF inside the fit window, around the brightest voxel
    [~, imax] = max(abs(IQ(:)));
    [px, py, pz] = ind2sub([nx, ny, nz], imax);
    wx = ceil(Rlat2/vox(1)); wy = ceil(Rlat2/vox(2)); wz = ceil(Rax2/vox(3));
    Ew = sum(abs(IQ(max(1, px-wx):min(nx, px+wx), max(1, py-wy):min(ny, py+wy), max(1, pz-wz):min(nz, pz+wz))).^2, 'all');
    info.energyFraction = Ew / E0;

    % Sampling
    info.undersampled = [vox(1) > sn, vox(2) > sn, vox(3) > sa];
    if any(info.undersampled)
        warning('fit_RCA_PSF_autocorr:undersampled', ...
            ['Voxel size [%.1f %.1f %.1f] um is larger than the fitted sigma_narrow / sigma_axial (%.1f / %.1f um). ' ...
             'The integer-lag autocorrelation is not trustworthy; reconstruct the PSF on a finer grid.'], ...
             vox(1)*1e6, vox(2)*1e6, vox(3)*1e6, sn*1e6, sa*1e6)
    end
    if any(atBound)
        warning('fit_RCA_PSF_autocorr:atBound', 'A lateral sigma hit its bound; the lateral fit is not well determined (window or sampling).')
    end

    % Separability, directly at a few combined lags
    info.separability = [];
    if opt.CheckSeparability
        ln = max(1, round(sn/vox(1)) + 1); lw = max(ln + 1, round(sw/vox(1)));
        lz = max(1, round(sa/vox(3)) + [1 2]);
        lagsLat = [ln 0; lw 0; ln ln]; rows = zeros(0, 6);
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

    % Compare against another sigma set (e.g., the one currently in use) over the same window
    info.compare = [];
    if ~isempty(opt.CompareSigma)
        cs = double(opt.CompareSigma(:)).';
        cl = A0lat*latModel(cs(1), cs(2), DX, DY) - C; ca = A0ax*exp(-zlag.^2/(4*cs(3)^2)) - Cax;
        info.compare = struct('sigma', cs, 'rmsLateral', sqrt(mean(cl(mLat).^2)), 'maxLateral', max(abs(cl(mLat))), ...
            'rmsAxial', sqrt(mean(ca(mAx).^2)), 'maxAxial', max(abs(ca(mAx))));
    end

    %% 5. Report
    if opt.Verbose
        fprintf('Autocorrelation fit of the RCA PSF (voxel %.1f x %.1f x %.1f um, volume %d x %d x %d)\n', ...
            vox(1)*1e6, vox(2)*1e6, vox(3)*1e6, nx, ny, nz);
        fprintf('  pilot sigma_narrow / sigma_wide / sigma_axial = %.1f / %.1f / %.1f um\n', 1e6*sn1, 1e6*sw1, 1e6*sa1);
        fprintf('  window: +-%.0f um lateral (%d lags) x +-%.0f um axial (%d lags); %.0f%% of the PSF energy is inside it\n', ...
            1e6*Rlat2, nnz(mLat), 1e6*Rax2, nnz(mAx), 100*info.energyFraction);
        fprintf('Final: sigma_narrow / sigma_wide / sigma_axial = %.1f / %.1f / %.1f um (field)\n', info.sigma_um);
        fprintf('  lateral RMS %.4f (max %.3f), axial RMS %.4f (max %.3f), %d/9 start points agree\n', ...
            info.rmsLateral, info.maxLateral, info.rmsAxial, info.maxAxial, nAgree);
        if opt.FitAmplitude
            fprintf('  fitted amplitude (1 = no noise floor): lateral %.3f, axial %.3f\n', A0lat, A0ax);
        end
        if ~isempty(info.separability)
            fprintf('  separability: max |direct 3D - lateral*axial| = %.4f over %d lags\n', info.maxSeparabilityError, height(info.separability));
        end
        if ~isempty(info.compare)
            fprintf('  for comparison, sigma = [%.1f %.1f %.1f] um: lateral RMS %.4f (max %.3f), axial RMS %.4f (max %.3f)\n', ...
                1e6*info.compare.sigma, info.compare.rmsLateral, info.compare.maxLateral, info.compare.rmsAxial, info.compare.maxAxial);
        end
    end

    %% 6. Plot
    if opt.Plot
        plotFit(info, sigma, opt.CompareSigma, Rlat2, Rax2, [A0lat, A0ax])
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
        error('fit_RCA_PSF_autocorr:params', 'params.mat not found in %s or one folder up.', datapath)
    end
    vox = D.PData.PDelta(1:3) * Pf.P.wl;                     % [m]
end

function [Rlat, Rax, E0] = psfAutocorr(IQ)
    % Lateral autocorrelation at dz = 0 (sum over z of the 2D slice autocorrelations) and axial autocorrelation
    % at dx = dy = 0 (sum over x, y of the 1D column autocorrelations), by FFT with zero padding. Chunked, so
    % the memory use stays around the size of IQ.
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

function Y = latModel(sn, sw, DX, DY)
    S = sn^2 + sw^2; c = 2*sn*sw/S;
    Y = ((exp(-DX.^2/(4*sn^2) - DY.^2/(4*sw^2)) + exp(-DX.^2/(4*sw^2) - DY.^2/(4*sn^2)))/2 ...
        + c*exp(-(DX.^2 + DY.^2)/(2*S))) / (1 + c);
end

function [sn, sw, A0, nAgree, atBound] = fitLateral(DX, DY, C, R, de, vox, fitAmp)
    r = hypot(DX, DY);
    m = r <= R & r > 0;
    if nnz(m) < 10, error('fit_RCA_PSF_autocorr:window', 'Only %d lateral lags in the fit window; use a finer grid or a larger window.', nnz(m)), end
    dx = DX(m); dy = DY(m); c = C(m);
    lb = [min(vox(1:2))/10, 1, 1]; ub = [R, 50, 1];
    if fitAmp, lb(3) = 0.3; end                       % amplitude is held at 1 unless fitted; a noise floor only lowers it
    s0 = de/2 * [0.5 1 2]; r0 = [2 5 10];
    [S0, R0] = ndgrid(s0, r0); P0 = [min(max(S0(:), lb(1)*1.01), ub(1)*0.99), R0(:), ones(numel(R0), 1)];
    o = optimoptions('lsqnonlin', 'Display', 'off', 'FunctionTolerance', 1e-14, 'StepTolerance', 1e-12, 'MaxIterations', 400);
    ps = zeros(size(P0)); costs = zeros(size(P0, 1), 1);
    for k = 1:size(P0, 1)
        [ps(k, :), costs(k)] = lsqnonlin(@(q) q(3)*latModel(q(1), q(1)*q(2), dx, dy) - c, P0(k, :), lb, ub, o);
    end
    [cost, best] = min(costs);
    nAgree = sum(abs(costs - cost) < 1e-6*cost + 1e-10);   % relative, with an absolute floor for near-perfect fits
    sn = ps(best, 1); sw = ps(best, 1)*ps(best, 2); A0 = ps(best, 3);
    atBound = [abs(sn - lb(1)) < 1e-3*lb(1), abs(sw - ub(1)) < 1e-3*ub(1)];
end

function [sa, A0] = fitAxial(zlag, Cax, R, de, dz, fitAmp)
    m = zlag <= R & zlag > 0;
    if nnz(m) < 4, error('fit_RCA_PSF_autocorr:window', 'Only %d axial lags in the fit window; use a finer grid or a larger window.', nnz(m)), end
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

function plotFit(info, sigma, cmpSigma, Rlat, Rax, A0)
    sn = sigma(1); sw = sigma(2); sa = sigma(3);
    M = info.measured; vox = info.voxelSize;
    um = 1e6;
    figure('Name', 'RCA PSF autocorrelation fit'); tiledlayout(2, 2, 'TileSpacing', 'compact');

    nexttile; ax = M.DY == 0 & M.DX <= 1.6*Rlat; d = M.DX(ax); dd = linspace(0, max(d), 300);
    plot(d*um, M.lateral(ax), 'ko', 'MarkerFaceColor', 'k'); hold on
    plot(dd*um, A0(1)*latModel(sn, sw, dd, 0*dd), 'LineWidth', 1.6)
    if ~isempty(cmpSigma), plot(dd*um, A0(1)*latModel(cmpSigma(1), cmpSigma(2), dd, 0*dd), '--', 'LineWidth', 1.4); end
    xline(Rlat*um, ':', 'HandleVisibility', 'off'); hold off; grid on; xlabel('lateral lag along x [\mum]'); ylabel('|g_{1,s}|'); title('x axis'); ylim([0 1.02])

    nexttile; dg = abs(M.DX - M.DY) < 0.5*min(vox(1:2)) & M.DX*sqrt(2) <= 1.6*Rlat; d = M.DX(dg)*sqrt(2); dd = linspace(0, max(d), 300);
    plot(d*um, M.lateral(dg), 'ko', 'MarkerFaceColor', 'k'); hold on
    plot(dd*um, A0(1)*latModel(sn, sw, dd/sqrt(2), dd/sqrt(2)), 'LineWidth', 1.6)
    if ~isempty(cmpSigma), plot(dd*um, A0(1)*latModel(cmpSigma(1), cmpSigma(2), dd/sqrt(2), dd/sqrt(2)), '--', 'LineWidth', 1.4); end
    xline(Rlat*um, ':'); hold off; grid on; xlabel('lateral lag along the diagonal [\mum]'); title('diagonal'); ylim([0 1.02])
    legend({'measured', sprintf('fit [%.0f %.0f %.0f] \\mum', sigma*um), 'comparison'}, 'Location', 'northeast')

    nexttile; z = M.zlag(M.zlag <= 1.6*Rax); c = M.axial(1:numel(z)); dd = linspace(0, max(z), 300);
    plot(z*um, c, 'ko', 'MarkerFaceColor', 'k'); hold on
    plot(dd*um, A0(2)*exp(-dd.^2/(4*sa^2)), 'LineWidth', 1.6)
    if ~isempty(cmpSigma), plot(dd*um, A0(2)*exp(-dd.^2/(4*cmpSigma(3)^2)), '--', 'LineWidth', 1.4); end
    xline(Rax*um, ':', 'HandleVisibility', 'off'); hold off; grid on; xlabel('axial lag [\mum]'); title('axial'); ylim([0 1.02])

    nexttile; m = hypot(M.DX, M.DY) <= Rlat; R = A0(1)*latModel(sn, sw, M.DX, M.DY) - M.lateral; R(~m) = NaN;
    imagesc(M.DX(:, 1)*um, M.DY(1, :)*um, R.'); axis image; colorbar; clim(max(abs(R(:)), [], 'omitnan')*[-1 1]);
    xlim([0 Rlat]*um); ylim([0 Rlat]*um);
    xlabel('dx [\mum]'); ylabel('dy [\mum]'); title('lateral residual (model - measured)')
end
