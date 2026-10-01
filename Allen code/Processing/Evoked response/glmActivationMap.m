%% Description: voxel-wise GLM activation analysis for superframe-sampled data (e.g., 3D fUS), using a high-rate stimulus signal
%   The stimulus (sampled at the DAQ rate, e.g. 1000 Hz) is convolved with an HRF at its native rate, and THEN averaged over each
%   superframe's acquisition window. This avoids the aliasing / phase problems of downsampling the stimulus before convolving.
%   A design matrix [stim regressor(s), (temporal derivatives), extra regressors, constant, polynomial drift] is fit to every voxel.
%   Works on any number of data types (e.g., PDI, CDI, vUS) at once, and on any spatial dimensionality.
%
% Inputs:
%   - data: the data to fit. One of:
%       - numeric array, with a time dimension of length numel(sfStarts) (by default the last dimension that matches). e.g., (# x, # y, # z, # superframes)
%       - cell array of the above (e.g., PDIallSF{1:3}); empty cells are passed through as []
%       - scalar struct whose fields are either of the above (e.g., data.PDI, data.CDI, data.vUS)
%   - sfStarts: superframe start times [s] (numel = # superframes), in the same time base as the stim
%   - stim: stimulus waveform sampled at 'Fs', e.g. the 0/1 vector of when the stim is on (cleanStim in evokedResponse3D.m).
%           (# samples x 1), or (# samples x # conditions) for one regressor per condition.
%           Sample n is at time StimStart + (n-1)/Fs.
%
% Name-value inputs:
%   - 'Fs': stim sampling rate [Hz] (default 1000)
%   - 'StimStart': time [s] of the first stim sample in the superframe time base (default 0, i.e. the stim was cropped to the start of acquisition)
%   - 'SFWidth': superframe acquisition duration [s]. The regressor is averaged over [sfStarts, sfStarts + SFWidth] (default 0, i.e. point-sampled)
%   - 'HRF': the hemodynamic response function. One of:
%       - 'gamma': gamma pdf with peak at HRFPeak and shape HRFShape (default)
%       - 'doublegamma': 'gamma' minus a delayed, smaller gamma (undershoot; peak at 2.7 * HRFPeak, 1/6 the amplitude)
%       - 'none': no convolution (regressor = the windowed stim)
%       - a function handle @(t) of time [s] since the stim, or a numeric vector sampled at Fs starting at t = 0
%     The HRF is normalized to sum to 1, so a beta is the response amplitude to a sustained stim in the units of the data.
%   - 'HRFPeak' [s] (default 3), 'HRFShape' (default 4), 'HRFDuration' [s] (default 30): parameters for the built-in HRFs
%   - 'TemporalDerivative': also include the temporal derivative of each regressor, to absorb latency mismatches (default false)
%   - 'DriftOrder': order of the polynomial (Legendre) drift terms (default 2; 0 = constant only)
%   - 'ExtraRegressors': (# superframes x # regressors) nuisance regressors, e.g. motion or pupil (default none)
%   - 'AR': 0 for ordinary least squares, 1 to prewhiten with an AR(1) model of the residuals (default 1).
%           rho is estimated once from the residuals of all voxels, so the same filter is applied to every voxel.
%   - 'Contrast': contrast vector over the design columns (see info.columnNames). It's padded with zeros, so a shorter vector only
%                 needs to cover the stim regressors. Default selects the first condition's regressor.
%   - 'Normalize': per-voxel normalization before fitting. 'none' (default) or 'percent' (100 * (x - mean) / |mean|, so betas are in % change).
%                  A string for all data types, or a struct with a field per data type name (plus an optional 'default' field).
%                  Use 'none' for signed data like CDI or velocity.
%   - 'TimeDim': the time dimension of the data (default 0 = auto: the last dimension whose length matches numel(sfStarts))
%   - 'ChunkSize': # voxels to fit at once (default 2e5), to limit memory use
%
% Outputs:
%   - results: same container type as data (array/cell/struct), but each array is replaced by a struct with fields (each the size of the
%              spatial dimensions of the data):
%       - beta: contrast estimate
%       - t: t-statistic of the contrast
%       - z: the t-statistic converted to a z-score (via the two-sided p-value)
%       - p: two-sided p-value
%       - q: Benjamini-Hochberg FDR-adjusted p-value (across the voxels of that array)
%     Voxels with non-finite data or no variance are NaN.
%   - info: struct with the design matrix (X, rows of unusable superframes are NaN), columnNames, contrast, hrf, tHRF, rho, dof,
%           validSF (superframes used in the fit), and params
%
% Example (see evokedResponse3D.m):
%   dataTypes.PDI = PDIallSF; dataTypes.CDI = CDIallSF; dataTypes.vUS = VallSF_SG;
%   sfWidth = P.numFramesPerBuffer / P.frameRate;
%   [res, info] = glmActivationMap(dataTypes, TD.sfStarts, cleanStim, 'Fs', P.daqrate, 'SFWidth', sfWidth, ...
%       'Normalize', struct('PDI', 'percent', 'default', 'none'));
%   % res.PDI{3}.t is (# x, # y, # z). Threshold, e.g. res.PDI{3}.q < 0.05 & res.PDI{3}.beta > 0

function [results, info] = glmActivationMap(data, sfStarts, stim, opts)
    arguments
        data
        sfStarts (:, 1) double
        stim double
        opts.Fs (1, 1) double = 1000
        opts.StimStart (1, 1) double = 0
        opts.SFWidth (1, 1) double = 0
        opts.HRF = 'gamma'
        opts.HRFPeak (1, 1) double = 3
        opts.HRFShape (1, 1) double = 4
        opts.HRFDuration (1, 1) double = 30
        opts.TemporalDerivative (1, 1) logical = false
        opts.DriftOrder (1, 1) double {mustBeInteger, mustBeNonnegative} = 2
        opts.ExtraRegressors double = []
        opts.AR (1, 1) double {mustBeMember(opts.AR, [0, 1])} = 1
        opts.Contrast (1, :) double = []
        opts.Normalize = 'none'
        opts.TimeDim (1, 1) double = 0
        opts.ChunkSize (1, 1) double = 2e5
    end

    nSF = numel(sfStarts);
    fs = opts.Fs;
    if isvector(stim), stim = stim(:); end
    [nStim, nCond] = size(stim);
    if any(~isfinite(stim), 'all'), error('stim contains non-finite values.'), end
    if ~isempty(opts.ExtraRegressors) && size(opts.ExtraRegressors, 1) ~= nSF
        error('ExtraRegressors must have one row per superframe (%d).', nSF)
    end

    %% Stim regressors: convolve with the HRF at the stim rate, then average over each superframe window
    [h, tHRF] = buildHRF(opts.HRF, fs, opts.HRFDuration, opts.HRFPeak, opts.HRFShape);

    i0 = round((sfStarts - opts.StimStart) * fs) + 1; % First stim sample in each superframe window
    i1 = max(round((sfStarts + opts.SFWidth - opts.StimStart) * fs), i0); % Last stim sample in each superframe window
    validSF = i0 >= 1 & i1 <= nStim;

    reg = nan(nSF, nCond);
    dreg = nan(nSF, nCond);
    for c = 1:nCond
        xhr = fftConv(stim(:, c), h);
        reg(validSF, c) = windowAvg(xhr, i0(validSF), i1(validSF));
        if opts.TemporalDerivative
            dreg(validSF, c) = windowAvg(gradient(xhr) * fs, i0(validSF), i1(validSF));
        end
    end

    %% Design matrix: [stim regressors, (derivatives), extra, constant, polynomial drift]
    colNames = compose("stim%d", 1:nCond);
    X = reg;
    if opts.TemporalDerivative
        X = [X, dreg];
        colNames = [colNames, compose("dstim%d", 1:nCond)];
    end
    if ~isempty(opts.ExtraRegressors)
        X = [X, opts.ExtraRegressors];
        colNames = [colNames, compose("extra%d", 1:size(opts.ExtraRegressors, 2))];
        validSF = validSF & all(isfinite(opts.ExtraRegressors), 2);
    end

    tv = sfStarts(validSF);
    if numel(tv) < 2, error('Fewer than 2 superframes overlap the stim recording.'), end
    xs = 2 * (sfStarts - tv(1)) / (tv(end) - tv(1)) - 1; % Time scaled to [-1, 1] for the Legendre polynomials
    drift = zeros(nSF, opts.DriftOrder + 1);
    drift(:, 1) = 1;
    if opts.DriftOrder >= 1, drift(:, 2) = xs; end
    for n = 1:opts.DriftOrder - 1 % Bonnet's recursion: (n+1) P_{n+1} = (2n+1) x P_n - n P_{n-1}
        drift(:, n + 2) = ((2*n + 1) * xs .* drift(:, n + 1) - n * drift(:, n)) / (n + 1);
    end
    X = [X, drift];
    colNames = [colNames, "const", compose("drift%d", 1:opts.DriftOrder)];
    X(~validSF, :) = NaN;
    p = size(X, 2);
    Xv = X(validSF, :);

    contrast = zeros(1, p);
    if isempty(opts.Contrast)
        contrast(1) = 1;
    else
        if numel(opts.Contrast) > p, error('Contrast has more elements than the %d design columns.', p), end
        contrast(1:numel(opts.Contrast)) = opts.Contrast;
    end

    ctx = struct('nSF', nSF, 'validSF', validSF, 'Xv', Xv, 'contrast', contrast, 'AR', opts.AR, ...
        'Normalize', {opts.Normalize}, 'TimeDim', opts.TimeDim, 'ChunkSize', opts.ChunkSize);

    %% Fit every data type
    [results, rhos] = walk(data, '', ctx);

    %% Info
    info.X = X;
    info.columnNames = colNames;
    info.contrast = contrast;
    info.hrf = h;
    info.tHRF = tHRF;
    info.rho = rhos;
    info.dof = nnz(validSF) - rank(Xv); % Without prewhitening losing any
    info.validSF = validSF;
    info.params = opts;
end

%% Recursively apply to arrays inside cells / structs
function [res, rho] = walk(x, name, ctx)
    if isnumeric(x) || islogical(x)
        if isempty(x)
            res = []; rho = [];
        else
            [res, rho] = fitArray(x, name, ctx);
        end
    elseif iscell(x)
        res = cell(size(x)); rho = cell(size(x));
        for k = 1:numel(x)
            [res{k}, rho{k}] = walk(x{k}, name, ctx);
        end
    elseif isstruct(x) && isscalar(x)
        res = struct(); rho = struct();
        for f = fieldnames(x).'
            fn = f{1};
            [rk, rhok] = walk(x.(fn), fn, ctx);
            res.(fn) = rk; rho.(fn) = rhok;
        end
    else
        error('Unsupported data type: %s. Use numeric arrays, cells, or scalar structs.', class(x))
    end
end

%% Fit the GLM to one numeric array (time along dimension td)
function [res, rho] = fitArray(x, name, ctx)
    mode = resolveNormalize(ctx.Normalize, name);

    td = ctx.TimeDim;
    if td == 0, td = find(size(x) == ctx.nSF, 1, 'last'); end
    if isempty(td) || size(x, td) ~= ctx.nSF
        error('Data ''%s'' (size %s) has no time dimension of length %d.', name, mat2str(size(x)), ctx.nSF)
    end
    order = [td, setdiff(1:max(ndims(x), td), td)];
    Xm = permute(x, order); % Time first
    spatialSz = size(Xm, 2:ndims(Xm));
    nVox = prod(spatialSz);
    Y = reshape(Xm, ctx.nSF, nVox);
    clear Xm

    Xv = ctx.Xv;
    [nT, p] = size(Xv);
    chunks = [1:ctx.ChunkSize:nVox; min((1:ctx.ChunkSize:nVox) + ctx.ChunkSize - 1, nVox)];

    %% Pass 1 (if prewhitening): OLS fit to estimate the residual lag-1 autocorrelation, pooled over all voxels
    rho = 0;
    if ctx.AR == 1
        num = 0; den = 0;
        for ch = chunks
            [Yc, ~] = loadChunk(Y, ch, ctx.validSF, mode);
            R = Yc - Xv * (Xv \ Yc);
            num = num + sum(R(2:end, :) .* R(1:end-1, :), 'all');
            den = den + sum(R(1:end-1, :).^2, 'all');
        end
        rho = min(max(num / den, 0), 0.99);
    end

    %% Pass 2: (prewhitened) fit
    Xw = whiten(Xv, rho);
    XtXi = pinv(Xw.' * Xw);
    dof = nT - rank(Xw);
    c = ctx.contrast;
    cVar = c * XtXi * c.'; % Variance of the contrast, in units of the residual variance
    cProj = c * XtXi * Xw.'; % Contrast estimate = cProj * Yw

    beta = nan(1, nVox); tstat = nan(1, nVox);
    for ch = chunks
        cols = ch(1):ch(2);
        [Yc, bad] = loadChunk(Y, ch, ctx.validSF, mode);
        Yw = whiten(Yc, rho);
        est = cProj * Yw;
        R = Yw - Xw * (XtXi * (Xw.' * Yw));
        s2 = sum(R.^2, 1) / dof;
        tc = est ./ sqrt(s2 * cVar); % NaN where there's no variance (0/0)
        est(bad) = NaN; tc(bad) = NaN;
        beta(cols) = est;
        tstat(cols) = tc;
    end

    % p-values from the t distribution (via the incomplete beta function), z-scores, and FDR q-values
    pval = betainc(dof ./ (dof + tstat.^2), dof / 2, 0.5);
    z = sign(tstat) .* sqrt(2) .* erfcinv(pval);
    q = nan(size(pval));
    ok = find(~isnan(pval));
    [ps, ord] = sort(pval(ok));
    qs = ps .* numel(ok) ./ (1:numel(ok));
    qs = min(flip(cummin(flip(qs))), 1);
    q(ok(ord)) = qs;

    shape = @(v) reshape(v, [spatialSz, 1]);
    res.beta = shape(beta);
    res.t = shape(tstat);
    res.z = shape(z);
    res.p = shape(pval);
    res.q = shape(q);
end

%% Get a chunk of voxels (valid superframes only), normalized, with non-finite voxels zeroed and flagged
function [Yc, bad] = loadChunk(Y, ch, validSF, mode)
    Yc = double(Y(validSF, ch(1):ch(2)));
    if strcmp(mode, 'percent')
        m = mean(Yc, 1);
        Yc = 100 * (Yc - m) ./ abs(m);
    end
    bad = any(~isfinite(Yc), 1);
    Yc(:, bad) = 0;
end

%% AR(1) prewhitening filter along the first dimension
function A = whiten(A, rho)
    if rho == 0, return, end
    A = [sqrt(1 - rho^2) * A(1, :); A(2:end, :) - rho * A(1:end-1, :)];
end

%% Mean of x over each window [i0(k), i1(k)] (sample indices)
function m = windowAvg(x, i0, i1)
    cs = [0; cumsum(x(:))];
    m = (cs(i1 + 1) - cs(i0)) ./ (i1 - i0 + 1);
end

%% Causal convolution of x with h, truncated to the length of x (via FFT, since the stim can have millions of samples)
function y = fftConv(x, h)
    n = numel(x) + numel(h) - 1;
    nfft = 2^nextpow2(n);
    y = real(ifft(fft(x(:), nfft) .* fft(h(:), nfft)));
    y = y(1:numel(x));
end

%% Build the HRF at sampling rate fs, normalized to sum to 1
function [h, t] = buildHRF(spec, fs, dur, peak, shape)
    t = (0:round(dur * fs)).' / fs;
    if isa(spec, 'function_handle')
        h = spec(t);
    elseif isnumeric(spec)
        h = spec(:);
        t = (0:numel(h) - 1).' / fs;
    else
        gam = @(tt, pk) (tt.^(shape - 1) .* exp(-tt / (pk / (shape - 1)))); % Gamma with peak at pk
        switch validatestring(char(spec), {'gamma', 'doublegamma', 'none'})
            case 'gamma'
                h = gam(t, peak);
                h = h / sum(h);
            case 'doublegamma'
                h1 = gam(t, peak); h1 = h1 / sum(h1);
                h2 = gam(t, 2.7 * peak); h2 = h2 / sum(h2);
                h = h1 - h2 / 6;
            case 'none'
                h = 1; t = 0;
        end
    end
    h = h(:);
    if sum(h) == 0, error('HRF sums to 0; cannot normalize.'), end
    h = h / sum(h);
end

%% Choose the normalization for a data type name
function mode = resolveNormalize(norm, name)
    if isstruct(norm)
        if ~isempty(name) && isfield(norm, name)
            mode = norm.(name);
        elseif isfield(norm, 'default')
            mode = norm.default;
        else
            mode = 'none';
        end
    else
        mode = norm;
    end
    mode = validatestring(char(mode), {'none', 'percent'});
end
