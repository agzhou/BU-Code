%% Description: peri-stimulus (evoked response) trial averaging of superframe-sampled data
%   Each trial's data is interpolated (in time) onto a common peri-stimulus time grid, optionally baseline-
%   normalized, and then averaged across trials. Superframe times don't need to be phase-locked to the trials.
%   Works on any number of data types (e.g., PDI, CDI, vUS) at once, and on any spatial dimensionality.
%
% Inputs:
%   - data: the data to average. One of:
%       - numeric array, with a time dimension of length numel(sfStarts) (by default the last dimension that matches). e.g., (# x, # y, # z, # superframes)
%       - cell array of the above (e.g., PDIallSF{1:3}); empty cells are passed through as []
%       - scalar struct whose fields are either of the above (e.g., data.PDI, data.CDI, data.vUS)
%   - sfStarts: superframe start times [s] (numel = # superframes), in the same time base as stimOnsets
%   - stimOnsets: stim onset times [s] for each trial
%
% Name-value inputs:
%   - 'SFWidth': superframe acquisition duration [s]. The superframe is assigned to its center time, sfStarts + SFWidth/2. (default 0, i.e. start times are used as-is)
%   - 'Pre', 'Post': time before / after the stim onset to include [s] (default 5, 30)
%   - 'Dt': peri-stimulus grid spacing [s] (default: median superframe period)
%   - 'BaselineWindow': [start, end) in time relative to the stim onset [s] used for the baseline (default [-Pre, 0])
%   - 'Normalize': how to normalize each trial to its baseline. Either a string for all data types or a struct with a field per data type
%                  name (plus an optional 'default' field). Options:
%                    'percent': 100 * (x - b) / |b|    (default; use for strictly-positive data like PDI)
%                    'diff':    x - b                  (use for signed data like CDI or velocity, or when the baseline can be ~0)
%                    'ratio':   x / b
%                    'zscore':  (x - b) / std(baseline samples)
%                    'none'
%   - 'Method': interp1 method for the temporal interpolation (default 'linear')
%   - 'MaxGap': don't interpolate across superframe gaps longer than this [s]; those grid points are set to NaN for that trial.
%               (default 2 * median superframe period, i.e. interpolates across at most one dropped superframe)
%   - 'TimeDim': the time dimension of the data (default 0 = auto: the last dimension whose length matches numel(sfStarts))
%   - 'TrialsToUse': indices of trials (into stimOnsets) to include (default all)
%   - 'KeepTrials': also return the individual trial data in info.trialData (default false; can be large)
%
% Outputs:
%   - evoked: trial-averaged data on the peri-stimulus grid. Same container type as data (array/cell/struct); the time dimension
%             is replaced by the peri-stimulus time grid (same dimension position). NaN where no trial has data.
%   - tGrid: peri-stimulus time grid [s] (column vector), where 0 = stim onset
%   - info: struct with fields
%       - sem: standard error of the mean across trials (same container type as evoked)
%       - trialData: individual trials, with trials as the last dimension (only if 'KeepTrials')
%       - nTrialsPerBin: number of trials contributing to each grid point
%       - trialsUsed: indices of trials that had usable data
%       - tGrid, params
%
% Example (see evokedResponse3D.m):
%   dataTypes.PDI = PDIallSF; dataTypes.CDI = CDIallSF; dataTypes.vUS = VallSF_SG;
%   sfWidth = P.numFramesPerBuffer / P.frameRate;
%   [evoked, tGrid, info] = periStimulusAverage(dataTypes, TD.sfStarts, TD.stimOnsetTimestamps, 'SFWidth', sfWidth, ...
%       'Pre', 5, 'Post', 30, 'Normalize', struct('PDI', 'percent', 'CDI', 'diff', 'vUS', 'diff'));
%   % evoked.PDI{3} is (# x, # y, # z, # grid times)

function [evoked, tGrid, info] = periStimulusAverage(data, sfStarts, stimOnsets, opts)
    arguments
        data
        sfStarts (:, 1) double
        stimOnsets (:, 1) double
        opts.SFWidth (1, 1) double = 0
        opts.Pre (1, 1) double = 5
        opts.Post (1, 1) double = 30
        opts.Dt (1, 1) double = NaN
        opts.BaselineWindow (1, 2) double = [NaN, NaN]
        opts.Normalize = 'percent'
        opts.Method (1, :) char = 'linear'
        opts.MaxGap (1, 1) double = NaN
        opts.TimeDim (1, 1) double = 0
        opts.TrialsToUse (:, 1) double = []
        opts.KeepTrials (1, 1) logical = false
    end

    %% Time grids
    sfCenters = sfStarts + opts.SFWidth/2;
    if numel(sfCenters) < 2 || any(diff(sfCenters) <= 0)
        error('sfStarts must have at least 2 elements and be strictly increasing.')
    end
    sfPeriod = median(diff(sfCenters));

    dt = opts.Dt;
    if isnan(dt), dt = sfPeriod; end
    maxGap = opts.MaxGap;
    if isnan(maxGap), maxGap = 2 * sfPeriod; end

    % Grid is a multiple of dt so it includes t = 0 exactly
    tGrid = (-floor(opts.Pre/dt + 1e-9) : floor(opts.Post/dt + 1e-9)).' * dt;

    baselineWindow = opts.BaselineWindow;
    if any(isnan(baselineWindow)), baselineWindow = [-opts.Pre, 0]; end
    baselineMask = tGrid >= baselineWindow(1) - 1e-9 & tGrid < baselineWindow(2) - 1e-9;

    %% Per-trial sampling info (shared across all data types, since they share the superframe times)
    trialsToUse = opts.TrialsToUse;
    if isempty(trialsToUse), trialsToUse = (1:numel(stimOnsets)).'; end
    if any(trialsToUse < 1 | trialsToUse > numel(stimOnsets) | trialsToUse ~= round(trialsToUse))
        error('TrialsToUse must be valid indices into stimOnsets.')
    end

    nTrials = numel(trialsToUse);
    trial = repmat(struct('idx', [], 'tRel', [], 'bad', true(size(tGrid))), nTrials, 1);
    for k = 1:nTrials
        tRel = sfCenters - stimOnsets(trialsToUse(k)); % Superframe times relative to this trial's stim onset
        idx = find(tRel >= tGrid(1) - maxGap & tRel <= tGrid(end) + maxGap); % Include superframes just outside the grid so edge points can be interpolated
        if numel(idx) < 2, continue, end % Not enough data for this trial; leave it all 'bad'
        tRelSel = tRel(idx);

        % Grid points to exclude: outside the sampled range, or in a gap of missing superframes
        bin = discretize(tGrid, tRelSel);
        bin(bin == numel(tRelSel)) = numel(tRelSel) - 1;
        bad = isnan(bin);
        gap = inf(size(bin));
        gap(~bad) = tRelSel(bin(~bad) + 1) - tRelSel(bin(~bad));
        bad = bad | gap > maxGap;

        trial(k).idx = idx;
        trial(k).tRel = tRelSel;
        trial(k).bad = bad;
    end

    ctx = struct('nSF', numel(sfStarts), 'tGrid', tGrid, 'trial', trial, 'baselineMask', baselineMask, ...
        'Normalize', {opts.Normalize}, 'Method', opts.Method, 'TimeDim', opts.TimeDim, 'KeepTrials', opts.KeepTrials);

    %% Process every data type
    [evoked, sem, trialData] = walk(data, '', ctx);

    %% Info
    validTrial = false(nTrials, 1);
    nTrialsPerBin = zeros(numel(tGrid), 1);
    for k = 1:nTrials
        validTrial(k) = any(~trial(k).bad & baselineMask);
        if validTrial(k), nTrialsPerBin = nTrialsPerBin + ~trial(k).bad; end
    end
    info.sem = sem;
    if opts.KeepTrials, info.trialData = trialData; end
    info.nTrialsPerBin = nTrialsPerBin;
    info.trialsUsed = trialsToUse(validTrial);
    info.tGrid = tGrid;
    info.params = opts;
    info.params.Dt = dt;
    info.params.MaxGap = maxGap;
    info.params.BaselineWindow = baselineWindow;
end

%% Recursively apply to arrays inside cells / structs
function [a, s, t] = walk(x, name, ctx)
    if isnumeric(x) || islogical(x)
        if isempty(x)
            a = []; s = []; t = [];
        else
            [a, s, t] = processArray(x, name, ctx);
        end
    elseif iscell(x)
        a = cell(size(x)); s = a; t = a;
        for k = 1:numel(x)
            [a{k}, s{k}, t{k}] = walk(x{k}, name, ctx);
        end
    elseif isstruct(x) && isscalar(x)
        a = struct(); s = struct(); t = struct();
        for f = fieldnames(x).'
            fn = f{1};
            [ak, sk, tk] = walk(x.(fn), fn, ctx);
            a.(fn) = ak; s.(fn) = sk; t.(fn) = tk;
        end
    else
        error('Unsupported data type: %s. Use numeric arrays, cells, or scalar structs.', class(x))
    end
end

%% Peri-stimulus average of one numeric array (time along dimension td)
function [avg, sem, trl] = processArray(x, name, ctx)
    mode = resolveNormalize(ctx.Normalize, name);
    if ~strcmp(mode, 'none') && ~any(ctx.baselineMask)
        error('BaselineWindow contains no grid points; cannot normalize ''%s''.', name)
    end

    td = ctx.TimeDim;
    if td == 0, td = find(size(x) == ctx.nSF, 1, 'last'); end
    if isempty(td) || size(x, td) ~= ctx.nSF
        error('Data ''%s'' (size %s) has no time dimension of length %d.', name, mat2str(size(x)), ctx.nSF)
    end
    nd = max(ndims(x), td);
    if ~isfloat(x), x = single(x); end
    X = moveDim(x, td, nd); % Time last
    spatialSz = size(X, 1:nd-1);
    nVox = prod(spatialSz);
    X = reshape(X, nVox, ctx.nSF);

    nGrid = numel(ctx.tGrid);
    nTrials = numel(ctx.trial);
    S = zeros(nGrid, nVox);
    S2 = zeros(nGrid, nVox);
    N = zeros(nGrid, nVox);
    if ctx.KeepTrials, T = nan(nGrid, nVox, nTrials, 'single'); end

    for k = 1:nTrials
        tr = ctx.trial(k);
        if isempty(tr.idx), continue, end

        Y = interp1(tr.tRel, X(:, tr.idx).', ctx.tGrid, ctx.Method); % (grid times x voxels)
        Y(tr.bad, :) = NaN;

        % Normalize to this trial's baseline
        if ~strcmp(mode, 'none')
            base = mean(Y(ctx.baselineMask, :), 1, 'omitnan');
            switch mode
                case 'percent', Y = 100 * (Y - base) ./ abs(base);
                case 'diff',    Y = Y - base;
                case 'ratio',   Y = Y ./ base;
                case 'zscore',  Y = (Y - base) ./ std(Y(ctx.baselineMask, :), 0, 1, 'omitnan');
            end
        end

        valid = ~isnan(Y);
        Y0 = double(Y); Y0(~valid) = 0;
        S = S + Y0;
        S2 = S2 + Y0.^2;
        N = N + valid;
        if ctx.KeepTrials, T(:, :, k) = single(Y); end
    end

    avg = S ./ N; % NaN where no trial contributes
    variance = (S2 - S.^2 ./ N) ./ (N - 1);
    sem = sqrt(max(variance, 0) ./ N);
    sem(N < 2) = NaN;

    % Back to the original array layout: (spatial..., grid times) -> time grid in dimension td
    avg = moveDim(reshape(avg.', [spatialSz, nGrid]), nd, td);
    sem = moveDim(reshape(sem.', [spatialSz, nGrid]), nd, td);
    if ctx.KeepTrials
        trl = moveDim(reshape(permute(T, [2, 1, 3]), [spatialSz, nGrid, nTrials]), nd, td);
    else
        trl = [];
    end
end

%% Choose the normalization for a data type name
function mode = resolveNormalize(norm, name)
    if isstruct(norm)
        if ~isempty(name) && isfield(norm, name)
            mode = norm.(name);
        elseif isfield(norm, 'default')
            mode = norm.default;
        else
            mode = 'percent';
        end
    else
        mode = norm;
    end
    mode = validatestring(char(mode), {'percent', 'diff', 'ratio', 'zscore', 'none'});
end

%% Move dimension 'from' of x to position 'to' (like numpy's moveaxis)
function x = moveDim(x, from, to)
    order = 1:max(ndims(x), max(from, to));
    order(from) = [];
    order = [order(1:to-1), from, order(to:end)];
    x = permute(x, order);
end
