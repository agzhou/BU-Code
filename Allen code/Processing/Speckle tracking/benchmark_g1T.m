% benchmark_g1T compares g1T.m (direct per-lag summation) against
% g1T_fft.m (FFT/Wiener-Khinchin based) for both speed and numerical
% accuracy, across a range of spatial dimensionalities (1D/2D/3D),
% voxel counts, and lag ranges (np relative to nf).
%
% Input (both optional):
%   nReps            - # of timed repetitions per case; the median is
%                       reported (default 3)
%   includeLargeCase - if true, also runs one case sized close to a real
%                       vUS volume (~620k voxels, 1000 frames). This can
%                       take a while for g1T (the direct method), which
%                       is the point. Default false.
%
% Output:
%   resultsTable - table with per-case timing (seconds), speedup, and
%                  accuracy (max absolute error, RMS error, relative
%                  error %) of g1T_fft vs. g1T
%
% Test data: synthetic complex speckle with an AR(1) temporal
% correlation (autocorrelation ~ rho^|tau|) along the frame dimension, so
% the comparison is done on data with realistic non-trivial correlation
% structure rather than pure white noise.

function resultsTable = benchmark_g1T(nReps, includeLargeCase)
    if nargin < 1 || isempty(nReps)
        nReps = 3;
    end
    if nargin < 2 || isempty(includeLargeCase)
        includeLargeCase = false;
    end

    rng(1); % reproducible test data
    rho = 0.985; % AR(1) temporal correlation coefficient for synthetic data

    cases = struct('name', {}, 'spatialSize', {}, 'nf', {}, 'np', {});
    cases(end+1) = struct('name', '1D, short lag range', 'spatialSize', 5000,          'nf', 1000, 'np', 50);
    cases(end+1) = struct('name', '1D, full lag range',  'spatialSize', 2000,          'nf', 1000, 'np', 1000);
    cases(end+1) = struct('name', '2D, short lag range', 'spatialSize', [100, 100],    'nf', 1000, 'np', 50);
    cases(end+1) = struct('name', '2D, full lag range',  'spatialSize', [40, 40],      'nf', 1000, 'np', 1000);
    cases(end+1) = struct('name', '3D, short lag range', 'spatialSize', [40, 40, 30],  'nf', 1000, 'np', 60);
    cases(end+1) = struct('name', '3D, full lag range',  'spatialSize', [20, 20, 15],  'nf', 500,  'np', 500);

    if includeLargeCase
        cases(end+1) = struct('name', '3D, near real-scale volume', 'spatialSize', [88, 88, 80], 'nf', 1000, 'np', 100);
    end

    nCases = numel(cases);
    name = strings(nCases, 1);
    spatialSizeStr = strings(nCases, 1);
    nf = zeros(nCases, 1);
    np = zeros(nCases, 1);
    tDirect = zeros(nCases, 1);
    tFFT = zeros(nCases, 1);
    maxAbsErr = zeros(nCases, 1);
    rmsErr = zeros(nCases, 1);
    relErrPercent = zeros(nCases, 1);

    for ci = 1:nCases
        c = cases(ci);
        fprintf('Case %d/%d: %s (spatial = %s, nf = %d, np = %d)\n', ...
            ci, nCases, c.name, mat2str(c.spatialSize), c.nf, c.np);

        data = genSyntheticSpeckle([c.spatialSize, c.nf], rho);

        g1_direct = g1T(data, c.np); % warm-up (JIT, memory paging) before timing
        ts = zeros(nReps, 1);
        for r = 1:nReps
            t0 = tic;
            g1_direct = g1T(data, c.np);
            ts(r) = toc(t0);
        end
        tDirect(ci) = median(ts);

        g1_fft = g1T_fft(data, c.np); % warm-up
        ts = zeros(nReps, 1);
        for r = 1:nReps
            t0 = tic;
            g1_fft = g1T_fft(data, c.np);
            ts(r) = toc(t0);
        end
        tFFT(ci) = median(ts);

        diffAbs = abs(g1_fft - g1_direct);
        maxAbsErr(ci) = max(diffAbs(:));
        rmsErr(ci) = sqrt(mean(diffAbs(:).^2));
        denomRms = sqrt(mean(abs(g1_direct(:)).^2)) + eps;
        relErrPercent(ci) = 100 * rmsErr(ci) / denomRms;

        name(ci) = c.name;
        spatialSizeStr(ci) = string(mat2str(c.spatialSize));
        nf(ci) = c.nf;
        np(ci) = c.np;
    end

    speedup = tDirect ./ tFFT;

    resultsTable = table(name, spatialSizeStr, nf, np, tDirect, tFFT, speedup, maxAbsErr, rmsErr, relErrPercent, ...
        'VariableNames', {'Case', 'SpatialSize', 'nf', 'np', 'g1T_seconds', 'g1T_fft_seconds', 'Speedup', 'MaxAbsErr', 'RMSErr', 'RelErrPercent'});

    disp(resultsTable)

    figure;
    bar(categorical(cellstr(name), cellstr(name)), [tDirect, tFFT]);
    set(gca, 'YScale', 'log');
    ylabel('Time [s] (log scale)');
    legend({'g1T (direct)', 'g1T_{fft}'}, 'Location', 'best');
    title('g1T vs g1T_{fft}: runtime by test case');
    xtickangle(20);
end

function data = genSyntheticSpeckle(dataSize, rho)
    % Complex speckle with AR(1) temporal correlation (~rho^|tau|) along
    % the last dimension, for use as realistic-ish comparison input.
    frameDim = numel(dataSize);
    w = (randn(dataSize) + 1i * randn(dataSize)) / sqrt(2);
    data = filter(sqrt(1 - rho^2), [1, -rho], w, [], frameDim);
end
