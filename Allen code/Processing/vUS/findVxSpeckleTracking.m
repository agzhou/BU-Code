% Description: Estimate the transverse group velocity v_xgp from spatial
% (lateral) speckle tracking -- cross-correlating a reference IQ line
% against a lagged IQ line and converting the sub-pixel correlation-peak
% displacement into a velocity -- using an ADAPTIVE tracking lag chosen
% from the already-available temporal g1(tau) magnitude curve.
%
% Why this exists: v_xgp and p (the intra-voxel axial velocity-spread
% parameter) are NOT jointly identifiable from a single voxel's g1(tau)
% curve alone -- they only ever enter the model through the same combined
% decay coefficient (v_xgp^2/(4*sigma_x^2) + p^2*(v_zgp*k0)^2). A Monte
% Carlo comparison (synthetic g1 built from the vUS model + realistic
% finite-sample estimator noise, matching Eqs. 2 and 15 of Tang et al.
% 2020) showed a direct joint fit (g1vUS2D_residJac.m) picks up a growing
% v_xgp bias as true p increases (~12mm/s bias at p=1 for a true v_xgp of
% 3mm/s). Spatial speckle tracking is a structurally different
% measurement -- it reads a coherent spatial SHIFT, not a decorrelation
% RATE -- so it doesn't inherit that confound: bias stayed near zero
% across the full p in [0,1] range in that same comparison.
%
% Why the lag needs to be adaptive: at a FIXED tracking lag, high p
% decorrelates the signal before that lag is reached, and tracking
% degrades from "unbiased but noisy" to "unreliable" -- in the same test,
% a fixed 20-frame (4 ms) lag gave 25.8 mm/s RMSE and 11.3% outright
% tracking failures at p=1 (rather than a small bias, a POOR correlation
% peak is just as likely to be noise as signal, so the peak-finding itself
% becomes unreliable). This function instead picks the LONGEST lag (for
% the best achievable displacement precision) at which the observed
% |g1(tau)| curve is still above a coherence threshold -- no knowledge of
% the true p is used, just the same g1(tau) curve already computed
% elsewhere in the pipeline (e.g. for findVzPhaseDiff.m). Verified in the
% same Monte Carlo: RMSE at p=1 dropped to 2.8 mm/s with 0% failures,
% while staying comparable to (occasionally better than) the fixed-lag
% version everywhere from p=0 to p=0.7.
%
% The correlation search window is derived ONLY from a physiologically
% plausible max-velocity bound (default 50 mm/s, matching the lb/ub used
% throughout this codebase's fits), centered at zero shift -- NOT centered
% on the true shift, since a real deployment doesn't know that in advance.
%
% Scope/limitations: operates on ONE lateral line (single depth row) at a
% time, not vectorized across many voxels/lines -- extending to a full 2D
% image would need looping this function per depth row (or per voxel with
% a small 2D patch) and has not been tested. 1D lateral tracking only; a
% real 2D block-matching implementation would also track axial (z) shift
% for extra precision, which this does not attempt (v_zgp is assumed to
% already come from findVzPhaseDiff.m, which does not have this
% v_xgp/p confound in the first place).
%
% Inputs:
%   IQline: [nx, nFrames] complex IQ samples along one lateral line across
%       time, frame 1 = the reference frame (t=0). nFrames must be >= Lmax+1.
%   g1curve: complex (or real magnitude) g1(tau) for this same voxel,
%       [1,nTauG1] or [nTauG1,1], index 1 = tau=0 (g1T.m/findVzPhaseDiff.m
%       convention). Used only to choose the tracking lag -- not tracked
%       itself.
%   dx: lateral pixel spacing of IQline [m]
%   PP: struct with PP.frameRate [Hz] (sIQ frame rate)
%   cohThreshold: (optional, default 0.35) minimum smoothed |g1(tau)| to
%       trust a candidate lag
%   Lmin, Lmax: (optional, default 3, min(60, nTauG1-1)) bounds on the
%       candidate lag, in frames
%   vxMaxExpected: (optional, default 50e-3 m/s) max plausible |v_xgp|,
%       used only to size the (zero-centered) correlation search window
%
% Outputs:
%   vxEst: estimated transverse group velocity [m/s]. NaN if tracking failed.
%   Lused: the tracking lag actually used [frames]
%   failed: true if the correlation peak fell on the search-window edge
%       (unreliable -- vxEst is NaN in this case)
function [vxEst, Lused, failed] = findVxSpeckleTracking(IQline, g1curve, dx, PP, cohThreshold, Lmin, Lmax, vxMaxExpected)
    if nargin < 5 || isempty(cohThreshold), cohThreshold = 0.35; end
    if nargin < 6 || isempty(Lmin), Lmin = 3; end
    if nargin < 8 || isempty(vxMaxExpected), vxMaxExpected = 50e-3; end
    nTauG1 = numel(g1curve);
    if nargin < 7 || isempty(Lmax), Lmax = min(60, nTauG1-1); end

    dt = 1/PP.frameRate;
    [nx, nFrames] = size(IQline);
    Lmax = min([Lmax, nTauG1-1, nFrames-1]);

    % ---- 1. Choose the tracking lag from the already-available g1(tau) curve ----
    ag1 = abs(g1curve(:)).';
    ag1_smooth = movmean(ag1, 5);
    idxRange = (Lmin+1):(Lmax+1); % index 1 corresponds to tau=0
    idxRange = idxRange(idxRange <= numel(ag1_smooth));
    aboveThresh = ag1_smooth(idxRange) >= cohThreshold;
    lastGoodIdx = find(aboveThresh, 1, 'last');
    if isempty(lastGoodIdx)
        Lused = Lmin; % nothing meets threshold; fall back to the shortest allowed lag
    else
        Lused = idxRange(lastGoodIdx) - 1;
    end
    elapsedTime = Lused*dt;

    % ---- 2. Cross-correlate the reference and lagged IQ lines ----
    frame0 = IQline(:,1).';
    frameL = IQline(:,Lused+1).';

    searchRange = max(3, ceil(vxMaxExpected*elapsedTime/dx));
    searchRange = min(searchRange, floor(nx/2)-1);
    lags = -searchRange:searchRange;
    cc = zeros(size(lags));
    for li = 1:numel(lags)
        L = lags(li);
        if L >= 0
            a = frame0(1:end-L); b = frameL(1+L:end);
        else
            a = frame0(1-L:end); b = frameL(1:end+L);
        end
        if isempty(a)
            cc(li) = 0; continue
        end
        cc(li) = abs(mean(conj(a).*b));
    end

    [~, imax] = max(cc);
    failed = (imax==1 || imax==numel(lags));
    if failed
        vxEst = NaN;
        return
    end

    % 3-point parabolic sub-pixel refinement of the correlation peak
    y1 = cc(imax-1); y2 = cc(imax); y3 = cc(imax+1);
    delta = 0.5*(y1-y3)/(y1-2*y2+y3+eps);
    peakLagSamples = lags(imax) + delta;
    estShift = peakLagSamples*dx;
    vxEst = estShift / elapsedTime;
end
