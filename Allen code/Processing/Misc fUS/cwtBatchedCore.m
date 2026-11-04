function S = cwtBatchedCore(X, P, outputType)
% CWTBATCHEDCORE  CWT of many real series at once (same math as cwt(x,'FilterBank',fb)).
%   X : [nT x nSeries] real; time down the columns. Can be a gpuArray.
%   P : struct from cwtPrepare (psiT must be same class / device as X)
%   outputType : "complex" (default) | "abs" | "abs2"
%   S : [nT x nSeries x nFreq]  (time-major: avoids a slow strided write; permute afterwards if needed)
if nargin < 3, outputType = "complex"; end
[nT, nS] = size(X);
pad = P.pad;
if pad > 0
    if P.reflect
        Xp = [flip(X(1:pad, :), 1); X; flip(X(end-pad+1:end, :), 1)];
    else
        Xp = [zeros(pad, nS, 'like', X); X; zeros(pad, nS, 'like', X)];
    end
else
    Xp = X;
end
Xf = fft(Xp, [], 1);
nF = size(P.psiT, 2);
if outputType == "complex"
    S = complex(zeros(nT, nS, nF, 'like', X));
else
    S = zeros(nT, nS, nF, 'like', real(X));
end
rows = pad+1:pad+nT;
for k = 1:nF
    C = ifft(Xf .* P.psiT(:, k), [], 1);
    C = C(rows, :);
    switch outputType
        case "abs",  C = abs(C);
        case "abs2", C = real(C).^2 + imag(C).^2;
    end
    S(:, :, k) = C;
end
end
