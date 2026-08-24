function [CM, SVs, V] = getSVs1D_SVD(IQ)
    [zp, xp, nf] = size(IQ);
    
    CM = reshape(IQ, [zp*xp, nf]); % Covariance matrix
    tic
    %     [U, S, V] = svd(PP); % Already sorted in decreasing order
    [~, S, V] = svd(CM, 'econ'); % Already sorted in decreasing order
    SVs = diag(S);
end