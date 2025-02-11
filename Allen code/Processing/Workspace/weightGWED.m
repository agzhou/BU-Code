% Calculate the NLM weight based on the Gaussian-weighted Euclidean distance between two matrices

function [weight] = weightGWED(mati, matj, h)
    veci = mati(:);
%     mu = mean(mati, 'all');
    mu = mean(veci);
%     gaus = @(vecj) exp(-abs(vecj - mu) .^ 2 ./ h^2)
    vecj = matj(:);
%     weight = sum((veci - vecj).^2 .* exp(-abs(vecj - mu).^2 ./ h^2)); % I think the abs is to account for complex values
    GWED = -sum((veci - vecj).^2 .* exp(-abs(vecj - mu).^2 ./ h^2));
    weight = exp(GWED ./ h^2); % I think the abs is to account for complex values

end