function L = goodFFTLength(n)
% Smallest integer >= n whose only prime factors are 2, 3, 5 and 7 (sizes cuFFT handles efficiently)
L = ceil(n);
while true
    m = L; % Test integer m
    % for p = [2 3 5 7] % Go through selected prime factors
    for p = [2 3 5] % Go through selected prime factors
        while mod(m, p) == 0 % Keep dividing m by the prime factors until they don't divide anymore
            m = m / p;
        end
    end
    if m == 1 % If m has been divided as many times as possible by the prime factors and is left with 1, it only has those prime factors, so return and choose that number
        return
    end
    L = L + 1;
end
end