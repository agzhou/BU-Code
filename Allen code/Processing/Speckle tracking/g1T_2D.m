% g1T calculation in 2D
% Input: IQ (coherently summed across angles)
% Output: g1

% Notes: 
% - MIGHT NOT WORK WITH MULTIPLE SUPERFRAMES

function [g1] = g1T_2D(IQ_coherent_sum)
    tstart = clock;
    [zp, xp, nf] = size(IQ_coherent_sum);

    g1 = zeros(zp, xp, nf); % g1 is calculated for each pixel
                                % z, x, tau step
    numer = zeros(zp, xp, nf);
    denom = mean((conj(IQ_coherent_sum) .* IQ_coherent_sum), 3); % temporal (frame) average
    for f = 1:nf % go through each subframe to get the g1 at each tau step
        numer = mean(conj(IQ_coherent_sum(:, :, 1:(nf - f + 1))) .* IQ_coherent_sum(:, :, f:end), 3);
        g1(:, :, f) = numer ./ denom;
    
    end

    tend = clock;
    disp(strcat("Temporal g_{1} processing done, elapsed time is ", num2str(etime(tend, tstart)), "s"))
end