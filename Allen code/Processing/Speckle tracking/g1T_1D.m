% Input: IQ (coherently summed across angles)
% Output: g1

% Notes: 
% - MIGHT NOT WORK WITH MULTIPLE SUPERFRAMES

function [g1] = g1T_1D(IQ_coherent_sum)
    tstart = clock;
%     [xp, yp, zp, nsubf, nsupf] = size(IQ_coherent_sum);
    [xp, yp, nsubf] = size(IQ_coherent_sum);

    g1 = zeros(xp, yp, nsubf); % g1 is calculated for each voxel
                                % x, y, z, tau step
    numer = zeros(xp, yp, nsubf);
    denom = mean((conj(IQ_coherent_sum) .* IQ_coherent_sum), 3); % temporal (frame) average
    for f = 1:nsubf % go through each subframe to get the g1 at each tau step
        numer = mean(conj(IQ_coherent_sum(:, :, 1:(nsubf - f + 1))) .* IQ_coherent_sum(:, :, f:end), 3);
        g1(:, :, f) = numer ./ denom;
    
    end

    tend = clock;
    disp(strcat("Temporal g_{1} processing done, elapsed time is ", num2str(etime(tend, tstart)), "s"))
end