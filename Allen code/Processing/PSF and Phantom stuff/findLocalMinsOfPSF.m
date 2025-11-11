function [lowerBound_ind, upperBound_ind] = findLocalMinsOfPSF(y)
    [~, mind] = max(y);

    yleft = y(1:mind);
    yright = y(mind:end);
%     xleft = x(1:mind);
%     xright = x(mind:end);

    dleft = diff(yleft);
    dright = diff(yright);

    lowerBound_ind = find(dleft <= 0, 1, 'last') + 1; % lower bound, add 1 bc diff takes away a value
    upperBound_ind = mind + find(dright >= 0, 1, 'first') - 1; % upper bound
    

end
