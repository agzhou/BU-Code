% find fwhm of a 1D data series
% linear interpolation in case the data is too sparse
% assumes the minimum y value is 0 (at least, theoretically)

function [width] = fwhm(x, y)
% x = -10:0.1:10;
% y = -testx .^ 2 + 20;
% plot(x, y)
    [m, mind] = max(y);
    hm = m/2;

    yleft = y(1:mind);
    yright = y(mind:end);
    xleft = x(1:mind);
    xright = x(mind:end);

    leftind_lb = find(yleft <= hm, 1, 'last'); % lower bound
%     leftind_ub = find(yleft >= hm, 1, 'first');
    leftind_ub = leftind_lb + 1;               % upper bound

    rightind_lb = find(yright <= hm, 1, 'first');
    rightind_ub = rightind_lb - 1;
    if yleft(leftind_ub) < hm | yright(rightind_ub) < hm
        warning('Data has weird behavior')
    end

%     leftind_lb
%     leftind_ub
%     rightind_ub
%     rightind_lb

    leftbound = interp1([yleft(leftind_lb), yleft(leftind_ub)], [xleft(leftind_lb), xleft(leftind_ub)], hm);
    rightbound = interp1([yright(rightind_lb), yright(rightind_ub)], [xright(rightind_lb), xright(rightind_ub)], hm);

    width = rightbound - leftbound;

end