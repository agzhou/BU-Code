function [z] = fisherTransform3D(r, ntp)
    z = sqrt(ntp - 3)/2 .* log((1 + r) ./ (1 - r));
end