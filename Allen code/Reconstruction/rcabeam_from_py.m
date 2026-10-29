%% Description

% Converts a CuPy volume returned by rcabeam_matlab.Session.beamform_* into a MATLAB gpuArray of size
% [nx, ny, nz, nt]. The adapter returns a C-contiguous (nt, nz, ny, nx) array. Whether gpuArray keeps the Python
% index order (size [nt nz ny nx], needs a permute) or the memory order (already [nx ny nz nt]) is not documented, so
% this accepts both and errors on anything else.

% Input: pyVol, the CuPy array; want, [nx, ny, nz, nt]
% Output: v, gpuArray [nx, ny, nz, nt]

function v = rcabeam_from_py(pyVol, want)

    v = gpuArray(pyVol);
    sz = size(v);
    sz(end+1 : 4) = 1; % trailing singleton dimensions (e.g., one frame)
    if isequal(sz, want)
        return
    elseif isequal(sz, fliplr(want))
        v = permute(v, [4, 3, 2, 1]);
    else
        error('rcabeam_from_py:size', 'Got a %s array from Python, expected %s (or %s)', ...
            mat2str(sz), mat2str(want), mat2str(fliplr(want)))
    end

end
