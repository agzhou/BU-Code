%% Description

% Checks the hand-over of arrays between MATLAB and rcabeam_matlab.py with pycupyarray and gpuArray, using arrays with
% known values. Run this once on the GPU machine (MATLAB R2026b), before trusting RCA_DAS_rcabeam. Neither
% pycupyarray nor gpuArray(<cupy array>) documents whether the MATLAB indices or the raw memory order are carried
% across, so this reports which one happens. The adapter and rcabeam_from_py accept both, so a PASS for either is fine.
% Set pythonExe beforehand, or the RCABEAM_PYTHON environment variable, to the python.exe with cupy and rcabeam.

if exist('pythonExe', 'var') ~= 1, pythonExe = ''; end
mod = rcabeam_setup(pythonExe);
fprintf('%s\n\n', char(mod.gpu_info()))

%% MATLAB -> Python: channel data [nt, na, nc, ns]
shape = [7, 5, 4, 3]; % workspace (ns, nc, na, nt), all different so a mix-up shows
nt = shape(4); na = shape(3); nc = shape(2); ns = shape(1);
v = single(1 : prod(shape));
A = gpuArray(complex(reshape(v, [nt, na, nc, ns]), -reshape(v, [nt, na, nc, ns])));
pyA = pycupyarray(A);
fprintf('pycupyarray gave: %s\n', char(mod.describe(pyA)))

pyW = mod.echo_workspace_layout(pyA, py.list(num2cell(double(shape))));
fprintf('workspace layout: %s\n', char(mod.describe(pyW)))
W = gather(gpuArray(pyW));
A = gather(A);
if isequal(size(W), shape) && isequal(W, permute(A, [4, 3, 2, 1]))
    fprintf('gpuArray(cupy) keeps the Python indices: size %s\n', mat2str(size(W)))
    inOK = true;
elseif isequal(size(W), fliplr(shape)) && isequal(W, A)
    fprintf('gpuArray(cupy) keeps the memory order (reverses the dimensions): size %s\n', mat2str(size(W)))
    inOK = true;
else
    fprintf('Unexpected array of size %s coming back\n', mat2str(size(W)))
    inOK = false;
end
fprintf('Input layout: %s\n\n', passFail(inOK))

%% Python -> MATLAB: volume (nx, nz, ny, nt) -> [nx, ny, nz, nt]
nx = 6; nz = 5; ny = 4; nT = 3;
pyV = mod.test_volume(nx, nz, ny, nT);
fprintf('test volume from Python: %s\n', char(mod.describe(pyV)))
V = gather(rcabeam_from_py(pyV, [nx, ny, nz, nT]));
[X, Y, Z, T] = ndgrid(1:nx, 1:ny, 1:nz, 1:nT);
expected = single(X + 1e3 * Z + 1e6 * Y) + 1j * single(T);
outOK = isequal(V, expected);
fprintf('Output layout: %s\n', passFail(outOK))


function s = passFail(tf)
    if tf
        s = 'PASS';
    else
        s = 'FAIL';
    end
end
