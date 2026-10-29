%% Description

% Sets up MATLAB's Python interface for the rcabeam adapter (rcabeam_matlab.py) and returns the imported module.
% Needs MATLAB R2026b or newer (pycupyarray), the Parallel Computing Toolbox, and a Python environment with cupy
% (CUDA 13) and rcabeam installed. Python runs out-of-process, which pycupyarray and gpuArray(<cupy array>) require.

% Input: pythonExe, optional full path to the python.exe of that environment. If empty, the RCABEAM_PYTHON environment
%        variable is used, and if that is empty too, whatever pyenv already points to.
% Output: mod, the Python module rcabeam_matlab (call its functions as mod.Session(...), mod.gpu_info(), ...)

% Notes:
%           - pyenv can only be changed before Python is first used in this MATLAB session (or after
%             terminate(pyenv)). If Python was already started in-process, restart MATLAB.

function mod = rcabeam_setup(pythonExe)

    if nargin < 1
        pythonExe = '';
    end
    if isempty(pythonExe)
        pythonExe = getenv('RCABEAM_PYTHON');
    end

    tooOld = true;
    try %#ok<TRYNC> isMATLABReleaseOlderThan errors on a release newer than the running one
        tooOld = isMATLABReleaseOlderThan("R2026b");
    end
    if tooOld
        error('rcabeam_setup:release', 'pycupyarray needs MATLAB R2026b or newer (this is %s).', matlabRelease.Release)
    end

    pe = pyenv;
    if string(pe.Status) == "NotLoaded"
        if isempty(pythonExe)
            pyenv(ExecutionMode="OutOfProcess");
        else
            pyenv(Version=pythonExe, ExecutionMode="OutOfProcess");
        end
    else
        if string(pe.ExecutionMode) ~= "OutOfProcess"
            error('rcabeam_setup:inProcess', ['Python is already loaded in-process. pycupyarray needs out-of-process ', ...
                'mode: restart MATLAB and call pyenv(ExecutionMode="OutOfProcess") before using Python.'])
        end
        if ~isempty(pythonExe) && ~strcmpi(char(pe.Executable), pythonExe)
            warning('rcabeam_setup:pyenv', 'Python is already running from %s, not %s. Use terminate(pyenv) to switch.', ...
                char(pe.Executable), pythonExe)
        end
    end

    % Make rcabeam_matlab.py importable
    adapterDir = fileparts(mfilename('fullpath'));
    paths = cellfun(@char, cell(py.sys.path), 'UniformOutput', false);
    if ~any(strcmp(paths, adapterDir))
        insert(py.sys.path, int32(0), adapterDir);
    end
    mod = py.importlib.import_module('rcabeam_matlab');

end
