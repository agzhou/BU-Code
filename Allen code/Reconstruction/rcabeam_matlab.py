"""Adapter between MATLAB (pyenv, out-of-process) and rcabeam's OpwWorkspace.

MATLAB (R2026b+) hands complex single gpuArrays over with pycupyarray and takes the result back with gpuArray, so
the channel data and the beamformed volumes stay on the GPU. This module only handles what sits between the two:
  - building and holding an OpwWorkspace (allocated once, reused for every chunk of frames),
  - array layout (MATLAB is column-major, the workspace wants C-contiguous),
  - copying the result out of the workspace's borrowed output buffer,
  - separate RC and CR volumes (the workspace itself only returns their coherent sum).

Used by RCA_DAS_rcabeam.m. All vectors arrive from MATLAB as Python lists of floats, and scalars as floats.

Array conventions (n = number of, s = samples, c = channels, a = angles, t = frames):
  channel data   MATLAB [nt, na, nc, ns]  ->  workspace C-contiguous (ns, nc, na, nt). The same bytes, so no copy
                 if pycupyarray keeps the MATLAB memory order.
  volume         workspace (nx, nz, ny, nt) -> returned C-contiguous (nt, nz, ny, nx), i.e. MATLAB [nx, ny, nz, nt]
Whether pycupyarray/gpuArray carry the MATLAB indices or the raw memory order across is not documented. The layout
helpers accept both for the input, and MATLAB-side rcabeam_from_py.m resolves both for the output. test_rcabeam_layout.m
checks them against known values.
"""

from __future__ import annotations

import time

import numpy as np


def _vector(v) -> np.ndarray:
    return np.asarray(v, dtype=np.float64).ravel()


def _fnumber(v):
    """None (full aperture) for None or a non-positive value, which is how MATLAB says 'off'."""
    if v is None:
        return None
    v = float(np.asarray(v, dtype=np.float64).ravel()[0]) if np.size(v) else 0.0
    return v if v > 0 else None


def to_workspace_layout(xp, a, shape):
    """Return `a` as a C-contiguous complex64 array of exactly `shape` = (ns, nc, na, nt).

    `a` is either the raw-memory view of the MATLAB array (shape == `shape`) or the index-preserving view of the
    MATLAB array [nt, na, nc, ns] (shape == `shape` reversed). xp is numpy or cupy. Copies only if needed.
    """
    shape = tuple(int(n) for n in shape)
    rev = shape[::-1]
    if a.ndim != 4:
        raise ValueError(f"Expected a 4-D array, got shape {a.shape}")
    if a.dtype != xp.complex64:
        raise ValueError(f"Expected complex64 (MATLAB complex single), got {a.dtype}")
    if shape == rev:
        # Palindromic shape: the dimensions can't tell the two views apart, so use the memory order
        flip = bool(a.shape == shape and a.flags.f_contiguous and not a.flags.c_contiguous)
    elif a.shape == rev:
        flip = True
    elif a.shape == shape:
        flip = False
    else:
        raise ValueError(f"Array shape {a.shape} is neither {shape} nor MATLAB's {rev}")
    if flip:
        a = a.T
    return xp.ascontiguousarray(a)


def to_matlab_xyz(xp, vol):
    """(nx, nz, ny, nt) -> C-contiguous (nt, nz, ny, nx), which is the memory of a MATLAB [nx, ny, nz, nt] array.

    Always a copy, which also frees the caller from the workspace's borrowed output buffer.
    """
    return xp.ascontiguousarray(vol.transpose(3, 1, 2, 0))


class Session:
    """One OpwWorkspace plus the bookkeeping MATLAB needs. Create once per acquisition geometry and reuse.

    shape: (ns, nc, na, nt)        samples per acquisition, channels, angles per orientation, frames per chunk
    angles: na steering angles (rad), shared by RC and CR
    t_start: na times (s) of the first sample, shared by RC and CR. The IQ must have been demodulated to absolute
             time: iq[n] = analytic(rf)[n] * exp(-2j*pi*f_demod*(t_start + n/fs))
    gx, gz, gy: voxel coordinates (m). Output is (nx, nz, ny, nt) in the workspace
    x_el, y_el: element positions (m) of the RC channels (along x) and the CR channels (along y)
    fnumber: receive f-number, <= 0 for the full aperture
    """

    def __init__(self, shape, angles, t_start, gx, gz, gy, x_el, y_el, c, fs, f_demod, fnumber=0.0,
                 iq_storage="float32"):
        import cupy as cp
        from rcabeam.gpu import OpwWorkspace
        from rcabeam.sim import RCAGeometry

        self._cp = cp
        self.shape = tuple(int(round(n)) for n in _vector(shape))
        geom = RCAGeometry(
            x_el=_vector(x_el), y_el=_vector(y_el), c=float(c), fs=float(fs), f_demod=float(f_demod),
            fnumber=_fnumber(fnumber),
        )
        grid = (_vector(gx), _vector(gz), _vector(gy))
        self._ws = OpwWorkspace(
            self.shape, _vector(angles), _vector(t_start), grid, geom, iq_storage=str(iq_storage),
        )
        self.grid_shape = tuple(len(g) for g in grid)  # (nx, nz, ny)
        self._zeros = None
        self.times = {}

    # Pack + reconstruct + copy out. The workspace output is borrowed, so it is transposed into a new array
    def _run(self, rc, cr):
        cp = self._cp
        t0 = time.perf_counter()
        self._ws.pack(rc, cr)
        t1 = time.perf_counter()
        vol = self._ws.reconstruct()
        t2 = time.perf_counter()
        out = to_matlab_xyz(cp, vol)
        cp.cuda.Stream.null.synchronize()
        t3 = time.perf_counter()
        self.times = {"pack": t1 - t0, "reconstruct": t2 - t1, "copy_out": t3 - t2}
        return out

    def _zero(self):
        if self._zeros is None:
            self._zeros = self._cp.zeros(self.shape, self._cp.complex64)
        return self._zeros

    def beamform_sum(self, rc, cr):
        """Coherent RC + CR sum (OPW), one pass. rc/cr are MATLAB [nt, na, nc, ns] complex single."""
        cp = self._cp
        return self._run(to_workspace_layout(cp, rc, self.shape), to_workspace_layout(cp, cr, self.shape))

    def beamform_rc(self, rc):
        """RC volume alone (zeros in place of CR). The sum is linear, so beamform_rc + beamform_cr == beamform_sum."""
        cp = self._cp
        return self._run(to_workspace_layout(cp, rc, self.shape), self._zero())

    def beamform_cr(self, cr):
        """CR volume alone (zeros in place of RC)."""
        cp = self._cp
        return self._run(self._zero(), to_workspace_layout(cp, cr, self.shape))

    def report(self) -> str:
        """Seconds spent in the last call, as text."""
        return ", ".join(f"{k} {v * 1e3:.1f} ms" for k, v in self.times.items())

    def release(self) -> None:
        """Drop the workspace and return its GPU memory to the driver."""
        self._ws = None
        self._zeros = None
        pool = self._cp.get_default_memory_pool()
        pool.free_all_blocks()


# Helpers for test_rcabeam_layout.m and for checking the install from MATLAB

def gpu_info() -> str:
    import cupy as cp

    props = cp.cuda.runtime.getDeviceProperties(cp.cuda.Device().id)
    free, total = cp.cuda.runtime.memGetInfo()
    try:
        import rcabeam

        rc = f"rcabeam {getattr(rcabeam, '__version__', '(no version)')} from {rcabeam.__file__}"
    except ImportError as e:
        rc = f"rcabeam NOT importable: {e}"
    name = props["name"].decode() if isinstance(props["name"], bytes) else props["name"]
    return f"{name}, {free / 2**30:.1f}/{total / 2**30:.1f} GiB free, cupy {cp.__version__}, {rc}"


def describe(a) -> str:
    return (f"shape {tuple(a.shape)}, strides {tuple(a.strides)}, dtype {a.dtype}, "
            f"c_contiguous {a.flags.c_contiguous}, f_contiguous {a.flags.f_contiguous}")


def echo_workspace_layout(a, shape):
    """Input-layout test: what the workspace would be handed for MATLAB array `a`, as a C-contiguous array."""
    import cupy as cp

    return to_workspace_layout(cp, a, [int(round(n)) for n in _vector(shape)])


def test_volume(nx, nz, ny, nt):
    """Output-layout test: a workspace-style (nx, nz, ny, nt) volume with value (x+1) + 1e3*(z+1) + 1e6*(y+1) + 1j*(t+1),
    returned the way beamform_* returns volumes. In MATLAB, out(x, y, z, t) must equal that for the 1-based x, y, z, t."""
    import cupy as cp

    nx, nz, ny, nt = (int(round(n)) for n in (nx, nz, ny, nt))
    x, z, y, t = np.meshgrid(np.arange(nx), np.arange(nz), np.arange(ny), np.arange(nt), indexing="ij")
    vol = ((x + 1) + 1e3 * (z + 1) + 1e6 * (y + 1)) + 1j * (t + 1)
    return to_matlab_xyz(cp, cp.asarray(vol.astype(np.complex64)))
