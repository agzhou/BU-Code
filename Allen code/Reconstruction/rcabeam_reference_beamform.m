%% Description

% MATLAB (CPU) port of rcabeam's per-channel delay-and-sum for one transmit orientation. It follows delay_sample in
% rcabeam's core.cu and the OPW path of ensemble.cu line by line, so that RCA_DAS_rcabeam's data preparation (angle and
% group mapping, t_start, demodulation phase) can be tested without a GPU or Python, and so that the CUDA kernels can be
% compared against a plain description of the same model. Slow, meant for small test volumes.

% The model, per voxel (x, z, y), channel e and angle a (config 0 = RC, 1 = CR):
%   u_tx, v_rx = (y, x) for RC, (x, y) for CR            TX is steered along u, the receive elements lie along v
%   dv     = v_rx - el(e)
%   active = fNumber <= 0 or |dv| <= z / (2 fNumber)      hard receive aperture
%   tau    = (z cos(th) + u_tx sin(th)) / c + sqrt(z^2 + dv^2) / c      two-way plane wave time of flight
%   n      = (tau - t0(a)) * fs                           0-based, fractional sample index
%   value  = linear interpolation of IQ between floor(n) and floor(n)+1 (0 outside 0 <= floor(n) <= ns-2) * exp(+2j pi fc tau)
% summed over e and a.

% Input: D       IQ [nT, na, nc, ns] (complex, MATLAB layout of the channel data given to the rcabeam adapter)
%        config  0 (RC) or 1 (CR)
%        angles  1 x na, t0  1 x na (s), gx, gz, gy  voxel coordinates in m, el  1 x nc element positions in m
%        c, fs, fc, fNumber
% Output: vol [nx, ny, nz, nT] single

function vol = rcabeam_reference_beamform(D, config, angles, t0, gx, gz, gy, el, c, fs, fc, fNumber)

    D = gather(D);
    [nT, na, nc, ns] = size(D);
    nz = numel(gz);
    if config == 0
        uv = gy(:); vv = gx(:); % TX steered along y, receive along x
    else
        uv = gx(:); vv = gy(:);
    end
    nu = numel(uv); nv = numel(vv);
    el = reshape(el, 1, []);
    angles = reshape(angles, 1, []);

    dv = vv - el; % (nv, nc)
    [~, C, ~, A] = ndgrid(1:nv, 1:nc, 1:nu, 1:na);
    vol = complex(zeros(nv, nu, nz, nT));

    for iz = 1:nz
        z = gz(iz);
        rxT = sqrt(z^2 + dv.^2) / c; % (nv, nc)
        if fNumber > 0
            active = abs(dv) <= z / (2 * fNumber);
        else
            active = true(nv, nc);
        end
        txT = (z * cos(angles) + uv * sin(angles)) / c; % (nu, na)
        tau = reshape(rxT, nv, nc, 1, 1) + reshape(txT, 1, 1, nu, na); % (nv, nc, nu, na)
        s = (tau - reshape(t0, 1, 1, 1, na)) * fs;
        i0 = floor(s);
        w = s - i0;
        ok = i0 >= 0 & i0 <= ns - 2 & repmat(active, 1, 1, nu, na);
        ph = exp(1j * 2*pi * fc * tau);
        idx0 = sub2ind([na, nc, ns], A(ok), C(ok), i0(ok) + 1);
        wok = w(ok);
        phok = ph(ok);
        for f = 1:nT
            Df = double(reshape(D(f, :, :, :), na, nc, ns));
            val = ((1 - wok) .* Df(idx0) + wok .* Df(idx0 + na * nc)) .* phok;
            contrib = complex(zeros(nv, nc, nu, na));
            contrib(ok) = val;
            vol(:, :, iz, f) = reshape(sum(sum(contrib, 4), 2), nv, nu);
        end
    end

    if config == 0 % (x, y, z, t) already
        vol = single(vol);
    else           % (y, x, z, t)
        vol = single(permute(vol, [2, 1, 3, 4]));
    end

end
