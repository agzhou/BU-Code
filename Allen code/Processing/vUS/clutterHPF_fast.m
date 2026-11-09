% SVD clutter filter + high pass filter in ONE gemm.
%
% Equivalent to getSVs2D -> applySVs2D_accel -> filter(HPF.b, HPF.a, IQf, [], 4), but:
%   - the nf-by-nf eigendecomposition runs on the CPU (it is tiny, and slow on a GPU)
%   - the SVD projector Q and the (causal, linear) HPF are combined into a single nf-by-nf matrix T = Q*Hm,
%     so the volume is only multiplied once and no intermediate clutter-filtered volume is formed
%   - Q is built from whichever of the kept / discarded singular vectors is fewer
%
% Input:  IQ      [x, y, z, frame], complex, CPU or gpuArray
%         svLo, svHi   kept singular value range (indices into the descending sort), as in applySVs2D_accel
%         HPF     struct with b, a (from butter)
% Output: Y       [x, y, z, frame], same class / location as IQ

function Y = clutterHPF_fast(IQ, svLo, svHi, HPF)

    [xp, yp, zp, nf] = size(IQ);
    if svHi > nf
        error('Upper threshold is larger than the number of eigenvalues (frames)')
    end

    PP = reshape(IQ, [], nf); % space x time

    % Covariance in time (pagemtimes with a transpose flag avoids materializing PP')
    CM = pagemtimes(PP, 'ctranspose', PP, 'none');
    CM = double(gather(CM));
    CM = (CM + CM') / 2; % exactly Hermitian, so eig takes the Hermitian path

    [V, ev] = eig(CM, 'vector'); % ascending, real
    [~, order] = sort(real(ev), 'descend');
    V = V(:, order);

    % Clutter-filter projector Q = V_sel*V_sel'. If most components are kept, build it as I - V_rm*V_rm' instead
    nKeep = svHi - svLo + 1;
    if nf - nKeep < nKeep
        V_rm = V(:, [1:svLo-1, svHi+1:nf]);
        Q = eye(nf) - V_rm * V_rm';
    else
        V_sel = V(:, svLo:svHi);
        Q = V_sel * V_sel';
    end

    % HPF as a matrix: row j = response to an impulse at frame j, so (x * Hm) == filter(b, a, x, [], 2) for a row x.
    % filter() has zero initial conditions and is causal, so this is exact
    Hm = filter(HPF.b, HPF.a, eye(nf), [], 2);

    T = Q * Hm; % row vector x -> x*Q (clutter filter) -> *Hm (HPF)
    T = cast(T, underlyingType(PP));
    if isa(PP, 'gpuArray')
        T = gpuArray(T);
    end

    Y = reshape(PP * T, [xp, yp, zp, nf]);
end
