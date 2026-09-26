function model = batchModelFromSingle(fun)
%% Description:
%   Adapter: turn a model written for ONE voxel at a time into the batched
%   interface of fitBatchedLM.m by looping over the voxels of each batch.
%   Use it for models that cannot be vectorised (e.g. anything built on
%   integral(), like vUS_3D_num_wrapper.m):
%       fun   = @(x, cols) vUS_3D_num_wrapper(x, tau(cols), k0, sigma);
%       model = batchModelFromSingle(fun);
%       x = fitBatchedLM(model, x0, g1, lb, ub, opts);
%   NOTHING IS VECTORISED HERE. The number of model calls is the same as in a serial
%   lsqnonlin loop, so with a slow per-voxel model the fit is only as fast as that model
%   (the solver overhead is still shared across voxels). For a real speed-up write the
%   model over a whole batch, as vUS_3D_quad_batchModel.m does.
%
% Inputs:
%   fun: function handle  y = fun(x, cols)  or  y = fun(x, cols, id)  for one voxel
%       x: [1, nP] parameters in physical units;  cols: [L, 1] data sample indices;  id: voxel index
%       y: L-element vector (complex or real). A second output dy [L, nP] (analytic
%       Jacobian) is used if fitBatchedLM is run with opts.jacobian = 'analytic'
%
% Outputs:
%   model: function handle  [Y, J] = model(X, cols, ids)  in the fitBatchedLM.m convention

    model = @(X, cols, ids) loopModel(fun, X, cols, ids);
end

function [Y, J] = loopModel(fun, X, cols, ids)
    [nP, Bc] = size(X);  L = numel(cols);
    withId = nargin(fun) >= 3;
    Y = zeros(L, Bc, 'like', X);
    if nargout > 1, J = zeros(L, nP, Bc, 'like', X); end
    for b = 1:Bc
        args = {X(:, b).', cols};
        if withId, args{3} = ids(b); end
        if nargout > 1
            [y, dy] = fun(args{:});
            J(:, :, b) = dy;
        else
            y = fun(args{:});
        end
        Y(:, b) = y(:);
    end
end
