function fcn=gaussFcn(paramCell)
%Generate a handle to the fitted Gaussian.
%
%   Fcn=gaussFcn({D,A,my,Sigma})  %input is the parameter cell output from gaussfitn
%   Fcn=gaussFcn( D,A,my,Sigma )  %parameters given as separate arguments
%
%Fcn is vectorized with calling syntax Fcn(X,Y,...,Z)

    
    
    [D,A,mu,Sigma]=deal(paramCell{:});
    
    fcn=@(varargin) Gfun(D,A,mu,Sigma,varargin{:});

end

function fvals=Gfun(D,A,mu,Sigma,varargin)

 N=numel(mu);

 assert(N==numel(varargin),'Too many inputs');

 sz=size(varargin{1});


 dX=cell2mat( cellfun(@(r)r(:), varargin,'uni',0) )-mu(:).';

 fvals=reshape(  A*exp( sum(  (dX/(-2*Sigma)).*dX ,2))  ,sz);

 if D, fvals=fvals+D; end

end