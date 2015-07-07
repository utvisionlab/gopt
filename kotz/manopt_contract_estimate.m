function [C,info] = manopt_contract_estimate(Data,initC,func,gfunc,method,options)
%
%   [C,info] = manopt_estimate(Data,initC,a,method,options)
%
%   This function calculates the covariance in log-contractive case
%
%   C - \sum_{i=1}^{n}  x_i @gfunc(x_i' C^-1 x_i) x_i' =0
%
%   Where C is the covariance matrix
%
%   gfunc is the derivative of func ans is a log non-expensive map 
%    for example: gfunc=@(x)c1*x.^a+c2./x+c3 s.t. c1,c2,c3>0 & -1<a<1
%
%
%  Inputs:
%    Data: Matrix with Ns columns each of which representing an
%          n-dimensional data point
%    initC : Te initial covariance matrix
%    (g)func : see above description
%    options :  Structure with following Fields
%          tolgradnorm: Tolerance for the norm of gradient 
%          maxiter: The maximum number of iterations
%
%  The function outputs the Covariance C and result of iteration
%
%  Original author: Reshad Hosseini, Nov, 02, 2014

if nargin<5
    method='LBFGS';
end
if nargin<4
    disp('Not Enough Input Arguments');
    return;
end
if nargin<6
    options=struct;
end
if ~isfield(options,'tolgradnorm')
    options.tolgradnorm=1e-7;
end
if ~isfield(options,'maxiter')
    options.maxiter=300;
end
if ~isfield(options,'verbosity')
    options.verbosity=0;
end

m = size(Data,1);

% Create the problem structure.
manifold = spdfactory(m);
problem.M = manifold;

% Define the problem cost function and its gradient.
problem.costgrad = @(x) fun(x,Data,func,gfunc,manifold);

% Numerically check gradient consistency.
if false
    checkgradient(problem);
end


warning('off', 'manopt:getHessian:approx');
disp(['Method = ' method]);
switch method
    case 'SD', 
        [xout xcost info] = steepestdescent(problem,initC,options);
    case 'LBFGS',
        [xout xcost info] = lbfgs(problem,initC,options);
    case 'TR',
        [xout xcost info] = trustregions(problem,initC,options);  
    case 'CG',
        [xout xcost info] = conjugategradient(problem,initC,options);
end
C=inv(xout);


function [f,gf]=fun(x,Data,func,gfunc,manifold)
[m n] = size(Data);
[Cchol, p] = chol(x,'upper');
if p > 0
    f = Inf;
    gf = Inf(size(x));
else
    fv = func(sum((Cchol*Data).^2,1));    
    f = -2*n* sum(log(diag(Cchol))) +  sum(fv);
    if nargout>1
        cons = gfunc(sum((Cchol*Data).^2,1));
        Datac = bsxfun(@times,cons,Data);
        gfC = -n*inv_posdef(x) + Datac* Data';
        gf = manifold.egrad2rgrad(x, gfC);
    end
end

