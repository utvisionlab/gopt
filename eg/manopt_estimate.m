function [C,info] = manopt_estimate(Data,initC,a,method,options)
%
%   [C,info] = manopt_estimate(Data,initC,a,method,options)
%
%  This function is for estimating the ML parameters of
%  elliptical gamma distribution (EGD) which has the following density:
%
%  f(X)=Gamma(0.5*n)/(pi^(0.5*n))*1/Gamma(a)*1/b^a*det(C)^-0.5*u^(a-n/2)*exp(-U/b)
%
%  where  U=X' C^{-1} X,  b= n/a  and  n  is the dimension.
%  The parameters of the distribution are the nxn Covariance matrix  C
%  and the shape parameter  a .
%
%  Inputs:
%    Data: Matrix with Ns columns each of which representing an
%          n-dimensional data point
%    initC : Te initial covariance matrix
%    a     : shape parameter of EG 
%    options :  Structure with following Fields
%          tolgradnorm: Tolerance for the norm of gradient 
%          maxiter: The maximum number of iterations
%
%  The function outputs the Covariance C and result of iteration
%
%  Original author: Reshad Hosseini, Nov, 02, 2014

if nargin<4
    method='LBFGS';
end
if nargin<3
    disp('Not Enough Input Arguments');
    return;
end
if nargin<5
    options=struct;
end
if ~isfield(options,'tolgradnorm')
    options.tolgradnorm=1e-5;
end
if ~isfield(options,'maxiter')
    options.maxiter=200;
end

d = size(Data,1);
b = d/a;
func = @(x)2/b*(x)+2*(d/2-a).*log(x);
gfunc = @(x)2/b+2*(d/2-a)./x;

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
    % = chol(x,'upper');
    fv = func(sum((Cchol*Data).^2,1));    
    f = -2*n* sum(log(diag(Cchol))) +  sum(fv);
    if nargout>1
        cons = gfunc(sum((Cchol*Data).^2,1));
        Datac = bsxfun(@times,cons,Data);
        gfC = -n*inv_posdef(x) + Datac* Data';
        gf = manifold.egrad2rgrad(x, gfC);
    end
end

