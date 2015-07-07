function [C, res] = ecd_contract_estimate(Data, initC, func, gfunc, options)
%
%   [C, res] = ecd_contract_estimate(Data, initC, func, gfunc)
%
%   This function calculates the covariance in log-contractive case
%
%   C - \sum_{i=1}^{n}  x_i @func(x_i' C^-1 x_i) x_i' =0
%
%   Where C is the covariance matrix
%
%   gfunc is the derivative of func ans is a log non-expensive map 
%    for example: func=@(x)c1*x.^a+c2./x+c3 s.t. c1,c2,c3>0 & -1<a<1
%
%  Inputs:
%    Data    : Matrix with Ns columns each of which representing an
%              n-dimensional data point
%    initC   : Te initial covariance matrix
%    (g)func : see above description
%    options : Structure with following Fields
%          TolX: The tolerence for covariance
%          TolY: Tolerance for the function 
%          MaxIter: The maximum number of iterations
%          Display: 1 if we want the results to be shown
%
%  The function outputs the Covariance C and result of iteration
%
%  Original author: Reshad Hosseini, Nov, 03, 2014

if nargin < 4
    disp('Not Enough Input Arguments');
    return;
end
if nargin < 5
    options=struct;
end
if ~isfield(options,'TolX')
    options.TolX = 1e-10;
end
if ~isfield(options,'TolY')
    options.TolY = 1e-10;
end
if ~isfield(options,'MaxIter')
    options.MaxIter = 300;
end
if ~isfield(options,'Display')
    options.Display='OFF';
end
if strcmp(options.Display,'ON')
    fprintf('%10s%15s%15s\n','Iteration','Diff','fval');
end
d = size(Data,1);
n = size(Data,2);
C = initC;
tocs = zeros(1, options.MaxIter);
fvals = zeros(1, options.MaxIter);
tic;
for counter = 1:options.MaxIter
    Cchol = chol(C,'lower');
    midval = sum((Cchol\Data).^2,1);
    cons = gfunc(midval);
    Datac = bsxfun(@times,cons,Data);
    fv = func(midval);
    f = 2*n* sum(log(diag(Cchol))) +  sum(fv);
    fvals(counter) = f;
    tocs(counter) = toc;
    C_new = 1/n*Datac* Data';
    diffC = max(max(abs(C_new-C)));
    if strcmp(options.Display,'ON')
        fprintf('%10d %15.5e %15.5e\n',counter,diffC,f);
    end
    if counter>=2
        diffF = -fvals(counter) + fvals(counter-1);
        if diffF < options.TolY
            disp('Log-likelihood difference smaller than TolY');
            break;
        end
    end
    if diffC < options.TolX
        disp('Variance difference smaller than TolX');
        break;
    end
    C = C_new;
end
if nargout > 1
    res.fvals = fvals(1:counter);
    res.tocs = tocs(1:counter);
end
