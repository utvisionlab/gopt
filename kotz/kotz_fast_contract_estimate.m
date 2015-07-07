function [C, res] = kotz_fast_contract_estimate(Data, initC, a, c, options)
%
%   C=contract_estimate(Data,initC,func)
%
%  This function is for estimating the ML parameters for non-concave case
%     (a<n/2)
%  of Kotz-type distribution (KTD) which has the following density:
%
%  f(X)=Gamma(0.5*n)/(pi^(0.5*n))*c/Gamma(a)*1/b^a*det(C)^-0.5*u^(a*c-n/2)*exp(-U^c/b)
%
%  where  U=X' C^{-1} X,  n  is the dimension.
%  The parameters of the distribution are the nxn Covariance matrix  C 
%  and the shape parameters  a, c .  Ns  is the number of samples
%  Since C is true covariance b is equal to 
%      b=(n*gamma(a)/gamma(a+1/c))^(c)
%  The mean value of the generalized gamma is equal to
%     gamma(a+1/c)/gamma(a)*b^(1/c)
%
%  The methodolpgy is explained in siam paper
%
%  Inputs:
%    Data: Matrix with Ns columns each of which representing an
%          n-dimensional data point
%    initC : Te initial covariance matrix
%    a & c : shape parameters of Kotz 
%    options :  Structure with following Fields
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
b = (d*gamma(a)/gamma(a+1/c))^c;
func = @(x)2*1/b*(x.^c)+2*(d/2-a.*c).*log(x);
gfunc = @(x)2*c/b*(x.^(c-1))+2*(d/2-a.*c)./x;
n = size(Data,2);
C = initC;
tocs = zeros(1, options.MaxIter);
fvals = zeros(1, options.MaxIter);
tic;
for counter = 1:options.MaxIter
    Cchol = chol(C,'lower');
    midval = sum((Cchol\Data).^2,1);
    val = (a*b/mean(midval.^c))^(1/c);
    fv = func(midval);
    f= 2*n* sum(log(diag(Cchol))) +  sum(fv);
    cons = gfunc(midval*val);     %C= 1/val*C;
    Datac = bsxfun(@times,cons,Data);
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
