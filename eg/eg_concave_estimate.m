function [C,res] = eg_concave_estimate(Data, initC, a)
%
%   [C,res] = eg_concave_estimate(Data, initC, a)
%
%  This function is for estimating the ML parameters for concave case
%     (a>=n/2)
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
%          TolX: The tolerence for covariance
%          TolY: Tolerance for the function 
%          MaxIter: The maximum number of iterations
%          Display: 1 if we want the results to be shown
%
%  The function outputs the Covariance C and result of iteration
%
%  Original author: Reshad Hosseini, Nov, 02, 2014


if nargin<2
    disp('Not Enough Input Arguments');
    return;
end
if nargin < 4
    options=struct;
end
if ~isfield(options,'TolX')
    options.TolX = 1e-10;
end
if ~isfield(options,'TolY')
    options.TolY = 1e-10;
end
if ~isfield(options,'MaxIter')
    options.MaxIter = 200;
end
if ~isfield(options,'Display')
    options.Display='ON';
end
if strcmp(options.Display,'ON')
    fprintf('%10s%15s%15s\n','Iteration','Diff','fval');
end
d = size(Data,1);
if a < d/2
    disp('Please Run the contractive case');
    return;  
end
b = d/a;
n = size(Data,2);
func = @(x)2/b*(x)+2*(d/2-a).*log(x);
%counter = 0;
tocs = zeros(1,options.MaxIter);
fvals = zeros(1,options.MaxIter);
tic;
% first the data are whitened
C = Data * Data' / n;
R = chol(C); % C=R' R ( R upper trangular)
Rinv = R \ eye(size(R,1)); % Faster version of inv_triu(R); 
% Transformation to whiten the data
Data = Rinv' * Data; % whitenned data

M = Rinv * inv_posdef(initC)  * Rinv';
logdetR = sum(log(diag(R)));

for counter = 1:options.MaxIter
    [M_new, midval] = updateConvex(Data, a, b, M);
    fv = func(midval);
    Mchol = chol(M);
    f = 2*n* (-sum(log(diag(Mchol)))+logdetR) +  sum(fv);
    fvals(counter) = f;
    tocs(counter) = toc;
    diffC = max(max(abs(M_new-M)));
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
    M = M_new;
end
if nargout > 1
    res.fvals = fvals(1:counter);
    res.tocs = tocs(1:counter);
end
C = R' * inv_posdef(M) * R;

function [Mout,u]=updateConvex(dataX, a, b, M)
q = size(dataX,1);
n = size(dataX,2);
sM = sqrtm(M);
u =  sum((sM * dataX).^2, 1);
Mout = bsxfun(@times, dataX, 1./u) * dataX';
Mout = b*(a-q/2)/n * sM * Mout * sM' + b/2 * eye(q);