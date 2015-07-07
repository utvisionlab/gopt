function [C, res] = ecd_ls_contract_estimate(Data, initC, func, gfunc, options)
%
%   [C, res] = ecd_ls_contract_estimate(Data, initC, func, gfunc)
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
%  The method uses line-search approach
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
val = 0.1;
for counter = 1:options.MaxIter
    Cchol = chol(C,'lower');
    midval = sum((Cchol\Data).^2,1);
    logd = 2*n* sum(log(diag(Cchol)));
    const = 1;
    flag = 0;
    count = 1;
    emsave = Inf;
    while 1
        cons = gfunc(1/const * midval);
        Datac = bsxfun(@times,cons,Data);
        if count == 1
            fv = func(1/const * midval);
            f =  n*log(const) + logd + sum(fv); 
            fvals(counter) = f;
        end
        tocs(counter) = toc;
        C_new = 1/n*Datac* Data';
        OPTS.disp = 0;
        emax = eigs(C_new,C*const,1,'LM',OPTS);
        emin = eigs(C_new,C*const,1,'SM',OPTS);
        em = max(emax, 1/emin);
        if em < emsave
            C_n = C_new;
        end
        if (emax>1 && emin>1)
            if flag == 2
                val = val/2;
            end
            flag = 1;
            const = const*(val+1);
            count = count+1;
        elseif (emax<1 && emin<1)
            if flag == 1
                val = val/2;
            end
            flag = 2;
            const = const/(val+1);
            count = count+1;
        else
            break;
        end
        if count == 30
            break;
        end
    end
    val = val * 2;
    diffC = max(max(abs(C_n-C)));
    if strcmp(options.Display,'ON')
        fprintf('%10d %15.5e %15.5e\n',counter,diffC,f);
    end
    if counter >= 2
        diffF = -fvals(counter) + fvals(counter-1);
        if abs(diffF) < options.TolY
            disp('Log-likelihood difference smaller than TolY');
            break;
        end
    end
    if diffC < options.TolX
        disp('Variance difference smaller than TolX');
        break;
    end
    C = C_n;
end
if nargout > 1
    res.fvals = fvals(1:counter);
    res.tocs = tocs(1:counter);
end
