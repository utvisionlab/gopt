function cons = paramA(d)
%%
% Calculating parameter A (given by the following expectation
%
%   A = E[ \log \sum d_i N_i^2] - E[\log \sum N_i^2]
%
%  where N_is are indepependent zero-mean & univ-variance Gaussians
%

% d should be a column vector
if size(d,2) > 1
    d = d.';
end

d = sort(d);

n = length(d);
% Compute the length of approximation using sum of gamma distributions
Lo = finde(1e-5, d);
if true
    % true is for accurate and fast numerical methods for evaluation
    if Lo < 5000
        cons = elsd(d,Lo);
    else
        % If length is large switch to numerical integration
        cc = max(d)*2;
        l = d/cc;
        is = 1/(2*(n+1)*max(l));
        ia = 1-2*is*l;
        ib = 2*l./ia;
        ia = ia.^2;
        isb = 1/is;
        isa = is^2;
        C = digamma(1);
        fun = @(x) inint(x, ia, ib, isa, isb, C, l);
        cons = quadgk(fun,0,Inf);
        cons = cons/pi + log(cc) - log(2) - digamma(n/2);
    end
else
    % false is for slow and inaccurate monte-carlo integration
    u1 = randn(n, 100000) .^ 2;
    cons = mean(log(sum(bsxfun(@times, u1, d))) - log(sum(u1)));
end

end

function y = inint(x,ia,ib,isa,isb,C,l)
%% Function for value inside integral
% x should have one row and others should have one column
it = 0.5 * sum(atan(ib*x),1) - atan(isb.*x);
% doing for numerical stability
im = prod(bsxfun(@plus, ia , (2*l*x).^2).^(-0.25),1);
y =  1./( sqrt(isa+x.^2)) .* im .* ...
    ((C - 0.5 * log(isa+x.^2)) .* cos (it) + atan(isb.*x) .* sin(it));
end

function Lo = finde(fm,l)
%% Function for calculating length given error fm
fm = log(fm);
q = length(l);
qh = q/2;
b = 2 * l(1) * l(q) / (l(1) + l(q));
dv = (b./l);
logc0 = 0.5 * sum( log(dv) );
e = max(abs(1-dv));
errp = @(x) gammaln(qh+x)-gammaln(qh+1)-gammaln(x+2)+...
    (x+1)*log(e)-1*log((1-(qh+x)/x*e))+ logc0;

x0 = ceil(q/2/(1-e));
f0 = errp(x0);
if f0 <= fm
    Lo = x0;
    return;
end
x1 = x0 * 2;
f1 = errp(x1);
while  f1 >= fm
    x0 = x1;
    f0 = f1;
    x1 = x1*2;
    f1 = errp(x1);
end

iter = 0;
while (x1-x0) > 0.5
    iter = iter + 1;
    if iter > 100
        break;
    end
    xn = ( (fm-f0)*(x1-x0)/(f1-f0) + x0);
    fn = errp(xn);
    if fn <= fm
        x1 = xn;
        f1 = fn;
    end
    if fn > fm
        x0 = xn;
        f0 = fn;
    end
    if abs(f1-fm) < 1e-6
        x0 = x1;
    end
    if abs(f0-fm) < 1e-6
        x1 = x0;
        f1 = f0;
    end
end
Lo = ceil(x1);
end