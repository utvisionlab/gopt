function Data = kotz_rand(a,c,C,Ns)
%
%  Data=kotz_rand(a,c,C,Ns)
%
%  Kotzrnd generates samples from a Kotz-type distribution (KTD)
%  which has the following density:
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
%  The function outputs a matrix with Ns columns each of which
%  representing an n-dimensional data point.
%

n = size(C,1);
U = randn(n,Ns);
U = bsxfun(@rdivide,U,sqrt(sum(U.^2))); % Unit data is generated
b = (n*gamma(a)/gamma(a+1/c))^c;
if exist('randg','file') == 3
    % First check for matlab randg function
    samplegamma = @(x,y)randg(x, [1 y]);
elseif exist('randgamma','file') == 3
    % Then Check for installed lightspeed randgamma
    samplegamma = @(x,y)randgamma(x * ones(1,y));
else
    % At last run randraw which is twice slower
    samplegamma = @(x,y)randraw('gamma', x, [1 y]);
end
R = samplegamma(a, Ns);
R = R*b;
R = R.^(1/c);
Ru = chol(C); 
Data = Ru'*bsxfun(@times,U,R);
