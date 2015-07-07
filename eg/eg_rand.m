function D = eg_rand(C, a, Ns)
%
%  D = eg_rand(C, a, Ns)
%
%  EG_RAND generates samples from an elliptical gamma distribution (EGD)
%  which has the following density:
%
%  f(X)=Gamma(0.5*n)/(pi^(0.5*n))*1/Gamma(a)*1/b^a*det(C)^-0.5*u^(a-n/2)*exp(-U/b)
%
%  where  U=X' C^{-1} X,  b= n/a  and  n  is the dimension.
%  The parameters of the distribution are the nxn Covariance matrix  C 
%  and the shape parameter  a .  Ns  is the number of samples 
%
%  The function outputs a matrix with Ns columns each of which
%  representing an n-dimensional data point.
%
%  Original author: Reshad Hosseini, Nov, 02, 2014
%
d = size(C,1);
b = d/a;
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
U = randn(d, Ns);
U = bsxfun(@rdivide, U, sqrt(sum(U.^2))); % Unit data is generated
Ru = chol(C); 
D = Ru' * bsxfun(@times, U, sqrt(R)); 