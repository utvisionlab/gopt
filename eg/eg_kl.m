function kl=eg_kl(a1,C1,a2,C2)
%%
%  kl=EGDkl(a1,C1,a2,C2)
%
%  This function is for calculating the KL-divergence between
%  two elliptical gamma distributions (EGDs)
%
%  The form of this distribution with parameters a and C is as below
%
%  f(X)=Gamma(0.5*n)/(pi^(0.5*n))*1/Gamma(a)*1/b^a*det(C)^-0.5*u^(a-n/2)*exp(-U/b)
%
%  where  U=X' C^{-1} X,  b= n/a  and  n  is the dimension.
%  The parameters of the distribution are the nxn Covariance matrix  C
%  and the shape parameter  a .
%
%  a1 and C1 are the parameters of the ESTIMATED distribtion
%  and
%  a2 and C2 are the parameters of the MODEL distribtion
%
%  The function outputs kl which is the resulting KL-divergence
%
%  Original author: Reshad Hosseini, Nov, 02, 2014
%


n = size(C1,1);
b1 = n / a1;
b2 = n / a2;

d = eig(C1, C2, 'chol'); % faster than schur and svd

cons = paramA(d);

kl = -0.5 .* sum(log(d)) + gammaln(a2) + a2.*log(b2) - gammaln(a1)-...
    a2.*log(b1) + (a1-a2).*psi(a1) - a1 + a1.*b1./(n.*b2).* sum(d)...
    - (a2 - n/2) * (cons);
end

