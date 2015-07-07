function out = sdivmean(A,maxit,verbose)
% SDIVMEAN -- compute s-mean of psd matrices
%
%    OUT = SDIVMEAN(A) 
%    OUT = SDIVMEAN(A,MAXIT)
%
%  A: cell array of psd matrices
%  maxit: optional param to control # iterations, hence accuracy
%
%  out: structure containing mean and info about the algorithm's run
%
%  Author: Suvrit Sra <suvrit@gmail.com>
%  (c) 2011 Suvrit Sra
% 
% See also: sdivmedian
%
  
  if (~exist('maxit','var'))
    maxit=40;
  end
  if (~exist('verbose', 'var'))
    verbose=0;
  end

 out.f = nan*ones(maxit,1);
 out.ft = nan*ones(maxit,1);
 out.ng = nan*ones(maxit,1);
   
 out.start=tic;
   
  % initialize
  x = 0.9*arithmeticMean(A);
  xi = inv(x);
  
  for i=1:maxit
    [x,f,g] = smeanGX(xi, A);
    xi = inv(x);

    out.ft(i)=toc(out.start);
    out.f(i)=f;
    ng=norm(g(:));
    out.ng(i)=ng;
    
    if (verbose > 0)
      fprintf('%d\t%d\t|g| = %d\n',i,f,ng);
    end
  end
  out.mean=x;
end

function s=sumCell(A)
  s = A{1};
  for i=2:numel(A)
    s = s + A{i};
  end
end

function [x,f,g] = smeanGX(xi, A)
  x = zeros(size(A{1}));
  f = 0;
  g = zeros(size(x));
  for i=1:numel(A)
    x = x + inv(xi+A{i});
    f = f + sdiv(x, A{i});
  end
  x = (2/numel(A))*x;
  g = x - inv(xi);
end

function [f, g] = smeanobj(x, A)
  f = 0;
  g = zeros(size(x));
  for i=1:numel(A)
    f = f + sdiv(x, A{i});
    if (nargout > 1)
      g = g + inv(x+A{i});
    end
  end

  if (nargout > 1)
    g = g - numel(A)*inv(x)/2;
  end
end

function d = sdiv(A,B)
  n=size(A,1);

  if (n <= 20)
    d = log(det(A+B)) - 0.5*log(det(A)) - 0.5*log(det(B)) - n*log(2);
  else
    r = chol(A+B);
    r1=chol(A);
    r2=chol(B);
    nr1=diag(r) .* diag(r);
    dr1=diag(r1) .* diag(r2);
    t = sum(log(nr1 ./dr1));
    d = t - n*log(2);
  end
end

function am = arithmeticMean(A)
   am = A{1};
   for i=2:numel(A)
      am = am + A{i};
   end
   am = am/numel(A);
end