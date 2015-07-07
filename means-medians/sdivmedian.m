function out=sdivmedian(A,maxit,verbose)
% SDIVMEDIAN -- compute s-median of psd matrices
%
%    OUT = SDIVMEDIAN(A) 
%    OUT = SDIVMEDIAN(A,MAXIT)
%
%  A: cell array of psd matrices
%  maxit: optional param to control # iterations, hence accuracy
%
%  out: structure containing the median and info about algo run
%
% 
%
% Author: Suvrit Sra <suvrit@gmail.com>
% (c) 2013 Suvrit Sra
% 
% See also: sdivmean
%

  if (~exist('maxit','var'))
    maxit=40;
  end
  if (~exist('verbose', 'var'))
    verbose=0;
  end
   
  m=numel(A);
  n=size(A{1});
  x = arithmeticMean(A);

  out.f = nan*ones(maxit,1);
  out.ft = nan*ones(maxit,1);
  out.ng = nan*ones(maxit,1);
   
  out.start=tic;
   
  for i=1:maxit
    [x,f,g]=smedianGX(x,A);
    out.ft(i)=toc(out.start);
    out.f(i)=f;
    ng=norm(g(:));
    out.ng(i)=ng;

    if (verbose > 0)
      fprintf('%d\t%d\t|g| = %d\n',i,f,ng);
    end
    
    if (ng < 5e-9)
      fprintf('|grad| < 1e-7....DONE\n');
      break;
    end
  end
  out.median=x;
end

function [x,f,g] = smedianGX(y,A)
  x = zeros(size(y));
  w = 1/numel(A);
  f = 0;

  for i=1:numel(A)
    d(i) = smet(y,A{i});
    x = x + 2*w*inv(y+A{i})/d(i);
  end
  f = w*sum(d);
  al = w*sum(1./d);
  g = x/2 - 0.5*al*inv(y);
  x = al*inv(x);
end

function d = smet(A,B)
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
  d=sqrt(d);
end

function am = arithmeticMean(A)
   am = A{1};
   for i=2:numel(A)
      am = am + A{i};
   end
   am = am/numel(A);
end