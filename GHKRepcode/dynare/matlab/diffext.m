function [J] = diffext(f,x,options,varargin)
%---------------------------------------------------------------------------
%DIFFEXT   Numerical approximation for hessian.
%          The method is Richardson`s extrapolation.
% Sample call
%   [D,err,relerr,n] = diffext('f',x,delta,toler)
% Inputs
%   f        name of the function
%   x        differentiation point
%   options  matrix of algorithm parameters
%   delta    error goal (1e-12)
%   toler    relative error goal (1e-12)
% Return
%   J        Jacobian
%
% NUMERICAL METHODS: MATLAB Programs, (c) John H. Mathews 1995
%
% Modified F. Collard, August 2001
%---------------------------------------------------------------------------
  global gstep_
  if nargin>3;
    if ~isempty(options);
      delta=options(1);
      toler=options(2);
    else
      delta=1e-12;
      toler=1e-12;
    end
  else
    delta=1e-12;
    toler=1e-12;
  end;
  ff=feval(f,x,varargin{:});
  nx=size(x,1);
  nf=size(ff,1);
  J=zeros(nf,nx);
  for fi=1:nf;
    for xi=1:nx;
      err = 1;
      relerr = 1;
      h=max(abs(x(xi)), gstep_)*eps^(1/3) ;
      dx = zeros(nx,1);
      dx(xi)=h;      
      j = 1;
      fs=feval(f,x+dx,varargin{:});
      fm=feval(f,x-dx,varargin{:});
      D(1,1) = (fs(fi) - fm(fi))/(2*h);
      while relerr>toler & err>delta & j<12
	h = h/2;
	dx(xi)=h;      
	fs=feval(f,x+dx,varargin{:});
	fm=feval(f,x-dx,varargin{:});
	D(j+1,1) = (fs(fi) - fm(fi))/(2*h);
	for k = 1:j,
	  D(j+1,k+1) = D(j+1,k) + (D(j+1,k)-D(j-1+1,k))/(4^k -1);
	end
	err = abs(D(j+1,j+1)-D(j,j));
	relerr = 2*err/(abs(D(j+1,j+1))+abs(D(j,j))+eps);
	j = j+1;
      end
      n=size(D,1);
      J(fi,xi)=D(n,n);
      clear D;
    end
  end

% 10/12/2001 MJ modified initial h







