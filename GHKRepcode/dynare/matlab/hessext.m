function [H] = hessext(f,x,options,varargin)
%---------------------------------------------------------------------------
%HESSEXT   Numerical approximation for hessian.
%          The method is Richardson`s extrapolation.
% Sample call
%   [H] = hessext(f,x,options,varargin)
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
  if nargin>3;
    if ~isempty(options);
      delta=options(1);
      toler=options(2);
    else
      delta=1e-12;
      toler=1e-9;
    end
  else
    delta=1e-12;
    toler=1e-9;
  end;
  ff=feval(f,x,varargin{:});
  nx=size(x,1);
  nf=size(ff,1);
  H=sparse(nx*nx,nf);
  for fi=1:nf;
    for xi=1:nx;
      for xj=xi:nx;
	err = 1;
	relerr = 1;
	h = 0.997;
	dx = zeros(nx,1);
	dy = zeros(nx,1);
	dx(xi)=h;      
	dy(xj)=h;      
	j = 1;
	fss=feval(f,x+dx+dy,varargin{:});
	fsm=feval(f,x+dx-dy,varargin{:});
	fms=feval(f,x-dx+dy,varargin{:});
	fmm=feval(f,x-dx-dy,varargin{:});
	D(1,1) = (fss(fi) -fsm(fi)-fms(fi)+ fmm(fi))/(4*h*h);
	while relerr>toler & err>delta
	  h = h/2;
	  dx = zeros(nx,1);
	  dy = zeros(nx,1);
	  dx(xi)=h;      
	  dy(xj)=h;      
	  fss=feval(f,x+dx+dy,varargin{:});
	  fsm=feval(f,x+dx-dy,varargin{:});
	  fms=feval(f,x-dx+dy,varargin{:});
	  fmm=feval(f,x-dx-dy,varargin{:});
	  D(j+1,1) = (fss(fi) -fsm(fi)-fms(fi)+ fmm(fi))/(4*h*h);
	  for k = 1:j,
	    D(j+1,k+1) = D(j+1,k) + (D(j+1,k)-D(j,k))/(4^k -1);
	  end
	  err = abs(D(j+1,j+1)-D(j,j));
	  relerr = 2*err/(abs(D(j+1,j+1))+abs(D(j,j))+eps);
	  j = j+1;
	  if j == 20
	    error(sprintf('Hessian evalutation didn''t converge. Equation %d, variables %d %d, relerr %e',fi,xi,xj,relerr))
	  end
	end
	n=size(D,1);
	H((xi-1)*nx+xj,fi)=D(n,n);
	H((xj-1)*nx+xi,fi)=D(n,n);
	D=0;
      end
    end
  end

% 10/12/2001 MJ modified initial h





