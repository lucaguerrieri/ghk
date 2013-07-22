function [x,status]=bicgstab(func,b,x,tole,kmax,varargin)
  status = 0;
  r=b-feval(func,x,varargin{:});
  rh_0 = r;
  rh = r;
  rho_0 = 1;
  alpha = 1;
  w = 1;
  v = 0;
  p = 0;
  k = 0;
  rho_1 = rh_0'*r;
  tolr = tole*norm(b);
  
  while norm(r) > tolr & k < kmax
    k = k+1;
    beta = (rho_1/rho_0)*(alpha/w);
    p = r+beta*(p-w*v);
    v = feval(func,p,varargin{:});
    alpha = rho_1/(rh_0'*v);
    r = r-alpha*v;
    t = feval(func,r,varargin{:});
    w = (t'*r)/(t'*t);
    rho_0 = rho_1;
    rho_1 = -w*(rh_0'*t);
    x = x+alpha*p+w*r;
    r = r-w*t;
  end
if k == kmax
  status = 1;
  warning(sprintf('BICSTABN didn''t converge after %d iterations: norm(r) = %g',kmax,norm(r)));
end