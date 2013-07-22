function b=f_var(x,a,nx)
  x=reshape(x,nx,nx);
  b=x-a*x*a';
  b=b(:);