function fs_dyn(a,b)
  nu=b;
  s=b*a^2;
  mu=sqrt(s/2)*gamma((nu-1)/2)/gamma(nu/2);
  sigma=sqrt(s/(nu-2)-mu^2);

  disp('     FS to DYNARE conversion for inverse gamma parameters')
  disp(sprintf('          Frank''s parameters: %f %f',a,b))
  disp(sprintf('          mu = %f       sigma = %f',mu,sigma))