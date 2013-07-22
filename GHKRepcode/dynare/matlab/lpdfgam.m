function  ldens = lpdfgam(x,a,b);
% log GAMMA PDF 
  ldens = -gammaln(a) -a*log(b)+ (a-1)*log(x) -x/b ;

% 10/11/03  MJ adapted from an earlier GAUSS version by F. Schorfeide,
%              translated to MATLAB by R. Wouters  
%              use MATLAB gammaln rather than lngam
