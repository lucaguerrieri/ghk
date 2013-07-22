function ldens = lpdfbeta(x,a,b);
% log BETA PDF
  ldens = gammaln(a+b) - gammaln(a) - gammaln(b) + (a-1)*log(x) + (b-1)*log(1-x);

% 10/11/03 MJ adapted from a GAUSS version by F. Schorfheide, translated
%             to Matlab by R. Wouters.  
%             use Matlab gammaln instead of lngam