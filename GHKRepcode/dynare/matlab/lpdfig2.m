function ldens = lpdfig2(x,s,nu)
% log INVERSE GAMMA (type 2)
%
% X ~ IG2(s,nu)
% X = inv(Z) where Z ~ G(nu/2,2/s) (Gamma distribution) 
% 
% See L. Bauwens, M. Lubrano and J-F. Richard [1999, appendix A] for more
% details.
%
% stephane.adjemian@cepremap.cnrs.fr [01/16/2004]

ldens = - gammaln(nu/2) - (nu/2)*log(2/s) - .5*(nu+2)*log(x) -.5*s./x;