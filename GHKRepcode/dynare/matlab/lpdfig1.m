function ldens = lpdfig1(x,s,nu)
% log INVERSE GAMMA (type 1)
%
% X ~ IG1(s,nu)
% X = sqrt(Y) where Y ~ IG2(s,nu) and Y = inv(Z) with Z ~ G(nu/2,2/s) (Gamma distribution) 
% 
% See L. Bauwens, M. Lubrano and J-F. Richard [1999, appendix A] for more
% details.
%
% stephane.adjemian@cepremap.cnrs.fr [01/16/2004]

ldens = log(2) - gammaln(nu/2) - (nu/2).*log(2/s) - (nu+1)*log(x) - .5*s./(x.^2);