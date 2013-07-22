function ldens = lpdfig(x,a,b)
% log INVERSE GAMMA
   ldens = log(2) - gammaln(b/2) + (b/2).*log(a/2) - ( (b+1)/2 ).*log(x.^2) - a./(2*x.^2);

% 10/11/03 MJ adapted from an earlier GAUSS version by F. Schorfeide,
%             then translated to Matlab by R. Wouters
%             uses Matlab gammaln rather than lngam
% 12/01/03 MJ modified according to MR suggestion