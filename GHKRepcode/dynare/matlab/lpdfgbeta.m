function ldens = lpdfgbeta(x,a,b,aa,bb);
% log (generalized) BETA PDF 
ldens = -betaln(a,b) + (a-1)*log(x-aa) + (b-1)*log(bb-x) - (a+b-1)*log(bb-aa);
%gammaln(a+b) - gammaln(a) - gammaln(b)
%betaln(a,b)
%pause
% 02/16/04 SA Interval [aa,bb] is the support of the probability density function. 
