function  f = lpdfnorm(x,m,s)
%LPDFNORM The log of the normal density function
%
%         f = lpdfnorm(x,Mean,StandardDeviation)

if nargin<3, s=1; end
if nargin<2, m=0; end
f = -log(s)-log(2*pi)/2-((x-m)./s).^2/2;

