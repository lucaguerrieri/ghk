% Copyright (C) 2001 Michel Juillard
%
function y=reshapel(x,m,n)
%RESHAPEL Change size.
%	RESHAPEL(X,M,N) returns the M-by-N matrix whose elements
%	are taken linewise from X.  An error results if X does
%	not have M*N elements.

[mm,nn] = size(x);

if mm*nn ~= m*n
    error('Matrix must have M*N elements.')
end

[i,j,s] = find(x') ;
if size(i,2) ~= 1,i = i';end
k = (j-1)*nn+i ;
j = rem(k-1,n)+1 ;
i = (k-j)/n+1 ;
y = full(sparse(i,j,s,m,n)) ;

return ;