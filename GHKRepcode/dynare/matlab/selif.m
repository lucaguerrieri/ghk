% Copyright (C) 2001 Michel Juillard
%
function x = selif(a,b)

if size(b,2) ~= 1
	error ('The second argument in SELIF must be à column vector') ;
end

x = a(find(b == 1),:) ;

return ;

