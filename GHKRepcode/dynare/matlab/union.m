function x = union(a,b)

xx = sort([a;b]) ;
x(1) = xx(1) ;
for i = 2:size(xx,1)
	if xx(i) ~= xx(i-1)
		x =[x;xx(i)];
	end
end
return ;


