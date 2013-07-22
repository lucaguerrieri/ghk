% Copyright (C) 2001 Michel Juillard
%
function x = bseastr(s1,s2)

m = size(s1,1) ;
x = zeros(m,1) ;
s1=upper(deblank(s1));
s2=upper(deblank(s2));

for im = 1:m
	key = s1(im,:) ;
	h = size(s2,1) ;
	l = 1 ;
	while l <= h
		mid = round((h+l)/2) ;
		temp = s2(mid,:) ;
		if ~ strcmp(key,temp)
 			for i = 1:min(length(key),length(temp))
				if temp(i) > key(i)
					h = mid - 1 ;
					break 
				else
					l = mid + 1 ;
					break 
				end
			end
		else
			x(im) = mid ;
			break 
		end
	end
end

