% Copyright (C) 2001 Michel Juillard
%
function x=fff(y)

global ykmin_ ykmax_ iy_

iyr = find(reshape(iy_'>0,(ykmin_+ykmax_+1)*size(iy_,2),1)) ;
y = kron(ones(ykmin_+ykmax_+1,1),y) ;
y = y(iyr,:) ;
x = ff_(y) ;

return ; 