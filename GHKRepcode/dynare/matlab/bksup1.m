% Copyright (C) 2001 Michel Juillard
%
function d = bksup1(ny,jcf,iyf,c)

global iter_

ir = [(iter_-2)*ny+1:ny+(iter_-2)*ny] ;
irf = iyf+(iter_-1)*ny ;
icf = [1:size(iyf,2)] ;
d = c(:,jcf) ;

for i = 2:iter_
	d(ir) = c(ir,jcf)-c(ir,icf)*d(irf) ;
	ir = ir-ny ;
	irf = irf-ny ;
end

