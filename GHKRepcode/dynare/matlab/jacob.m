% Copyright (C) 2001 Michel Juillard
%
function jacob(func,z)
		
global d1_ jacobia_ gstep_

func = str2func(func);
z = reshape(z,max(size(z)),1) ;
d1_ = feval(func,z) ;
nz = size(z,1) ;
jacobia_ = zeros(size(d1_,1),nz) ;
dh = max(abs(z), gstep_*ones(nz,1))*eps^(1/3) ;
xdh = z ;
for j = 1:nz
	xdh(j) = xdh(j)+dh(j) ;
	h = xdh(j)-z(j) ;
	jacobia_(:,j) = ((feval(func,xdh)-d1_)/h) ;
	xdh(j) = z(j) ;
end

% 8/27/2000 modified with function name MJ
