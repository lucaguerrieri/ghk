% Copyright (C) 2001 Michel Juillard
%
% symmetric formula to compute the Jacobian
%
function jacobia_=jacob_a(func,z)
		
global d1_ gstep_

func = str2func(func);
z = reshape(z,max(size(z)),1) ;
d1_ = feval(func,z) ;
nz = size(z,1) ;
jacobia_ = zeros(size(d1_,1),nz) ;
dh = max(abs(z), gstep_*ones(nz,1))*eps^(1/3) ;
xdh1 = z ;
xdh2 = z ;
for j = 1:nz
	xdh1(j) = z(j)-dh(j) ;
	xdh2(j) = z(j)+dh(j) ;
	h = xdh2(j)-xdh1(j) ;
	jacobia_(:,j) = ((feval(func,xdh2)-feval(func,xdh1))/h) ;
	xdh1(j) = z(j) ;
	xdh2(j) = z(j) ;
end

% 10/03/02 MJ creation
