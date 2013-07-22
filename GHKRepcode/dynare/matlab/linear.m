% Copyright (C) 2001 Michel Juillard
%
function linear(x)
%	LINEAR : LINEAR (x,'filename')

global iy_ ykmin_ ykmax_ xkmin_ xkmax_ jacobia_
global lgy_ lgx_ iter_ fname_

x = reshape(x,max(size(x)),1) ;
xx = x ;
nn = size(iy_,2) ;
nv = max(iy_(size(iy_,1),:)) ;
indic = zeros(0,0) ;

if size(x,1) == nv
	r = find(reshape(iy_'>0,(ykmin_+ykmax_+1)*nn,1)) ;
	n = abs(lgy_) ;
	n = kron(ones(ykmin_+ykmax_+1,1),n) ;
	n = setstr(n(r,:)) ;
	m = [-ykmin_:ykmax_]'*ones(1,nn) ;
	m = reshape(m',(ykmin_+ykmax_+1)*nn,1) ;
	m = m(r) ;
elseif size(x,1) == nn
	l = iy_>0 ;
	l = l.*kron(ones(size(iy_,1),1),[1:nn]) ;
	i = nonzeros(reshape(iy_,size(iy_,1)*nn,1)) ;
	j = nonzeros(reshape(l,size(iy_,1)*nn,1)) ;
	s = ones(1,max(iy_(size(iy_,1),:))) ;
	indic = sparse(i,j,s,nv,nn) ;
	x = indic*x ;
	n = lgy_ ;
	m = zeros(nn,1) ;
else
	error ('Wrong number of arguments in LINEAR.') ;
end

jacob([fname_ '_ff'],x) ;

if ~ isempty(indic)
	jacobia_=jacobia_*indic ;
	clear indic ;
end

fprintf (1,'Periods  :  ') ;
fprintf (1,'%4g \n',iter_) ;
fprintf (1,'Endogenous variables : ') ; fprintf (1,'\n') ;
fprintf (1,lgy_) ;
fprintf (1,'\n') ;
fprintf (1,'Exogenous variables : ') ; fprintf (1,'\n') ;
fprintf (1,lgx_) ;
fprintf (1,'\n') ;
fprintf (1,'Linearization around :') ;
fprintf (1,'\n') ;

for i = 1:size(n,1) ;
	fprintf (1,n(i,:)) ; fprintf (1,'(%1g)',m(i)) ;
	fprintf (1,' = %15.6f \n',xx(i)) ;
end

fprintf (1,'\n') ;

for i = 1:size(jacobia_,1)
	for j = 1:size(jacobia_,2)
		if jacobia_(i,j) ~= 0
			if jacobia_(i,j) == 1
				if j == 1
					fprintf (1,n(j,:)) ; fprintf (1,'(%1g)',m(j)) ;
				else
					fprintf (1, ' + ') ; fprintf (1,n(j,:)) ; fprintf (1,'(%1g)',m(j)) ;
				end
			elseif jacobia_(i,j) == -1
				fprintf (1,' - ') ; fprintf (1,n(j,:)) ; fprintf (1,'(%1g)',m(j)) ;
			elseif jacobia_(i,j) > 0
				if j>1, fprintf (1,' + ') ; end
				fprintf (1,'%15.6g',jacobia_(i,j)) ;
				fprintf (1,'*') ; fprintf (1,n(j,:)) ; fprintf (1,'(%1g)',m(j)) ;
			else
				fprintf (1,'%15.6g',jacobia_(i,j)) ;
				fprintf (1,'*') ; fprintf (1,n(j,:)) ; fprintf (1,'(%1g)',m(j)) ;
			
			end
		end
	end
	fprintf (1,'\n') ;
end

return ;
