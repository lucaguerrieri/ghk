function d1 = bksupk(ny,fid,jcf,icc1)

global ykmax_ iter_

icf = [1:jcf-1] ;
ir = [(iter_-1)*ny+1:ny*iter_] ;
irf = icc1+(iter_-1)*ny ;
d1 = zeros(iter_*ny,1) ;

ofs = (((iter_-1)*ny+1)-1)*jcf*8 ;
junk = fseek(fid,ofs,-1) ;
c = fread(fid,[jcf,ny],'float64')' ;

d1(ir) = c(:,jcf) ;
ir = ir-ny ;

i = 2 ;

while i <= ykmax_ | i <= iter_
	irf1 = selif(irf,irf<=iter_*ny) ;

	ofs = (((iter_-i)*ny+1)-1)*jcf*8 ;
	junk = fseek(fid,ofs,-1) ;
	c = fread(fid,[jcf,ny],'float64')' ;

	d1(ir) = c(:,jcf) - c(:,1:size(irf1,1))*d1(irf1) ;
	ir = ir - ny ;
	irf = irf - ny ;
	i = i + 1 ;
end

while i <= iter_

	ofs = (((iter_-i)*ny+1)-1)*jcf*8 ;
	junk = fseek(fid,ofs,-1) ;
	c = fread(fid,[jcf,ny],'float64')' ;

	d1(ir) = c(:,jcf)-c(:,icf)*d1(irf) ;
	ir = ir-ny ;			
	irf = irf-ny ;
	i = i+1;
end

