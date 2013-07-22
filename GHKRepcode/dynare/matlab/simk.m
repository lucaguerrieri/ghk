% Copyright (C) 2001 Michel Juillard
%
function simk

global dynatol_ maxit_ slowc_ scalv_
global iy_ ykmin_ ykmax_ iyr0 ct_ jacobia_ d1_
global gstep_ ys_ y_ it_ iter_
global ncc timing_ broyden_ fname_

func_name = [fname_ '_ff'];
nk = ykmin_ + ykmax_ + 1 ;
ny = size(iy_,2) ;
icc1 = iy_(nk,:) > 0;

for i = 1:ykmax_ -1
  icc1 = [iy_(nk-i,:) | icc1(1,:); icc1] ;
end

icc1 = find(icc1') ;
iy = iy_ > 0 ;
isc = cumsum(sum(iy',1))' ;
iyr0 = find(iy_') ;
ncc1 = size(icc1,1) ;
ncc = ncc1 + 1 ;
ncs = size(iyr0,1) ;

ky = zeros(ny,nk) ;            % indices of variables at each lead or lag
lky = zeros(nk,1) ;
for i = 1:nk
  j = find(iy_(i,:))' ;
  if isempty(j)
    lky(i) = 0;
  else
    lky(i) = size(j,1) ;
    ky(1:lky(i),i) = j ;
  end
end

jwc = find(iy(2:ykmax_+1,:)') ; % indices of columns for
                                % triangularization
				% as many rows as lags in model

if isempty(jwc)
  jwc = 0 ;
  ljwc = 0 ;
  temp = icc1 ;
else
  ljwc = size(jwc,1) ;          % length of each row in jwc
  temp = union(jwc,icc1) ;      % prepares next iteration
end

j1 = ky(1:lky(1),1) ;
lj1 = lky(1) ;

for i = 2:ykmin_
  [j1,lj1] = ffill(j1,lj1,selif(temp+(i-1)*ny,temp <= ny)) ;
  if ykmax_ == 1
    if lky(i+ykmax_) > 0
      [jwc,ljwc] = ffill(jwc,ljwc, ky(1:lky(i+ykmax_),i+ykmax_)+(ykmax_-1)*ny) ;
      if ljwc(i) == 0
 	temp = icc1;
      else
 	temp = union(jwc(1:ljwc(i),i),icc1) ;
      end
    else
      [jwc,ljwc] = ffill(jwc,ljwc,[]) ;
      temp = icc1 ;
    end
  else
    temp = temp(lj1(i)+1:size(temp,1),:) - ny ;
    if lky(i+ykmax_) > 0
      [jwc,ljwc] = ffill(jwc,ljwc,[temp;ky(1:lky(i+ykmax_),i+ykmax_)+(ykmax_-1)*ny]);
    else
      [jwc,ljwc] = ffill(jwc,ljwc,temp) ;
    end
    temp = union(jwc(1:ljwc(i),i),icc1) ;
  end
end

[j1,lj1] = ffill(j1,lj1,selif(temp+ykmin_*ny, temp <= ny)) ;
ltemp = zeros(ykmin_,1) ;
jwc1 = zeros(ncc1,ykmin_) ;

for i = 1:ykmin_
  temp = union(jwc(1:ljwc(i),i),icc1) ;
  ltemp(i) = size(temp,1) ;
  if ljwc(i) > 0
    jwc(1:ljwc(i),i) = indnv(jwc(1:ljwc(i),i),temp) ;
  end
  jwc1(:,i) = indnv(icc1,temp) ;
end

h1 = clock ;

disp (['-----------------------------------------------------']) ;
disp ('MODEL SIMULATION') ;
fprintf ('\n') ;

for iter = 1:maxit_
  h2 = clock ;
  y_ = y_(:);
  err_f = 0;
  
  fid = fopen('dynare.swp','w+') ;

  it_ = 1+ykmin_ ;
  ic = [1:ny] ;
  iyr = iyr0 ;
  i = ykmin_+1 ;
  while (i>1) & (it_<=iter_+ykmin_)
    h3 = clock ;
    if broyden_ & iter > 1
      d1_ = -feval(fh,y_(iyr));
    else
      jacob(func_name,y_(iyr)) ;
      d1_ = -d1_ ;
    end
    err_f = max(err_f,max(abs(d1_)));
    if lky(i) ~= 0
      j1i = ky(1:lky(i),i) ;
      w0 = jacobia_(:,isc(i-1)+1:isc(i)) ;
    else
      w0 = [];
    end
    ttemp = iy(i+1:i+ykmax_,:)' ;
    jwci = find(ttemp) ;
    if ~ isempty(jwci)
      w = jacobia_(:,isc(i)+1:isc(i+ykmax_)) ;
    end
    j = i ;
    while j <= ykmin_
      if ~isempty(w0)

	ofs = ((it_-ykmin_-ykmin_+j-2)*ny)*ncc*8 ;
	junk = fseek(fid,ofs,-1) ;
	c = fread(fid,[ncc,ny],'float64') ;
	c = c' ;

	if isempty(jwci)
	  w = -w0*c(j1i,1:ncc1) ;
	  jwci = icc1 ;
	else
	  iz = union(jwci,icc1) ;
	  ix = indnv(jwci,iz) ;
	  iy__ = indnv(icc1,iz) ;
	  temp = zeros(size(w,1),size(iz,1)) ;
	  temp(:,ix) = w ;
	  temp(:,iy__) = temp(:,iy__)-w0*c(j1i,1:ncc1) ;
	  w = temp ;
	  jwci = iz ;
	  clear temp iz ix iy__ ;
	end
	d1_ = d1_-w0*c(j1i,ncc) ;
	clear c ;
      end
      j = j + 1 ;
      if isempty(jwci)
	j1i = [];
	if lky(j+ykmax_) ~= 0
	  jwci = ky(1:lky(j+ykmax_),j+ykmax_) + (ykmax_-1)*ny ;
	  w = jacobia_(:,isc(j+ykmax_-1)+1:isc(j+ykmax_)) ;
	else
	  jwci = [] ;
	end
      else
	j1i = selif(jwci,jwci<(ny+1)) ;
	w0 = w(:,1:size(j1i,1)) ;
	if size(jwci,1) == size(j1i,1)
	  if lky(j+ykmax_) ~= 0
	    jwci = ky(1:lky(j+ykmax_),j+ykmax_)+(ykmax_-1)*ny ;
	    w = jacobia_(:,isc(j+ykmax_-1)+1:isc(j+ykmax_)) ;
	  else
	    jwci = [] ;
	  end
	else
	  jwci = jwci(size(j1i,1)+1:size(jwci,1),:)-ny ;
	  w = w(:,size(j1i,1)+1:size(w,2)) ; 
	  if lky(j+ykmax_) ~= 0
	    jwci = [ jwci; ky(1: lky(j+ykmax_),j+ykmax_)+(ykmax_-1)*ny] ;
	    w = [w jacobia_(:,isc(j+ykmax_-1)+1:isc(j+ykmax_))] ;
%	  else
%	    jwci = [] ;
	  end
	end
      end
    end
    jwci = [indnv(jwci,icc1);ncc] ;
    w = [w d1_] ;
    c = zeros(ny,ncc) ;
    c(:,jwci) = w0\w ;
    clear w w0 ;

    junk = fseek(fid,0,1) ;
    fwrite(fid,c','float64') ;
    clear c ;

    it_ = it_ + 1;
    ic = ic + ny ;
    iyr = iyr + ny ;
    i = i - 1 ;
  end
  icr0 = (it_-ykmin_-ykmin_ -1)*ny ;
  while it_ <= iter_+ykmin_
    if broyden_
      d1_ = -feval(fh,y_(iyr));
    else
      jacob(func_name,y_(iyr)) ;
      d1_ = -d1_ ;
    end
    err_f = max(err_f,max(abs(d1_)));
    w0 = jacobia_(:,1:isc(1)) ;
    w = jacobia_(:,isc(1)+1:isc(1+ykmax_)) ;
    j = 1 ;
    while j <= ykmin_
      icr = j1(1:lj1(j),j)-(j-1)*ny ;

      ofs = ((icr0+(j-1)*ny+1)-1)*ncc*8 ;
      junk = fseek(fid,ofs,-1) ;
      c = fread(fid,[ncc,ny],'float64') ;
      c = c' ;

      temp = zeros(ny,ltemp(j)) ;
      if ljwc(j) > 0
	temp(:,jwc(1:ljwc(j),j)) = w ;
      end
      temp(:,jwc1(:,j))=temp(:,jwc1(:,j))-w0*c(icr,1:ncc1) ;
      w = temp ;
      clear temp ;
      d1_ = d1_-w0*c(icr,ncc) ;
      clear c ;
      j = j + 1 ;
      w0 = w(:,1:lj1(j)) ;
      if ykmax_ == 1
	w = jacobia_(:,isc(j+ykmax_-1)+1:isc(j+ykmax_)) ;
      else
	w = w(:,lj1(j)+1:size(w,2)) ;

	if lky(j+ykmax_) > 0
	  w = [w jacobia_(:,isc(j+ykmax_-1)+1:isc(j+ykmax_))] ;
	end
      end
    end
    c = w0\[w d1_] ;
    d1_ = [] ;
    clear w w0 ;
    junk = fseek(fid,0,1) ;
    fwrite(fid,c','float64') ;
    clear c ;
    it_ = it_ + 1 ;
    ic = ic + ny ;
    iyr = iyr + ny ;
    icr0 = icr0 + ny ;
  end
  if ct_ == 1

    ofs = (((it_-ykmin_-2)*ny+1)-1)*ncc*8 ;
    junk = fseek(fid,ofs,-1) ;
    c = fread(fid,[ncc,ny],'float64') ;
    c = c' ;

    for i = 1:ykmax_
      w = tril(triu(ones(ny,ny+ncc1))) ;
      w(:,jwc1(:,ykmin_)) = w(:,jwc1(:,ykmin_))+c(:,1:ncc1) ;
      c = [w(:,ny+1:size(w,2))' c(:,ncc)]/w(:,1:ny) ;

      junk = fseek(fid,0,1) ;
      fwrite(fid,c','float64') ;

      it_ = it_+1 ;
      ic = ic + ny ;

    end
  end
  y_ = reshape(y_,ny,iter_+ykmin_+ykmax_) ;
  if ct_ == 1
    hbacsup = clock ;
    c = bksupk(ny,fid,ncc,icc1) ;
    hbacsup = etime(clock,hbacsup) ;
    c = reshape(c,ny,iter_+ykmax_)' ;
    y(:,1+ykmin_:(iter_+ykmax_+ykmin_)) = y(:,1+ykmin_:(iter_+ykmax_+ykmin_))+slowc_*c' ;
  else
    hbacsup = clock ;
    c = bksupk(ny,fid,ncc,icc1) ;
    hbacsup = etime(clock,hbacsup) ;
    c = reshape(c,ny,iter_)' ;
    y_(:,1+ykmin_:(iter_+ykmin_)) = y_(:,1+ykmin_:(iter_+ykmin_))+slowc_*c' ;
  end

  fclose(fid) ;

  h2 = etime(clock,h2) ;
  [junk,i1] = max(abs(c));
  [junk,i2] = max(junk);
  global lgy_
  disp(['variable ' lgy_(i2,:) ' period ' num2str(i1(i2))])
  err = max(max(abs(c./scalv_'))) ;
  disp ([num2str(iter) '-	err = ' num2str(err)]) ;
  disp (['err_f = ' num2str(err_f)])
  disp (['	Time of this iteration	: ' num2str(h2)]) ;
  if timing_
    disp (['	Back substitution		: ' num2str(hbacsup)]) ;
  end
  if err < dynatol_
    h1 = etime(clock,h1) ;
    fprintf ('\n') ;
    disp (['	Total time of simulation	: ' num2str(h1)]) ;
    fprintf ('\n') ;
    disp (['	Convergence achieved.']) ;
    disp (['-----------------------------------------------------']) ;
    fprintf ('\n') ;
    return ;
  end
end
disp(['WARNING : the maximum number of iterations is reached.']) ;
fprintf ('\n') ;
disp (['-----------------------------------------------------']) ;
return ;

% 2/11/99 MJ took out reshapel






