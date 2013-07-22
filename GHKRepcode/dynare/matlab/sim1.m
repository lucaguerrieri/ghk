% Copyright (C) 2001 Michel Juillard
%
function sim1

global iyp iyf dynatol_ slowc_ maxit_ scalv_ it_
global iy_ ykmin_ ykmax_ xkmin_ xkmax_ ct_ jacobia_ d1_
global iter_ y_ start_simul fname_

ny = size(y_,1) ;
nyp = nnz(iy_(1,:)) ;
nyf = nnz(iy_(3,:)) ;
nrs = ny+nyp+nyf+1 ;
nrc = nyf+1 ;
iyf = find(iy_(3,:)>0) ;
iyp = find(iy_(1,:)>0) ;
isp = [1:nyp] ;
is = [nyp+1:ny+nyp] ;
isf = iyf+nyp ;
isf1 = [nyp+ny+1:nyf+nyp+ny+1] ;
stop = 0 ;

disp (['-----------------------------------------------------']) ;
disp (['MODEL SIMULATION :']) ;
fprintf('\n') ;

if isempty(start_simul)
  it_init = 2 ;
else
  it_init = start_simul;
end

h1 = clock ;
for iter = 1:maxit_
	h2 = clock ;

	if ct_ == 0
		c = zeros(ny*iter_,nrc) ;
	else
		c = zeros(ny*(iter_+1),nrc) ;
	end

	it_ = it_init ;
	z = [y_(iyp,it_-1) ; y_(:,it_) ; y_(iyf,it_+1)] ;
	jacob ([fname_ '_ff'],z)
	jacobia_ = [jacobia_ -d1_] ;
	ic = [1:ny] ;
	icp = iyp ;
	c (ic,:) = jacobia_(:,is)\jacobia_(:,isf1) ;
	for it_ = it_init+1:iter_+1
		z = [y_(iyp,it_-1) ; y_(:,it_) ; y_(iyf,it_+1)] ;
		jacob ([fname_ '_ff'],z)
		jacobia_ = [jacobia_ -d1_] ;
		jacobia_(:,[isf nrs]) = jacobia_(:,[isf nrs])-jacobia_(:,isp)*c(icp,:) ;
		ic = ic + ny ;
		icp = icp + ny ;
		c (ic,:) = jacobia_(:,is)\jacobia_(:,isf1) ;
	end

	if ct_ == 1
		s = eye(ny) ;
		s(:,isf) = s(:,isf)+c(ic,1:nyf) ;
		ic = ic + ny ;
		c(ic,nrc) = s\c(:,nrc) ;
		c = bksup1(ny,nrc,iyf,c) ;
		c = reshape(c,ny,iter_+1) ;
		y_(:,it_init:iter_+2) = y_(:,it_init:iter_+2)+slowc_*c ;
	else
		c = bksup1(ny,nrc,iyf,c) ;
		c = reshape(c,ny,iter_) ;
		y_(:,it_init:iter_+1) = y_(:,it_init:iter_+1)+slowc_*c ;
	end

	err = max(max(abs(c./scalv_')));
	disp([num2str(iter) ' -	err = ' num2str(err)]) ;
	disp(['	Time of iteration 	:' num2str(etime(clock,h2))]) ;

	if err < dynatol_
		stop = 1 ;
		fprintf('\n') ;
		disp(['	Total time of simulation 	:' num2str(etime(clock,h1))]) ;
		fprintf('\n') ;
		disp(['	Convergency obtained.']) ;
		fprintf('\n') ;
		break
	end
end

if ~ stop
	fprintf('\n') ;
	disp(['	Total time of simulation 	:' num2str(etime(clock,h1))]) ;
	fprintf('\n') ;
	disp(['WARNING : maximum number of iterations is reached (modify maxit_).']) ;
	fprintf('\n') ;
end
disp (['-----------------------------------------------------']) ;
return ;

% 08/24/01 MJ added start_simul