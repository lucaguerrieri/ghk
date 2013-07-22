% Copyright (C) 2005 Michel Juillard
%

function [y_,int_width]=simultexdet(y0,dr,ex_,ex_det, iorder,var_list)
  global endo_nbr ykmin_ xkmin_ it_ options_ iy_ M_ exe_det_ Sigma_e_ lgy_

  iter = size(ex_,1)-xkmin_;
  nx = size(dr.ghu,2);
  y_ = zeros(size(y0,1),iter+ykmin_);
  y_(:,1:ykmin_) = y0;
  k1 = [ykmin_:-1:1];
  k2 = dr.kstate(find(dr.kstate(:,2) <= ykmin_+1),[1 2]);
  k2 = k2(:,1)+(ykmin_+1-k2(:,2))*endo_nbr;
  k3 = iy_(1:ykmin_,:)';
  k3 = find(k3(:));
  k4 = dr.kstate(find(dr.kstate(:,2) < ykmin_+1),[1 2]);
  k4 = k4(:,1)+(ykmin_+1-k4(:,2))*endo_nbr;
  
  old_options = options_;
  options_ = set_default_option(options_,'simul_algo',0);
  options_ = set_default_option(options_,'conf_sig',0.9);
  
  if options_.simul_algo == 1
    o1 = dr.nstatic+1;
    o2 = dr.nstatic+dr.npred;
    o3 = o2-dr.nboth+1;
    [junk, k5] = sort(dr.order_var(o1:o2));
    [junk, k6] = sort(dr.order_var(o3:end));
  end
  
  nvar = size(var_list,1);
  if nvar == 0
    nvar = endo_nbr;
    ivar = [1:nvar];
  else
    ivar=zeros(nvar,1);
    for i=1:nvar
      i_tmp = strmatch(var_list(i,:),lgy_,'exact');
      if isempty(i_tmp)
	disp(var_list(i,:));
	error (['One of the variable specified does not exist']) ;
      else
	ivar(i) = i_tmp;
      end
    end
  end

  if iorder == 1
    for i = ykmin_+1: iter+ykmin_
      tempx1 = y_(dr.order_var,k1);
      tempx2 = tempx1-repmat(dr.ys(dr.order_var),1,ykmin_);
      tempx = tempx2(k2);
      if options_.simul_algo == 0
	y_(dr.order_var,i) = dr.ys(dr.order_var)+dr.ghx*tempx+dr.ghu* ...
	    ex_(i+xkmin_-ykmin_,:)';
	for j=1:min(iter+ykmin_-i-M_.ex_det_length+1,M_.ex_det_length)
	  y_(dr.order_var,i) = y_(dr.order_var,i) + dr.ghud{j}*(ex_det(i+j-1,:)'-exe_det_');
	end
      elseif options_.simul_algo == 1
	it_ = i;
	m = dr.ys(dr.order_var);
	[y_(:,i), check] = dynare_solve('ff_simul1',y_(:,i-1),tempx1(k3), ...
					m(o3:end),tempx(k4),o1,o2,o3,k6);
      end
	
      k1 = k1+1;
    end
  elseif iorder == 2
    for i = ykmin_+1: iter+ykmin_
      tempx1 = y_(dr.order_var,k1);
      tempx2 = tempx1-repmat(dr.ys(dr.order_var),1,ykmin_);
      tempx = tempx2(k2);
      tempu = ex_(i+xkmin_-ykmin_,:)';
      tempuu = kron(tempu,tempu);
      if options_.simul_algo == 0
	tempxx = kron(tempx,tempx);
	tempxu = kron(tempx,tempu);
	y_(dr.order_var,i) = dr.ys(dr.order_var)+dr.ghs2/2+dr.ghx*tempx+ ...
	    dr.ghu*tempu+0.5*(dr.ghxx*tempxx+dr.ghuu*tempuu)+dr.ghxu* ...
	    tempxu;
	for j=1:min(iter+ykmin_-i-M_.ex_det_length+1,M_.ex_det_length)
	  tempud = ex_det(i+j-1,:)'-exe_det_';
	  tempudud = kron(tempud,tempud);
	  tempxud = kron(tempx,tempud);
	  tempuud = kron(tempu,tempud);
	  y_(dr.order_var,i) = y_(dr.order_var,i) + dr.ghud{j}*tempud + ...
	      dr.ghxud{j}*tempxud + dr.ghuud{j}*tempuud + ...
	      0.5*dr.ghudud{j,j}*tempudud;
	  for k=1:j-1
	    tempudk = ex_det(i+k-1,:)'-exe_det_';
	    tempududk = kron(tempudk,tempud);
	    y_(dr.order_var,i) = y_(dr.order_var,i) + ...
		dr.ghudud{k,j}*tempududk;
	  end
	end
      elseif options_.simul_algo == 1
	it_ = i;
	m = dr.ys(dr.order_var)+dr.ghs2/2;
	tempx1 = y_(:,k1);
	[y_(:,i), check] = dynare_solve('ff_simul2',y_(:,i-1),tempx1(k3), ...
					m(o3:end),tempx(k4),o1,o2,o3,k6);
      end
      k1 = k1+1;
    end
  end

  [A,B] = kalman_transition_matrix(dr);
  
  sigma_u = B*Sigma_e_*B';
  sigma_y = 0;
  
  for i=1:iter
    sigma_y = sigma_y+sigma_u;
    var_yf(i,dr.order_var) = diag(sigma_y(1:endo_nbr,1:endo_nbr))';
    if i == iter
      break
    end
    sigma_u = A*sigma_u*A';
  end

  fact = qnorm((1-options_.conf_sig)/2,0,1);
  
  int_width = zeros(iter,endo_nbr);
  for i=1:endo_nbr
    int_width(:,i) = fact*sqrt(var_yf(:,i));
  end
  
  for i=1:nvar
    my_subplot(i,nvar,2,3,'Forecasts');
    
    plot([-ykmin_+1:0],y0(ivar(i),1:ykmin_),'b-',...
	 [1:iter],y_(ivar(i),ykmin_+1:end),'g-',...
	 [1:iter],y_(ivar(i),ykmin_+1:end)'+int_width(:,ivar(i)),'g:',...
	 [1:iter],y_(ivar(i),ykmin_+1:end)'-int_width(:,ivar(i)),'g:',...
	 [1 iter],repmat(dr.ys(ivar(i)),1,2),'r-');
    title(lgy_(ivar(i),:));
  end

  options_ = old_options;