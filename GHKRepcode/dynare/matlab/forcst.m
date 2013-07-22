function [yf,int_width] = forcst(dr,y0,k,var_list)
  global endo_nbr exo_nbr ykmin_ Sigma_e_ ex_ options_ lgy_
  
  old_options = options_;
  options_.periods = k;
  options_ = set_default_option(options_,'conf_sig',0.9);

  make_ex_;
  yf = simult_(y0,dr,ex_(1:k,:),1);

  [A,B] = kalman_transition_matrix(dr);
  
  sigma_u = B*Sigma_e_*B';
  sigma_y = 0;
  
  for i=1:k
    sigma_y = sigma_y+sigma_u;
    var_yf(i,dr.order_var) = diag(sigma_y(1:endo_nbr,1:endo_nbr))';
    if i == k
      break
    end
    sigma_u = A*sigma_u*A';
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

  fact = qnorm((1-options_.conf_sig)/2,0,1);
  
  int_width = zeros(k,endo_nbr);
  for i=1:endo_nbr
    int_width(:,i) = fact*sqrt(var_yf(:,i));
  end
  

  for i=1:nvar
    my_subplot(i,nvar,2,3,'Forecasts');
    
    plot([-ykmin_+1:0],y0(ivar(i),1:ykmin_),'b-',...
	 [1:k],yf(ivar(i),ykmin_+1:end),'g-',...
	 [1:k],yf(ivar(i),ykmin_+1:end)'+int_width(:,ivar(i)),'g:',...
	 [1:k],yf(ivar(i),ykmin_+1:end)'-int_width(:,ivar(i)),'g:',...
	 [1 k],repmat(dr.ys(ivar(i)),1,2),'r-');
    title(lgy_(ivar(i),:));
  end

  options_ = old_options;








