% Copyright (C) 2001 Michel Juillard
%
function osr(var_list,params,W)
  global iter_ ys_ dr_ y_ dr_ exo_nbr lgy_ lgx_ Sigma_e_ ykmin_ ykmax_ ...
      endo_nbr options_ lgx_orig_ord_
  
  options_old = options_;
  options_.order = 1;
  options_.linear = 1;
  options_ = set_default_option(options_,'ar',5);
  options_ = set_default_option(options_,'irf',40);
  options_ = set_default_option(options_,'dr_algo',0);
  options_ = set_default_option(options_,'simul_algo',0);
  options_ = set_default_option(options_,'drop',100);
  options_ = set_default_option(options_,'replic',1);
  options_ = set_default_option(options_,'nomoments',0);
  options_ = set_default_option(options_,'nocorr',0);
  options_ = set_default_option(options_,'simul_seed',[]);
  options_ = set_default_option(options_,'hp_filter',0);
  options_ = set_default_option(options_,'hp_ngrid',512);
  options_ = set_default_option(options_,'simul',0);
  options_ = set_default_option(options_,'periods',1);
  options_ = set_default_option(options_,'TeX',0);
  options_ = set_default_option(options_,'noprint',0);
  if options_.simul & ~isempty(iter_) & options_.periods == 0
    options_.periods = iter_;
  end

  iter_ = max(options_.periods,1);
  options_.periods = iter_;
  
  make_ex_;

  disp(' ')
  disp('OPTIMAL SIMPLE RULE')
  disp(' ')
  [dr_,info] = osr1(params,W);
  if info > 0
    options_ = options_old;
    return
  elseif options_.noprint ~= 1
      disp('MODEL SUMMARY')
      disp(' ')
      disp(['  Number of variables:         ' int2str(endo_nbr)])
      disp(['  Number of stochastic shocks: ' int2str(exo_nbr)])
      disp(['  Number of state variables:   ' ...
	    int2str(length(find(dr_.kstate(:,2) <= ykmin_+1)))])
      disp(['  Number of jumpers:           ' ...
	    int2str(length(find(dr_.kstate(:,2) == ykmin_+2)))])
      disp(['  Number of static variables:  ' int2str(dr_.nstatic)])
      my_title='MATRIX OF COVARIANCE OF EXOGENOUS SHOCKS';
      labels = deblank(lgx_);
      headers = strvcat('Variables',labels);
      lh = size(labels,2)+2;
      table(my_title,headers,labels,Sigma_e_,lh,10,6);
      disp(' ')
      disp_dr(dr_,options_.order,var_list);
      if options_.order == 1 & options_.simul == 0 & options_.nomoments == 0
	disp_th_moments(dr_,var_list);
      elseif options_.simul == 1 | options_.nomoments == 0
	if options_.periods == 0
	  error('OSR error: number of periods for the simulation isn''t specified')
	end
	if options_.periods < options_.drop
	  disp(['OSR error: The horizon of simulation is shorter' ...
		' than the number of observations to be DROPed'])
	  return
	end
	
	y_ = simult(repmat(dr_.ys,1,ykmin_),dr_);
	dyn2vec;
	if options_.nomoments == 0
	  disp_moments(y_,var_list);
	end
      end
      
      n = size(var_list,1);
      if n == 0
	n = length(dr_.order_var);
	ivar = [1:n]';
	var_list = lgy_;
      else
	ivar=zeros(n,1);
	for i=1:n
	  i_tmp = strmatch(var_list(i,:),lgy_,'exact');
	  if isempty(i_tmp)
	    error (['One of the specified variables does not exist']) ;
	  else
	    ivar(i) = i_tmp;
	  end
	end
      end

      if n < 13 & options_.irf > 0
	if n == 1
	  nr = 1;
	  nc = 1;
	elseif n == 2
	  nr = 1;
	  nc = 2;
	elseif n <= 4
	  nr = 2;
	  nc = 2;
	elseif n <= 6
	  nr = 2;
	  nc = 3;
	elseif n <= 9
	  nr = 3;
	  nc = 3;
	elseif n <= 12
	  nr = 3;
	  nc = 4;
	end
	olditer = iter_;
	if options_.order == 1
	  options_.replic = 1;
	else
	  if options_.replic == 0
	    options_.replic = 50;
	  end
	end
	SS(lgx_orig_ord_,lgx_orig_ord_)=Sigma_e_+1e-14*eye(exo_nbr);
	cs = transpose(chol(SS));

	for i = 1:exo_nbr
	  if SS(i,i) > 1e-13
	    figure('Name',['Shock to ' lgx_(i,:)]);
	    y=irf(dr_,cs(lgx_orig_ord_,i), options_.irf, options_.drop, options_.replic, options_.order);
	    for j = 1:n
	      subplot(nr,nc,j);
	      plot([y(ivar(j),:)']);
	      title(var_list(j,:),'Interpreter','none');
	      assignin('base',[deblank(var_list(j,:)) '_' deblank(lgx_(i,:))],y(ivar(j),:)');
	    end
	  end
	  
	end
	iter_ = olditer;
      end
  end
  options_ = options_old;
  
% 01/10/01 FC dr_ and y_ made global
% 02/20/01 MJ ys_ removed from calling sequence for simult (all in dr_)
% 02/23/01 MJ added dyn2vec()
% 06/24/01 MJ steady -> steady_
% 09/24/01 MJ dr_ made global
% 08/28/02 MJ added var_list
% 10/09/02 MJ no simulation and theoretical moments for order 1 
% 10/14/02 MJ added plot of IRFs
% 10/30/02 MJ options_ are now a structure
% 01/01/03 MJ added dr_algo
% 01/09/03 MJ set default values for options_ (correct absence of autocorr
%             when order == 1)
% 01/12/03 MJ removed call to steady_ as already checked in resol()
% 01/31/03 MJ make IRF global with varname_shockname
% 02/09/03 MJ ys_ reset with value declared in initval after computations
% 02/18/03 MJ removed above change. ys_ shouldn't be affected by
%             computations in this function
%             new option SIMUL computes a stochastic simulation and save
%             results in y_ and via dyn2vec
% 04/03/03 MJ corrected bug for simulation with ykmin_ > 1