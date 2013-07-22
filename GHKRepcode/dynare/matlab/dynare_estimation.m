function dynare_estimation(varlist)
global options_ bayestopt_ estim_params_ lgy_ lgx_ Sigma_e_ exo_nbr endo_nbr
global lgy_ lgx_ fname_ oo_  trend_coeff_  iy_ ykmin_ lgx_TeX_ lgy_TeX_
global dsge_prior_weight ys_



options_.varlist = varlist;

options_.lgyidx2varobs = zeros(size(lgy_,1),1);
for i = 1:size(lgy_,1)
  tmp = strmatch(deblank(lgy_(i,:)),options_.varobs,'exact');
  if ~isempty(tmp)
    options_.lgyidx2varobs(i,1) = tmp;
  end           
end

options_ = set_default_option(options_,'first_obs',1);
options_ = set_default_option(options_,'prefilter',0);
options_ = set_default_option(options_,'presample',0);
options_ = set_default_option(options_,'lik_algo',1);
options_ = set_default_option(options_,'lik_init',1);
options_ = set_default_option(options_,'nograph',0);
options_ = set_default_option(options_,'mh_conf_sig',0.90);
options_ = set_default_option(options_,'mh_replic',20000);
options_ = set_default_option(options_,'mh_drop',0.5);
options_ = set_default_option(options_,'mh_jscale',0.2);
options_ = set_default_option(options_,'mh_init_scale',2*options_.mh_jscale);
options_ = set_default_option(options_,'mode_file','');
options_ = set_default_option(options_,'mode_compute',4);
options_ = set_default_option(options_,'mode_check',0);
options_ = set_default_option(options_,'prior_trunc',1e-10);
options_ = set_default_option(options_,'load_mh_file',0);
options_ = set_default_option(options_,'loglinear',0);
options_ = set_default_option(options_,'unit_root_vars',[]);
options_ = set_default_option(options_,'XTick',[]);
options_ = set_default_option(options_,'XTickLabel',[]);
options_ = set_default_option(options_,'bayesian_irf',0);
options_ = set_default_option(options_,'bayesian_th_moments',0);
options_ = set_default_option(options_,'TeX',0);
options_ = set_default_option(options_,'irf',40);
options_ = set_default_option(options_,'relative_irf',0);
options_ = set_default_option(options_,'order',1);
options_ = set_default_option(options_,'ar',5);
options_ = set_default_option(options_,'dr_algo',0);
options_ = set_default_option(options_,'linear',1);
options_ = set_default_option(options_,'drop',0);
options_ = set_default_option(options_,'replic',1);
options_ = set_default_option(options_,'hp_filter',0);
options_ = set_default_option(options_,'forecast',0);
options_ = set_default_option(options_,'smoother',0);
options_ = set_default_option(options_,'moments_varendo',0);
options_ = set_default_option(options_,'filtered_vars',0);
options_ = set_default_option(options_,'kalman_algo',1);
options_ = set_default_option(options_,'kalman_tol',10^(-12));
options_ = set_default_option(options_,'diffuse_d',[]);
options_ = set_default_option(options_,'nk',1);
options_ = set_default_option(options_,'varlag',4);
options_ = set_default_option(options_,'Opt6Iter',3);
options_ = set_default_option(options_,'Opt6Numb',40000);
options_ = set_default_option(options_,'ExcludedParams',[]);
options_ = set_default_option(options_,'ParamSubSet','All');

optim_options = optimset('display','iter','LargeScale','off',...
             'MaxFunEvals',100000,'TolFun',1e-8,'TolX',1e-6);

if isfield(options_,'optim_opt')
  eval(['optim_options = optimset(optim_options,' options_.optim_opt ');']);
end

n_varobs = size(options_.varobs,1);
dr = set_state_space([]);

% load data
rawdata = read_variables(options_.datafile,options_.varobs,[]);

k = [];
k1 = [];
for i=1:n_varobs
  k = [k strmatch(deblank(options_.varobs(i,:)),lgy_(dr.order_var,:), ...
          'exact')];
  k1 = [k1 strmatch(deblank(options_.varobs(i,:)),lgy_, 'exact')];
end

% check initial parameter values and creates bayestopt_
pnames=['     ';'beta ';'gamm ';'norm ';'invg ';'unif ';'invg2'];
[xparam1,estim_params_,bayestopt_,lb,ub]=set_prior(estim_params_);
if any(bayestopt_.pshape > 0)
  plot_priors
else
  options_.mh_replic = 0;
end

bayestopt_.mf = k;
bayestopt_.mfys = k1;

% set options for recursive forecast if necessary
options_ = set_default_option(options_,'nobs',size(rawdata,1)-options_.first_obs+1);
if options_.nobs(1) == options_.nobs(end) % no recursive estimation
  options_ = set_default_option(options_,'mh_nblck',2); 
  options_ = set_default_option(options_,'nodiagnostic',0);
else
  options_ = set_default_option(options_,'mh_nblck',1); 
  options_ = set_default_option(options_,'nodiagnostic',1);
end

bounds = prior_bounds(bayestopt_);
bounds(:,1)=max(bounds(:,1),lb);
bounds(:,2)=min(bounds(:,2),ub);
if any(xparam1 < bounds(:,1)) | any(xparam1 > bounds(:,2))
  find(xparam1 < bounds(:,1))
  find(xparam1 > bounds(:,2))
  error('Initial parameter values are outside parameter bounds')
end
lb = bounds(:,1);
ub = bounds(:,2);
bayestopt_.lb = lb;
bayestopt_.ub = ub;

if isempty(trend_coeff_)
  bayestopt_.with_trend = 0;
else
  bayestopt_.with_trend = 1;
  bayestopt_.trend_coeff = {};
  for i=1:n_varobs
    if i > length(trend_coeff_) | isempty(trend_coeff_{i})
      bayestopt_.trend_coeff{i} = '0';
    else
      bayestopt_.trend_coeff{i} = trend_coeff_{i};
    end
  end
end

bayestopt_.penalty=1e8;             % penalty 


nvx = estim_params_.nvx;
nvn = estim_params_.nvn;
ncx = estim_params_.ncx;
ncn = estim_params_.ncn;
np  = estim_params_.np ;
nx = nvx+nvn+ncx+ncn+np;

check_model;

%static solver
if exist([fname_ '_steadystate'])
  bayestopt_.static_solve = [fname_ '_steadystate'];
else
  bayestopt_.static_solve = 'dynare_solve';
end

% Initialization with unit-root variables
if ~isempty(options_.unit_root_vars)
  n_ur = length(options_.unit_root_vars);
  i_ur = zeros(n_ur,1);
  for i=1:n_ur
    i1 = strmatch(options_.unit_root_vars{i},lgy_(dr.order_var,:),'exact');
    if isempty(i1)
      error('Undeclared variable in unit_root_vars statement')
    end
    i_ur(i) = i1;
  end
  i_stable = ones(endo_nbr,1);
  i_stable(i_ur) = zeros(n_ur,1);
  i_stable = find(i_stable);
  if ykmin_ > 1
    l1 = flipud([cumsum(iy_(1:ykmin_-1,dr.order_var),1);ones(1, ...
                          endo_nbr)]);
    n1 = nnz(l1);
    bayestopt_.Pinf = zeros(n1,n1);
    l2 = find(l1');
    l3 = zeros(endo_nbr,ykmin_);
    l3(i_ur,:) = l1(:,i_ur)';
    l3 = l3(:);
    i_ur1 = find(l3(l2));
    bayestopt_.Pinf(i_ur1,i_ur1) = diag(ones(1,length(i_ur1)));
    l3 = zeros(endo_nbr,ykmin_);
    l3(i_stable,:) = l1(:,i_stable)';
    l3 = l3(:);
    bayestopt_.i_T_var_stable = find(l3(l2));
  else
    n1 = endo_nbr;
    bayestopt_.Pinf = zeros(n1,n1);
    bayestopt_.Pinf(i_ur,i_ur) = diag(ones(1,length(i_ur)));
    l1 = ones(endo_nbr,1);
    l1(i_ur,:) = zeros(length(i_ur),1);
    bayestopt_.i_T_var_stable = find(l1);
  end
  options_.lik_init = 3;
else
  i_stable = [1:endo_nbr]';
end % if ~isempty(options_.unit_root_vars)
bayestopt_.i_var_stable = i_stable;

if isempty(options_.datafile)
  error('ESTIMATION: datafile option is missing')
end

if isempty(options_.varobs)
  error('ESTIMATION: VAROBS is missing')
end


% if jscale isn't specified for an estimated parameter, use
% global option options_.jscale, set to 0.2, by default
k = find(isnan(bayestopt_.jscale));
bayestopt_.jscale(k) = options_.mh_jscale;

% transform data 
if options_.loglinear == 1
  rawdata = log(rawdata);
end
if ~isreal(rawdata)
  error(['There are complex values in the data. Probably  a wrong' ...
     ' transformation'])
end

% load a previous posterior mode if any
if length(options_.mode_file) > 0
  eval(['load ' options_.mode_file ';']);
  if exist('Scale')
    options_.mh_jscale = Scale;
    bayestopt_.jscale = ones(length(xparam1),1)*Scale;
  end
  if exist('PostMean')
    MeanPar = PostMean;
  end  
end

% compute sample moments if needed (var-dsge)
if ~isempty(strmatch('dsge_prior_weight',estim_params_.param_names)) | ~isempty(dsge_prior_weight)
    evalin('base',['[mYY,mXY,mYX,mXX,Ydata,Xdata] = ' ...
	   'VarSampleMoments(options_.first_obs,options_.first_obs+options_.nobs-1,options_.varlag,-1);'])
end

iter1 = 1;
for iter1 = 1:length(options_.nobs)
  gend = options_.nobs(iter1);
  if options_.prefilter == 1
    bayestopt_.mean_varobs = mean(rawdata(options_.first_obs+(0:gend-1),:),1);
    data = rawdata(options_.first_obs+(0:gend-1),:)'-...
	   repmat(bayestopt_.mean_varobs',1,gend);
  else
    data = rawdata(options_.first_obs+(0:gend-1),:)';
  end  
  initial_estimation_checks(xparam1,gend,data);
  if options_.mode_compute > 0
    if isempty(strmatch('dsge_prior_weight',estim_params_.param_names)) & isempty(dsge_prior_weight)
      fh=str2func('DsgeLikelihood');
    else
      fh=str2func('DsgeVarLikelihood');        
    end
    if options_.mode_compute == 1
      if isempty(strmatch('dsge_prior_weight',estim_params_.param_names)) & isempty(dsge_prior_weight)
	[xparam1,fval,exitflag,output,lamdba,grad,hessian_fmincon] = ...
            fmincon(fh,xparam1,[],[],[],[],lb,ub,[],optim_options,gend,data);
      else
	[xparam1,fval,exitflag,output,lamdba,grad,hessian_fmincon] = ...
            fmincon(fh,xparam1,[],[],[],[],lb,ub,[],optim_options,gend);
      end
    elseif options_.mode_compute == 2
      %   asamin('set','maximum_cost_repeat',0);
      if isempty(strmatch('dsge_prior_weight',estim_params_.param_names)) & isempty(dsge_prior_weight)
	[fval,xparam1,grad,hessian_asamin,exitflag] = ...
	    asamin('minimize','DsgeLikelihood',xparam1,lb,ub,-ones(size(xparam1)),gend,data);   
      else
	[fval,xparam1,grad,hessian_asamin,exitflag] = ...
	    asamin('minimize','DsgeVarLikelihood',xparam1,lb,ub,-ones(size(xparam1)),gend);   
      end  
    elseif options_.mode_compute == 3
      if isempty(strmatch('dsge_prior_weight',estim_params_.param_names)) & isempty(dsge_prior_weight)
	[xparam1,fval,exitflag] = fminunc(fh,xparam1,optim_options,gend,data);
      else
	[xparam1,fval,exitflag] = fminunc(fh,xparam1,optim_options,gend);
      end  
    elseif options_.mode_compute == 4
      H0 = 1e-4*eye(nx);
      crit = 1e-7;
      nit = 1000;
      verbose = 2;
      if isempty(strmatch('dsge_prior_weight',estim_params_.param_names)) & isempty(dsge_prior_weight)
	[fval,xparam1,grad,hessian_csminwel,itct,fcount,retcodehat] = ...
	    csminwel('DsgeLikelihood',xparam1,H0,[],crit,nit,gend,data);
	disp(sprintf('Objective function at mode: %f',fval))
	disp(sprintf('Objective function at mode: %f',DsgeLikelihood(xparam1,gend,data)))
      else
	[fval,xparam1,grad,hessian_csminwel,itct,fcount,retcodehat] = ...
	    csminwel('DsgeVarLikelihood',xparam1,H0,[],crit,nit,gend);
	disp(sprintf('Objective function at mode: %f',fval))
	disp(sprintf('Objective function at mode: %f',DsgeVarLikelihood(xparam1,gend)))
      end
    elseif options_.mode_compute == 5
      flag = 0;
      if isempty(strmatch('dsge_prior_weight',estim_params_.param_names)) & isempty(dsge_prior_weight)
	[xparam1, hh, gg, fval] = newrat('DsgeLikelihood',xparam1,[],[],flag,gend,data);
	eval(['save ' fname_ '_mode xparam1 hh gg fval;']);
      else
	[xparam1, hh, gg, fval] = newrat('DsgeVarLikelihood',xparam1,[],[],flag,gend);
	eval(['save ' fname_ '_mode xparam1 hh gg fval;']);
      end
    elseif options_.mode_compute == 6% New routine!
      if isempty(strmatch('dsge_prior_weight',estim_params_.param_names)) & isempty(dsge_prior_weight)
	fval = DsgeLikelihood(xparam1,gend,data);
      else
	fval = DsgeVarLikelihood(xparam1,gend);
      end      
      OldMode = fval;
      if ~exist('MeanPar')
	MeanPar = xparam1;
      end  
      for i=1:options_.Opt6Iter
	if exist('hh')
	  CovJump = inv(hh);
	else% The covariance matrix is initialized with the prior
	  % covariance (a diagonal matrix)
	  stdev = bayestopt_.pstdev;
	  indx = find(isinf(stdev));
	  stdev(indx) = 10*ones(length(indx),1);
	  CovJump = diag(stdev).^2;
	end
	OldPostVar = CovJump;
	if i == 1;
	  if options_.Opt6Iter > 1
	    [xparam1,PostVar,Scale,PostMean] = ...
		metropolis99(xparam1,gend,data,bounds,options_.Opt6Numb,options_.mh_jscale,'',MeanPar,CovJump);
	  else
	    [xparam1,PostVar,Scale,PostMean] = ...
		metropolis99(xparam1,gend,data,bounds,options_.Opt6Numb,options_.mh_jscale,'LastCall',MeanPar,CovJump);
	  end
	  options_.mh_jscale = Scale;
	  mouvement = max(max(abs(PostVar-OldPostVar)));
	  if isempty(strmatch('dsge_prior_weight',estim_params_.param_names)) & isempty(dsge_prior_weight)
	    fval = DsgeLikelihood(xparam1,gend,data);
	  else
	    fval = DsgeVarLikelihood(xparam1,gend);
	  end
	  disp(['Change in the covariance matrix = ' num2str(mouvement) '.'])
	  disp(['Mode improvement = ' num2str(abs(OldMode-fval))])
	  OldMode = fval;
	else
	  OldPostVar = PostVar;
	  if i<options_.Opt6Iter
	    [xparam1,PostVar,Scale,PostMean] = ...
		metropolis99(xparam1,gend,data,bounds,options_.Opt6Numb,options_.mh_jscale,'',PostMean,PostVar);
	  else
	    [xparam1,PostVar,Scale,PostMean] = ...
		metropolis99(xparam1,gend,data,bounds,options_.Opt6Numb,options_.mh_jscale,'LastCall',PostMean,PostVar);
	  end
	  options_.mh_jscale = Scale;
	  mouvement = max(max(abs(PostVar-OldPostVar)));
	  if isempty(strmatch('dsge_prior_weight',estim_params_.param_names)) & isempty(dsge_prior_weight)
	    fval = DsgeLikelihood(xparam1,gend,data);
	  else
	    fval = DsgeVarLikelihood(xparam1,gend);
	  end
	  disp(['Change in the covariance matrix = ' num2str(mouvement) '.'])
	  disp(['Mode improvement = ' num2str(abs(OldMode-fval))])
	  OldMode = fval;
	end
	options_.mh_jscale = Scale;
	bayestopt_.jscale = ones(length(xparam1),1)*Scale;
      end
      hh = inv(PostVar);
      if isempty(strmatch('dsge_prior_weight',estim_params_.param_names)) & isempty(dsge_prior_weight)
	fval = DsgeLikelihood(xparam1,gend,data);
      else
	fval = DsgeVarLikelihood(xparam1,gend);
      end
    end
    if options_.mode_compute ~= 5 & options_.mode_compute ~= 6 
      if isempty(strmatch('dsge_prior_weight',estim_params_.param_names)) & isempty(dsge_prior_weight)
	hh = reshape(hessian('DsgeLikelihood',xparam1,gend,data),nx,nx);
	eval(['save ' fname_ '_mode xparam1 hh fval;']);
      else
	hh = reshape(hessian('DsgeVarLikelihood',xparam1,gend),nx,nx);
	eval(['save ' fname_ '_mode xparam1 hh fval;']);
      end
    end
    if options_.mode_compute ~= 6
      eval(['save ' fname_ '_mode xparam1 hh;']);
    elseif options_.mode_compute == 6
      eval(['save ' fname_ '_mode xparam1 hh Scale PostMean;']);
    end
  end

  if options_.mode_check == 1
    mode_check(xparam1,0,hh,gend,data,lb,ub);
  end

  if options_.mode_compute ~= 6
    hh = generalized_cholesky(hh);
    invhess = inv(hh);
  else
    invhess = PostVar;
  end
  stdh = sqrt(diag(invhess));
  if any(bayestopt_.pshape > 0)
    disp(' ')
    disp('RESULTS FROM POSTERIOR MAXIMIZATION')
    tstath = zeros(nx,1);
    for i = 1:nx
      tstath(i) = abs(xparam1(i))/stdh(i);
    end
    tit1 = sprintf('%10s %7s %8s %7s %6s %4s %6s\n',' ','prior mean', ...
           'mode','s.d.','t-stat','prior','pstdev');
    if np
      ip = nvx+nvn+ncx+ncn+1;
      disp('parameters')
      disp(tit1)
      for i=1:np
        disp(sprintf('%12s %7.3f %8.4f %7.4f %7.4f %4s %6.4f', ...
                     deblank(estim_params_.param_names(i,:)), ...
                     bayestopt_.pmean(ip),xparam1(ip),stdh(ip),tstath(ip), ...
                     pnames(bayestopt_.pshape(ip)+1,:), ...
                     bayestopt_.pstdev(ip)));
        eval(['oo_.posterior_mode.parameters.' deblank(estim_params_.param_names(i,:)) ' = xparam1(ip);']);
        eval(['oo_.posterior_std.parameters.' deblank(estim_params_.param_names(i,:)) ' = stdh(ip);']); 
        ip = ip+1;
      end
    end
    if nvx
      ip = 1;
      disp('standard deviation of shocks')
      disp(tit1)
      for i=1:nvx
        k = estim_params_.var_exo(i,1);
        disp(sprintf('%12s %7.3f %8.4f %7.4f %7.4f %4s %6.4f', ...
                     deblank(lgx_(k,:)),bayestopt_.pmean(ip),xparam1(ip), ...
                     stdh(ip),tstath(ip),pnames(bayestopt_.pshape(ip)+1,:), ...
                     bayestopt_.pstdev(ip))); 
        Sigma_e_(k,k) = xparam1(ip)*xparam1(ip);
        eval(['oo_.posterior_mode.shocks_std.' deblank(lgx_(k,:)) ' = xparam1(ip);']);
        eval(['oo_.posterior_std.shocks_std.' deblank(lgx_(k,:)) ' = stdh(ip);']); 
        ip = ip+1;
      end
    end
    if nvn
      disp('standard deviation of measurement errors')
      disp(tit1)
      ip = nvx+1;
      for i=1:nvn
        disp(sprintf('%12s %7.3f %8.4f %7.4f %7.4f %4s %6.4f', ...
                     deblank(options_.varobs(estim_params_.var_endo(i,1),: ...
                                             )),bayestopt_.pmean(ip), ...
                     xparam1(ip),stdh(ip),tstath(ip), ...
                     pnames(bayestopt_.pshape(ip)+1,:), ...
                     bayestopt_.pstdev(ip)));
        eval(['oo_.posterior_mode.measurement_errors_std.' deblank(options_.varobs(estim_params_.var_endo(i,1),:)) ' = xparam1(ip);']);
        eval(['oo_.posterior_std.measurement_errors_std.' deblank(options_.varobs(estim_params_.var_endo(i,1),:)) ' = stdh(ip);']); 
        ip = ip+1;
      end
    end
    if ncx
      disp('correlation of shocks')
      disp(tit1)
      ip = nvx+nvn+1;
      for i=1:ncx
        k1 = estim_params_.corrx(i,1);
        k2 = estim_params_.corrx(i,2);
        name = [deblank(lgx_(k1,:)) ',' deblank(lgx_(k2,:))];
        disp(sprintf('%12s %7.3f %8.4f %7.4f %7.4f %4s %6.4f', name, ...
                     bayestopt_.pmean(ip),xparam1(ip),stdh(ip),tstath(ip),  ...
                     pnames(bayestopt_.pshape(ip)+1,:), bayestopt_.pstdev(ip)));
        Sigma_e_(k1,k2) = xparam1(ip)*sqrt(Sigma_e_(k1,k1)*Sigma_e_(k2,k2));
        Sigma_e_(k2,k1) = Sigma_e_(k1,k2);
        eval(['oo_.posterior_mode.shocks_corr.' deblank(lgx_(k1,:)) '_' deblank(lgx_(k2,:)) ' = xparam1(ip);']);
        eval(['oo_.posterior_std.shocks_corr.' deblank(lgx_(k1,:)) '_' deblank(lgx_(k2,:)) ' = stdh(ip);']); 
        ip = ip+1;
      end
    end
    if ncn
      disp('correlation of measurement errors')
      disp(tit1)
      ip = nvx+nvn+ncx+1;
      for i=1:ncn
        k1 = estim_params_.corrn(i,1);
        k2 = estim_params_.corrn(i,2);
        name = [deblank(lgy_(k1,:)) ',' deblank(lgy_(k2,:))];
        disp(sprintf('%12s %7.3f %8.4f %7.4f %7.4f %4s %6.4f', name, ...
                     bayestopt_.pmean(ip),xparam1(ip),stdh(ip),tstath(ip), ...
                     pnames(bayestopt_.pshape(ip)+1,:), bayestopt_.pstdev(ip)));
        eval(['oo_.posterior_mode.measurement_errors_corr.' deblank(lgy_(k1,:)) '_' deblank(lgy_(k2,:)) ' = xparam1(ip);']);
        eval(['oo_.posterior_std.measurement_errors_corr.' deblank(lgy_(k1,:)) '_' deblank(lgy_(k2,:)) ' = stdh(ip);']); 
        ip = ip+1;
      end
    end
    %%% Laplace approximation to the marginal log density:
    if isempty(strmatch('dsge_prior_weight',estim_params_.param_names)) & isempty(dsge_prior_weight)
        md_Laplace = .5*size(xparam1,1)*log(2*pi) + .5*log(det(invhess)) ...
            - DsgeLikelihood(xparam1,gend,data);
    else
      md_Laplace = .5*size(xparam1,1)*log(2*pi) + .5*log(det(invhess)) ...
          - DsgeVarLikelihood(xparam1,gend);
    end
    oo_.MarginalDensity.LaplaceApproximation = md_Laplace;    
    disp(' ')
    disp(sprintf('Log data density [Laplace approximation] is %f.',md_Laplace))
    disp(' ')
  else
    disp(' ')
    disp('RESULTS FROM MAXIMUM LIKELIHOOD')
    tstath = zeros(nx,1);
    for i = 1:nx
      tstath(i) = abs(xparam1(i))/stdh(i);
    end
    tit1 = sprintf('%10s %10s %7s %6s\n',' ', ...
                   'Estimate','s.d.','t-stat');
    if np
      ip = nvx+nvn+ncx+ncn+1;
      disp('parameters')
      disp(tit1)
      for i=1:np
        disp(sprintf('%12s %8.4f %7.4f %7.4f', ...
                     deblank(estim_params_.param_names(i,:)), ...
                     xparam1(ip),stdh(ip),tstath(ip)));
        eval(['oo_.mle_mode.parameters.' deblank(estim_params_.param_names(i,:)) ' = xparam1(ip);']);
        eval(['oo_.mle_std.parameters.' deblank(estim_params_.param_names(i,:)) ' = stdh(ip);']); 
        ip = ip+1;
      end
    end
    if nvx
      ip = 1;
      disp('standard deviation of shocks')
      disp(tit1)
      for i=1:nvx
        k = estim_params_.var_exo(i,1);
        disp(sprintf('%12s %8.4f %7.4f %7.4f', ...
                     deblank(lgx_(k,:)),xparam1(ip), ...
                     stdh(ip),tstath(ip)));
        Sigma_e_(k,k) = xparam1(ip)*xparam1(ip);
        eval(['oo_.mle_mode.shocks_std.' deblank(lgx_(k,:)) ' = xparam1(ip);']);
        eval(['oo_.mle_std.shocks_std.' deblank(lgx_(k,:)) ' = stdh(ip);']); 
        ip = ip+1;
      end
    end
    if nvn
      disp('standard deviation of measurement errors')
      disp(tit1)
      ip = nvx+1;
      for i=1:nvn
        disp(sprintf('%12s %8.4f %7.4f %7.4f', ...
                     deblank(options_.varobs(estim_params_.var_endo(i,1),: ...
                                             )), ...
                     xparam1(ip),stdh(ip),tstath(ip)))
        eval(['oo_.mle_mode.measurement_errors_std.' deblank(options_.varobs(estim_params_.var_endo(i,1),:)) ' = xparam1(ip);']);
        eval(['oo_.mle_std.measurement_errors_std.' deblank(options_.varobs(estim_params_.var_endo(i,1),:)) ' = stdh(ip);']);      
        ip = ip+1;
      end
    end
    if ncx
      disp('correlation of shocks')
      disp(tit1)
      ip = nvx+nvn+1;
      for i=1:ncx
        k1 = estim_params_.corrx(i,1);
        k2 = estim_params_.corrx(i,2);
        name = [deblank(lgx_(k1,:)) ',' deblank(lgx_(k2,:))];
        disp(sprintf('%12s %8.4f %7.4f %7.4f', name, ...
                     xparam1(ip),stdh(ip),tstath(ip)));
        Sigma_e_(k1,k2) = xparam1(ip)*sqrt(Sigma_e_(k1,k1)*Sigma_e_(k2,k2));
        Sigma_e_(k2,k1) = Sigma_e_(k1,k2);
        eval(['oo_.mle_mode.shocks_corr.' deblank(lgx_(k1,:)) '_' deblank(lgx_(k2,:)) ' = xparam1(ip);']);
        eval(['oo_.mle_std.shocks_corr.' deblank(lgx_(k1,:)) '_' deblank(lgx_(k2,:)) ' = stdh(ip);']);      
        ip = ip+1;
      end
    end
    if ncn
      disp('correlation of measurement errors')
      disp(tit1)
      ip = nvx+nvn+ncx+1;
      for i=1:ncn
        k1 = estim_params_.corrn(i,1);
        k2 = estim_params_.corrn(i,2);
        name = [deblank(lgy_(k1,:)) ',' deblank(lgy_(k2,:))];
        disp(sprintf('%12s %8.4f %7.4f %7.4f', name, ...
                     xparam1(ip),stdh(ip),tstath(ip)));
        eval(['oo_.mle_mode.measurement_error_corr.' deblank(lgy_(k1,:)) '_' deblank(lgy_(k2,:)) ' = xparam1(ip);']);
        eval(['oo_.mle_std.measurement_error_corr.' deblank(lgy_(k1,:)) '_' deblank(lgy_(k2,:)) ' = stdh(ip);']);      
        ip = ip+1;
      end
    end
  end
  
  if any(bayestopt_.pshape > 0) & options_.TeX %%%% Bayesian estimation (posterior mode) Latex output
    if np
      filename = [fname_ '_Posterior_Mode_1.TeX'];
      fidTeX = fopen(filename,'w');
      fprintf(fidTeX,'%% TeX-table generated by dynare_estimation (Dynare).\n');
      fprintf(fidTeX,'%% RESULTS FROM POSTERIOR MAXIMIZATION (parameters)\n');
      fprintf(fidTeX,['%% ' datestr(now,0)]);
      fprintf(fidTeX,' \n');
      fprintf(fidTeX,' \n');
      fprintf(fidTeX,'{\\tiny \n')
      fprintf(fidTeX,'\\begin{table}[H]\n');
      fprintf(fidTeX,'\\centering\n');
      fprintf(fidTeX,'\\begin{tabular}{l|lcccc} \n');
      fprintf(fidTeX,'\\hline\\hline \\\\ \n');
      fprintf(fidTeX,'  & Prior distribution & Prior mean  & Prior s.d. & Posterior mode & s.d. \\\\ \n');
      fprintf(fidTeX,'\\hline \\\\ \n');
      ip = nvx+nvn+ncx+ncn+1;
      for i=1:np
        fprintf(fidTeX,'$%s$ & %s & %7.3f & %6.4f & %8.4f & %7.4f \\\\ \n',...
                deblank(estim_params_.tex(i,:)),...
                deblank(pnames(bayestopt_.pshape(ip)+1,:)),...
                bayestopt_.pmean(ip),...
                bayestopt_.pstdev(ip),...
                xparam1(ip),...
                stdh(ip));
        ip = ip + 1;    
      end
      fprintf(fidTeX,'\\hline\\hline \n');
      fprintf(fidTeX,'\\end{tabular}\n ');    
      fprintf(fidTeX,'\\caption{Results from posterior parameters (parameters)}\n ');
      fprintf(fidTeX,'\\label{Table:Posterior:1}\n');
      fprintf(fidTeX,'\\end{table}\n');
      fprintf(fidTeX,'} \n')
      fprintf(fidTeX,'%% End of TeX file.\n');
      fclose(fidTeX);
    end
    if nvx
      TeXfile = [fname_ '_Posterior_Mode_2.TeX'];
      fidTeX = fopen(TeXfile,'w');
      fprintf(fidTeX,'%% TeX-table generated by dynare_estimation (Dynare).\n');
      fprintf(fidTeX,'%% RESULTS FROM POSTERIOR MAXIMIZATION (standard deviation of structural shocks)\n');
      fprintf(fidTeX,['%% ' datestr(now,0)]);
      fprintf(fidTeX,' \n');
      fprintf(fidTeX,' \n');
      fprintf(fidTeX,'{\\tiny \n');
      fprintf(fidTeX,'\\begin{table}[H]\n');
      fprintf(fidTeX,'\\centering\n');
      fprintf(fidTeX,'\\begin{tabular}{l|lcccc} \n');
      fprintf(fidTeX,'\\hline\\hline \\\\ \n');
      fprintf(fidTeX,'  & Prior distribution & Prior mean  & Prior s.d. & Posterior mode & s.d. \\\\ \n')
      fprintf(fidTeX,'\\hline \\\\ \n');
      ip = 1;
      for i=1:nvx
        k = estim_params_.var_exo(i,1);
        fprintf(fidTeX,[ '$%s$ & %4s & %7.3f & %6.4f & %8.4f & %7.4f \\\\ \n'],...
                deblank(lgx_TeX_(k,:)),...
                deblank(pnames(bayestopt_.pshape(ip)+1,:)),...
                bayestopt_.pmean(ip),...
                bayestopt_.pstdev(ip),...
                xparam1(ip),...
                stdh(ip)); 
        ip = ip+1;
      end
      fprintf(fidTeX,'\\hline\\hline \n');
      fprintf(fidTeX,'\\end{tabular}\n ');    
      fprintf(fidTeX,'\\caption{Results from posterior parameters (standard deviation of structural shocks)}\n ');
      fprintf(fidTeX,'\\label{Table:Posterior:2}\n');
      fprintf(fidTeX,'\\end{table}\n');
      fprintf(fidTeX,'} \n')
      fprintf(fidTeX,'%% End of TeX file.\n');
      fclose(fidTeX);
    end
    if nvn
      TeXfile = [fname_ '_Posterior_Mode_3.TeX'];
      fidTeX  = fopen(TeXfile,'w');
      fprintf(fidTeX,'%% TeX-table generated by dynare_estimation (Dynare).\n');
      fprintf(fidTeX,'%% RESULTS FROM POSTERIOR MAXIMIZATION (standard deviation of measurement errors)\n');
      fprintf(fidTeX,['%% ' datestr(now,0)]);
      fprintf(fidTeX,' \n');
      fprintf(fidTeX,' \n');
      fprintf(fidTeX,'\\begin{table}[H]\n');
      fprintf(fidTeX,'\\centering\n');
      fprintf(fidTeX,'\\begin{tabular}{l|lcccc} \n');
      fprintf(fidTeX,'\\hline\\hline \\\\ \n');
      fprintf(fidTeX,'  & Prior distribution & Prior mean  & Prior s.d. &  Posterior mode & s.d. \\\\ \n')
      fprintf(fidTeX,'\\hline \\\\ \n');
      ip = nvx+1;
      for i=1:nvn
        fprintf(fidTeX,'$%s$ & %4s & %7.3f & %6.4f & %8.4f & %7.4f \\\\ \n',...
                deblank(options_.varobs_TeX(estim_params_.var_endo(i,1),:)), ...
                deblank(pnames(bayestopt_.pshape(ip)+1,:)), ...        
                bayestopt_.pmean(ip), ...
                bayestopt_.pstdev(ip),...        
                xparam1(ip),...
                stdh(ip)); 
        ip = ip+1;
      end
      fprintf(fidTeX,'\\hline\\hline \n');
      fprintf(fidTeX,'\\end{tabular}\n ');    
      fprintf(fidTeX,'\\caption{Results from posterior parameters (standard deviation of measurement errors)}\n ');
      fprintf(fidTeX,'\\label{Table:Posterior:3}\n');
      fprintf(fidTeX,'\\end{table}\n');
      fprintf(fidTeX,'%% End of TeX file.\n');
      fclose(fidTeX);
    end
    if ncx
      TeXfile = [fname_ '_Posterior_Mode_4.TeX'];
      fidTeX = fopen(TeXfile,'w');
      fprintf(fidTeX,'%% TeX-table generated by dynare_estimation (Dynare).\n');
      fprintf(fidTeX,'%% RESULTS FROM POSTERIOR MAXIMIZATION (correlation of structural shocks)\n');
      fprintf(fidTeX,['%% ' datestr(now,0)]);
      fprintf(fidTeX,' \n');
      fprintf(fidTeX,' \n');
      fprintf(fidTeX,'\\begin{table}[H]\n');
      fprintf(fidTeX,'\\centering\n');
      fprintf(fidTeX,'\\begin{tabular}{l|lcccc} \n');
      fprintf(fidTeX,'\\hline\\hline \\\\ \n');
      fprintf(fidTeX,'  & Prior distribution & Prior mean  & Prior s.d. &  Posterior mode & s.d. \\\\ \n')
      fprintf(fidTeX,'\\hline \\\\ \n');
      ip = nvx+nvn+1;
      for i=1:ncx
        k1 = estim_params_.corrx(i,1);
        k2 = estim_params_.corrx(i,2);
        fprintf(fidTeX,[ '$%s$ & %s & %7.3f & %6.4f & %8.4f & %7.4f \\\\ \n'],...
                [deblank(lgx_TeX_(k1,:)) ',' deblank(lgx_TeX_(k2,:))], ...
                deblank(pnames(bayestopt_.pshape(ip)+1,:)), ...
                bayestopt_.pmean(ip), ...
                bayestopt_.pstdev(ip), ...
                xparam1(ip), ...
                stdh(ip));
        ip = ip+1;
      end
      fprintf(fidTeX,'\\hline\\hline \n');
      fprintf(fidTeX,'\\end{tabular}\n ');    
      fprintf(fidTeX,'\\caption{Results from posterior parameters (correlation of structural shocks)}\n ');
      fprintf(fidTeX,'\\label{Table:Posterior:4}\n');
      fprintf(fidTeX,'\\end{table}\n');
      fprintf(fidTeX,'%% End of TeX file.\n');
      fclose(fidTeX);
    end
    if ncn
      TeXfile = [fname_ '_Posterior_Mode_5.TeX'];
      fidTeX = fopen(TeXfile,'w');
      fprintf(fidTeX,'%% TeX-table generated by dynare_estimation (Dynare).\n');
      fprintf(fidTeX,'%% RESULTS FROM POSTERIOR MAXIMIZATION (correlation of measurement errors)\n');
      fprintf(fidTeX,['%% ' datestr(now,0)]);
      fprintf(fidTeX,' \n');
      fprintf(fidTeX,' \n');
      fprintf(fidTeX,'\\begin{table}[H]\n');
      fprintf(fidTeX,'\\centering\n');
      fprintf(fidTeX,'\\begin{tabular}{l|lcccc} \n');
      fprintf(fidTeX,'\\hline\\hline \\\\ \n');
      fprintf(fidTeX,'  & Prior distribution & Prior mean  & Prior s.d. &  Posterior mode & s.d. \\\\ \n')
      fprintf(fidTeX,'\\hline \\\\ \n');
      ip = nvx+nvn+ncx+1;
      for i=1:ncn
        k1 = estim_params_.corrn(i,1);
        k2 = estim_params_.corrn(i,2);
        fprintf(fidTeX,'$%s$ & %s & %7.3f & %6.4f & %8.4f & %7.4f \\\\ \n',...
                [deblank(lgy_(k1,:)) ',' deblank(lgy_(k2,:))], ...
                pnames(bayestopt_.pshape(ip)+1,:), ...
                bayestopt_.pmean(ip), ...
                bayestopt_.pstdev(ip), ...
                xparam1(ip), ...
                stdh(ip));
        ip = ip+1;
      end
      fprintf(fidTeX,'\\hline\\hline \n');
      fprintf(fidTeX,'\\end{tabular}\n ');    
      fprintf(fidTeX,'\\caption{Results from posterior parameters (correlation of measurement errors)}\n ');
      fprintf(fidTeX,'\\label{Table:Posterior:5}\n');
      fprintf(fidTeX,'\\end{table}\n');
      fprintf(fidTeX,'%% End of TeX file.\n');
      fclose(fidTeX);
    end 
  end    
  

  if any(bayestopt_.pshape  > 0)  & (options_.mh_replic | options_.load_mh_file)  % not ML estimation
    bounds = prior_bounds(bayestopt_);
    bayestopt_.lb = bounds(:,1);
    bayestopt_.ub = bounds(:,2);
    if any(xparam1 < bounds(:,1)) | any(xparam1 > bounds(:,2))
      find(xparam1 < bounds(:,1))
      find(xparam1 > bounds(:,2))
      error('Mode values are outside prior bounds. Reduce prior_trunc.')
    end
    metropolis(xparam1,invhess,gend,data,rawdata,bounds);
    get_posterior_parameters('posterior_mean');
    dr_ = resol(ys_,0);
  end

  % computes smoothed variables for posterior_mode or posterior_mean, if
  % metropolis). Except if options_.smoother or recursive forecast
  
  if ~options_.smoother & options_.nobs(1) == options_.nobs(end) 
    [atT,innov,measurement_error,filtered_state_vector,ys,trend_coeff] ...
	= DsgeSmoother(xparam1,gend,data);
    for i=1:endo_nbr
      eval(['oo_.SmoothedVariables.' deblank(lgy_(dr.order_var(i),:)) ' = atT(i,:)'';']);
      eval(['oo_.FilteredVariables.' deblank(lgy_(dr.order_var(i),:)) ' = filtered_state_vector(i,:)'';']);
    end
    [nbplt,nr,nc,lr,lc,nstar] = pltorg(exo_nbr);
    if options_.TeX
      fidTeX = fopen([fname_ '_SmoothedShocks.TeX'],'w');
      fprintf(fidTeX,'%% TeX eps-loader file generated by dynare_estimation.m (Dynare).\n');
      fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
      fprintf(fidTeX,' \n');
    end    
    if nbplt == 1
      hh = figure('Name','Smoothed shocks');
      NAMES = [];
      if options_.TeX, TeXNAMES = [], end
      for i=1:exo_nbr
        subplot(nr,nc,i);
        plot(1:gend,innov(i,:),'-k','linewidth',1)
        hold on
        plot([1 gend],[0 0],'-r','linewidth',.5)
        hold off
        xlim([1 gend])
        name    = deblank(lgx_(i,:));
        NAMES   = strvcat(NAMES,name);
        if ~isempty(options_.XTick)
          set(gca,'XTick',options_.XTick)
          set(gca,'XTickLabel',options_.XTickLabel)
        end
        if options_.TeX
          texname = lgx_TeX_(i,:);
          TeXNAMES   = strvcat(TeXNAMES,['$ ' deblank(texname) ' $']);
        end
        title(name,'Interpreter','none')
        eval(['oo_.SmoothedShocks.' deblank(lgx_(i,:)) ' = innov(i,:)'';']);
      end
      eval(['print -depsc2 ' fname_ '_SmoothedShocks' int2str(1)]);
      eval(['print -dpdf ' fname_ '_SmoothedShocks' int2str(1)]);
      saveas(hh,[fname_ '_SmoothedShocks' int2str(1) '.fig']);
      if options_.nograph, close(hh), end
      if options_.TeX
        fprintf(fidTeX,'\\begin{figure}[H]\n');
        for jj = 1:exo_nbr
          fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
        end    
        fprintf(fidTeX,'\\centering \n');
        fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_SmoothedShocks%s}\n',fname_,int2str(1));
        fprintf(fidTeX,'\\caption{Smoothed shocks.}');
        fprintf(fidTeX,'\\label{Fig:SmoothedShocks:%s}\n',int2str(1));
        fprintf(fidTeX,'\\end{figure}\n');
        fprintf(fidTeX,'\n');
        fprintf(fidTeX,'%% End of TeX file.\n');
        fclose(fidTeX);
      end    
    else
      for plt = 1:nbplt-1
        hh = figure('Name','Smoothed shocks');
        set(0,'CurrentFigure',hh)
        NAMES = [];
        if options_.TeX, TeXNAMES = [], end
        for i=1:nstar
          k = (plt-1)*nstar+i;
          subplot(nr,nc,i);
          plot([1 gend],[0 0],'-r','linewidth',.5)
          hold on
          plot(1:gend,innov(k,:),'-k','linewidth',1)
          hold off
          name = deblank(lgx_(k,:));
          NAMES = strvcat(NAMES,name);
          if ~isempty(options_.XTick)
            set(gca,'XTick',options_.XTick)
            set(gca,'XTickLabel',options_.XTickLabel)
          end
          xlim([1 gend])
          if options_.TeX
            texname = lgx_TeX_(k,:);
            TeXNAMES = strvcat(TeXNAMES,['$ ' deblank(texname) ' $']);
          end    
          title(name,'Interpreter','none')
          eval(['oo_.SmoothedShocks.' deblank(name) ' = innov(k,:)'';']);
        end
        eval(['print -depsc2 ' fname_ '_SmoothedShocks' int2str(plt)]);
        eval(['print -dpdf ' fname_ '_SmoothedShocks' int2str(plt)]);
        saveas(hh,[fname_ '_SmoothedShocks' int2str(plt) '.fig']);
        if options_.nograph, close(hh), end
        if options_.TeX
          fprintf(fidTeX,'\\begin{figure}[H]\n');
          for jj = 1:nstar
            fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
          end    
          fprintf(fidTeX,'\\centering \n');
          fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_SmoothedShocks%s}\n',fname_,int2str(plt));
          fprintf(fidTeX,'\\caption{Smoothed shocks.}');
          fprintf(fidTeX,'\\label{Fig:SmoothedShocks:%s}\n',int2str(plt));
          fprintf(fidTeX,'\\end{figure}\n');
          fprintf(fidTeX,'\n');
        end    
      end
      hh = figure('Name','Smoothed shocks');
      set(0,'CurrentFigure',hh)
      NAMES = [];
      if options_.TeX, TeXNAMES = [], end
      for i=1:exo_nbr-(nbplt-1)*nstar
        k = (nbplt-1)*nstar+i;
        if lr ~= 0
          subplot(lr,lc,i);
        else
          subplot(nr,nc,i);
        end    
        plot([1 gend],[0 0],'-r','linewidth',0.5)
        hold on
        plot(1:gend,innov(k,:),'-k','linewidth',1)
        hold off
        name     = deblank(lgx_(k,:));
        NAMES    = strvcat(NAMES,name);
        if ~isempty(options_.XTick)
          set(gca,'XTick',options_.XTick)
          set(gca,'XTickLabel',options_.XTickLabel)
        end
        xlim([1 gend])
        if options_.TeX
          texname  = lgx_TeX_(k,:);
          TeXNAMES = strvcat(TeXNAMES,['$ ' deblank(texname) ' $']);
        end
        title(name,'Interpreter','none')
        eval(['oo_.SmoothedShocks.' deblank(name) ' = innov(k,:)'';']);
      end
      eval(['print -depsc2 ' fname_ '_SmoothedShocks' int2str(nbplt)]);
      eval(['print -dpdf ' fname_ '_SmoothedShocks' int2str(nbplt)]);
      saveas(hh,[fname_ '_SmoothedShocks' int2str(nbplt) '.fig']);
      if options_.nograph, close(hh), end
      if options_.TeX
        fprintf(fidTeX,'\\begin{figure}[H]\n');
        for jj = 1:size(NAMES,1);
          fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
        end    
        fprintf(fidTeX,'\\centering \n');
        fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_SmoothedShocks%s}\n',fname_,int2str(nbplt));
        fprintf(fidTeX,'\\caption{Smoothed shocks.}');
        fprintf(fidTeX,'\\label{Fig:SmoothedShocks:%s}\n',int2str(nbplt));
        fprintf(fidTeX,'\\end{figure}\n');
        fprintf(fidTeX,'\n');
        fprintf(fidTeX,'%% End of TeX file.\n');
        fclose(fidTeX);
      end    
    end
    %%
    %%  Smooth observational errors...
    %%
    yf = zeros(gend,n_varobs);
    if options_.prefilter == 1
      yf = atT(bayestopt_.mf,:)+repmat(bayestopt_.mean_varobs',1,gend);
    elseif options_.loglinear == 1
      yf = atT(bayestopt_.mf,:)+repmat(log(ys(bayestopt_.mfys)),1,gend)+...
       trend_coeff*[1:gend];
    else
      yf = atT(bayestopt_.mf,:)+repmat(ys(bayestopt_.mfys),1,gend)+...
       trend_coeff*[1:gend];
    end
    if nvn
      number_of_plots_to_draw = 0;
      index = [];
      for i=1:n_varobs
    if max(abs(measurement_error(10:end))) > 0.000000001
      number_of_plots_to_draw = number_of_plots_to_draw + 1;
      index = cat(1,index,i);
    end
    eval(['oo_.SmoothedMeasurementErrors.' deblank(options_.varobs(i,:)) ...
          ' = measurement_error(i,:)'';']);
      end
      [nbplt,nr,nc,lr,lc,nstar] = pltorg(number_of_plots_to_draw);
      if options_.TeX
    fidTeX = fopen([fname_ '_SmoothedObservationErrors.TeX'],'w');
    fprintf(fidTeX,'%% TeX eps-loader file generated by dynare_estimation.m (Dynare).\n');
    fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
    fprintf(fidTeX,' \n');
      end    
      if nbplt == 1
    hh = figure('Name','Smoothed observation errors');
    set(0,'CurrentFigure',hh)
    NAMES = [];
    if options_.TeX, TeXNAMES = [], end
    for i=1:number_of_plots_to_draw
      subplot(nr,nc,i);
      plot(1:gend,measurement_error(index(i),:),'-k','linewidth',1)
      hold on
      plot([1 gend],[0 0],'-r','linewidth',.5)
      hold off
      name    = deblank(options_.varobs(index(i),:));
      NAMES   = strvcat(NAMES,name);
      if ~isempty(options_.XTick)
        set(gca,'XTick',options_.XTick)
        set(gca,'XTickLabel',options_.XTickLabel)
      end
      if options_.TeX
        texname = options_.varobs_TeX(index(i),:);
        TeXNAMES   = strvcat(TeXNAMES,['$ ' deblank(texname) ' $']);
      end
      title(name,'Interpreter','none')
    end
    eval(['print -depsc2 ' fname_ '_SmoothedObservationErrors' int2str(1)]);
    eval(['print -dpdf ' fname_ '_SmoothedObservationErrors' int2str(1)]);
    saveas(hh,[fname_ '_SmoothedObservationErrors' int2str(1) '.fig']);
    if options_.nograph, close(hh), end
    if options_.TeX
      fprintf(fidTeX,'\\begin{figure}[H]\n');
      for jj = 1:number_of_plots_to_draw
        fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
      end    
      fprintf(fidTeX,'\\centering \n');
      fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_SmoothedObservationErrors%s}\n',fname_,int2str(1));
      fprintf(fidTeX,'\\caption{Smoothed observation errors.}');
      fprintf(fidTeX,'\\label{Fig:SmoothedObservationErrors:%s',int2str(1));
      fprintf(fidTeX,'\\end{figure}\n');
      fprintf(fidTeX,'\n');
      fprintf(fidTeX,'%% End of TeX file.\n');
      fclose(fidTeX);
    end    
      else
    for plt = 1:nbplt-1
      hh = figure('Name','Smoothed observation errors');
      set(0,'CurrentFigure',hh)
      NAMES = [];
      if options_.TeX, TeXNAMES = [], end
      for i=1:nstar
        k = (plt-1)*nstar+i;
        subplot(nr,nc,i);
        plot([1 gend],[0 0],'-r','linewidth',.5)
        hold on
        plot(1:gend,measurement_error(index(k),:),'-k','linewidth',1)
        hold off
        name = deblank(options_.varobs(index(k),:));
        NAMES = strvcat(NAMES,name);
        if ~isempty(options_.XTick)
          set(gca,'XTick',options_.XTick)
          set(gca,'XTickLabel',options_.XTickLabel)
        end
        if options_.TeX
          texname = options_.varobs_TeX(k,:);
          TeXNAMES = strvcat(TeXNAMES,['$ ' deblank(texname) ' $']);
        end    
        title(name,'Interpreter','none')
      end
      eval(['print -depsc2 ' fname_ '_SmoothedObservationErrors' int2str(plt)]);
      eval(['print -dpdf ' fname_ '_SmoothedObservationErrors' int2str(plt)]);
      saveas(hh,[fname_ '_SmoothedObservationErrors' int2str(plt) '.fig']);
      if options_.nograph, close(hh), end
      if options_.TeX
        fprintf(fidTeX,'\\begin{figure}[H]\n');
        for jj = 1:nstar
          fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
        end    
        fprintf(fidTeX,'\\centering \n');
        fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_SmoothedObservationErrors%s}\n',fname_,int2str(plt));
        fprintf(fidTeX,'\\caption{Smoothed observation errors.}');
        fprintf(fidTeX,'\\label{Fig:SmoothedObservationErrors:%s}\n',int2str(plt));
        fprintf(fidTeX,'\\end{figure}\n');
        fprintf(fidTeX,'\n');
      end    
    end
    hh = figure('Name','Smoothed observation errors');
    set(0,'CurrentFigure',hh)
    NAMES = [];
    if options_.TeX, TeXNAMES = [], end
    for i=1:number_of_plots_to_draw-(nbplt-1)*nstar
      k = (nbplt-1)*nstar+i;
      if lr ~= 0
        subplot(lr,lc,i);
      else
        subplot(nr,nc,i);
      end    
      plot([1 gend],[0 0],'-r','linewidth',0.5)
      hold on
      plot(1:gend,measurement_error(index(k),:),'-k','linewidth',1)
      hold off
      name     = deblank(options_.varobs(index(k),:));
      NAMES    = strvcat(NAMES,name);
      if ~isempty(options_.XTick)
        set(gca,'XTick',options_.XTick)
        set(gca,'XTickLabel',options_.XTickLabel)
      end
      if options_.TeX
        texname  = options_.varobs_TeX(index(k),:);
        TeXNAMES = strvcat(TeXNAMES,['$ ' deblank(texname) ' $']);
      end
      title(name,'Interpreter','none');
    end
    eval(['print -depsc2 ' fname_ '_SmoothedObservationErrors' int2str(nbplt)]);
    eval(['print -dpdf ' fname_ '_SmoothedObservationErrors' int2str(nbplt)]);
    saveas(hh,[fname_ '_SmoothedObservationErrors' int2str(nbplt) '.fig']);
    if options_.nograph, close(hh), end
    if options_.TeX
      fprintf(fidTeX,'\\begin{figure}[H]\n');
      for jj = 1:size(NAMES,1);
        fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
      end    
      fprintf(fidTeX,'\\centering \n');
      fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_SmoothedObservedErrors%s}\n',fname_,int2str(nbplt));
      fprintf(fidTeX,'\\caption{Smoothed observed errors.}');
      fprintf(fidTeX,'\\label{Fig:SmoothedObservedErrors:%s}\n',int2str(nbplt));
      fprintf(fidTeX,'\\end{figure}\n');
      fprintf(fidTeX,'\n');
      fprintf(fidTeX,'%% End of TeX file.\n');
      fclose(fidTeX);
    end    
      end
    end 
    %%
    %%  Historical and smoothed variabes
    %%
    if nvn 
      [nbplt,nr,nc,lr,lc,nstar] = pltorg(n_varobs);
      if options_.TeX
    fidTeX = fopen([fname_ '_HistoricalAndSmoothedVariables.TeX'],'w');
    fprintf(fidTeX,'%% TeX eps-loader file generated by dynare_estimation.m (Dynare).\n');
    fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
    fprintf(fidTeX,' \n');
      end    
      if nbplt == 1
    hh = figure('Name','Historical and smoothed variables');
    NAMES = [];
    if options_.TeX, TeXNAMES = [], end
    for i=1:n_varobs
      subplot(nr,nc,i);
      plot(1:gend,yf(i,:),'-r','linewidth',1)
      hold on
      plot(1:gend,rawdata(options_.first_obs+(0:gend-1),i),'-k','linewidth',1)
      hold off
      name    = deblank(options_.varobs(i,:));
      NAMES   = strvcat(NAMES,name);
      if ~isempty(options_.XTick)
        set(gca,'XTick',options_.XTick)
        set(gca,'XTickLabel',options_.XTickLabel)
      end
      xlim([1 gend])
      if options_.TeX
        texname = options_.varobs_TeX(i,1);
        TeXNAMES   = strvcat(TeXNAMES,['$ ' deblank(texname) ' $']);
      end
      title(name,'Interpreter','none')
    end
    eval(['print -depsc2 ' fname_ '_HistoricalAndSmoothedVariables' int2str(1)]);
    eval(['print -dpdf ' fname_ '_HistoricalAndSmoothedVariables' int2str(1)]);
    saveas(hh,[fname_ '_HistoricalAndSmoothedVariables' int2str(1) '.fig']);
    if options_.nograph, close(hh), end
    if options_.TeX
      fprintf(fidTeX,'\\begin{figure}[H]\n');
      for jj = 1:n_varobs
        fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
      end    
      fprintf(fidTeX,'\\centering \n');
      fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_HistoricalAndSmoothedVariables%s}\n',fname_,int2str(1));
      fprintf(fidTeX,'\\caption{Historical and smoothed variables.}');
      fprintf(fidTeX,'\\label{Fig:HistoricalAndSmoothedVariables:%s}\n',int2str(1));
      fprintf(fidTeX,'\\end{figure}\n');
      fprintf(fidTeX,'\n');
      fprintf(fidTeX,'%% End of TeX file.\n');
      fclose(fidTeX);
    end    
      else
    for plt = 1:nbplt-1
      hh = figure('Name','Historical and smoothed variables');
      set(0,'CurrentFigure',hh)
      NAMES = [];
      if options_.TeX, TeXNAMES = [], end
      for i=1:nstar
        k = (plt-1)*nstar+i;
        subplot(nr,nc,i);
        plot(1:gend,yf(k,:),'-r','linewidth',1)
        hold on
        plot(1:gend,rawdata(options_.first_obs+(0:gend-1),k),'-k','linewidth',1)
        hold off
        name = deblank(options_.varobs(k,:));
        NAMES = strvcat(NAMES,name);
        if ~isempty(options_.XTick)
          set(gca,'XTick',options_.XTick)
          set(gca,'XTickLabel',options_.XTickLabel)
        end
        xlim([1 gend])
        if options_.TeX
          texname = options_.varobs_TeX(k,:);
          TeXNAMES = strvcat(TeXNAMES,['$ ' deblank(texname) ' $']);
        end    
        title(name,'Interpreter','none')
      end
      eval(['print -depsc2 ' fname_ '_HistoricalAndSmoothedVariables' int2str(plt)]);
      eval(['print -dpdf ' fname_ '_HistoricalAndSmoothedVariables' int2str(plt)]);
      saveas(hh,[fname_ '_HistoricalAndSmoothedVariables' int2str(plt) '.fig']);
      if options_.nograph, close(hh), end
      if options_.TeX
        fprintf(fidTeX,'\\begin{figure}[H]\n');
        for jj = 1:nstar
          fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
        end    
        fprintf(fidTeX,'\\centering \n');
        fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_HistoricalAndSmoothedVariables%s}\n',fname_,int2str(plt));
        fprintf(fidTeX,'\\caption{Historical and smoothed variables.}');
        fprintf(fidTeX,'\\label{Fig:HistoricalAndSmoothedVariables:%s}\n',int2str(plt));
        fprintf(fidTeX,'\\end{figure}\n');
        fprintf(fidTeX,'\n');
      end    
    end
    hh = figure('Name','Historical and smoothed variables');
    set(0,'CurrentFigure',hh)
    NAMES = [];
    if options_.TeX, TeXNAMES = [], end
    for i=1:n_varobs-(nbplt-1)*nstar
      k = (nbplt-1)*nstar+i;
      if lr ~= 0
        subplot(lr,lc,i);
      else
        subplot(nr,nc,i);
      end    
      plot(1:gend,yf(k,:),'-r','linewidth',1)
      hold on
      plot(1:gend,rawdata(options_.first_obs+(0:gend-1),k),'-k','linewidth',1)
      hold off
      name = deblank(options_.varobs(k,:));
      NAMES    = strvcat(NAMES,name);
      if ~isempty(options_.XTick)
        set(gca,'XTick',options_.XTick)
        set(gca,'XTickLabel',options_.XTickLabel)
      end
      xlim([1 gend])
      if options_.TeX
        texname  = options_.varobs_TeX(k,:);
        TeXNAMES = strvcat(TeXNAMES,['$ ' deblank(texname) ' $']);
      end
      title(name,'Interpreter','none');
    end
    eval(['print -depsc2 ' fname_ '_HistoricalAndSmoothedVariables' int2str(nbplt)]);
    eval(['print -dpdf ' fname_ '_HistoricalAndSmoothedVariables' int2str(nbplt)]);
    saveas(hh,[fname_ '_HistoricalAndSmoothedVariables' int2str(nbplt) '.fig']);
    if options_.nograph, close(hh), end
    if options_.TeX
      fprintf(fidTeX,'\\begin{figure}[H]\n');
      for jj = 1:size(NAMES,1);
        fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
      end    
      fprintf(fidTeX,'\\centering \n');
      fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_HistoricalAndSmoothedVariables%s}\n',fname_,int2str(nbplt));
      fprintf(fidTeX,'\\caption{Historical and smoothed variables.}');
      fprintf(fidTeX,'\\label{Fig:HistoricalAndSmoothedVariables:%s}\n',int2str(nbplt));
      fprintf(fidTeX,'\\end{figure}\n');
      fprintf(fidTeX,'\n');
      fprintf(fidTeX,'%% End of TeX file.\n');
      fclose(fidTeX);
    end    
      end % if nbplt == 1
    end   % if nvn
  end     % if ~((any(bayestopt_.pshape > 0) & options_.mh_replic) 
          % | (any(bayestopt_.pshape > 0) & options_.load_mh_file)) | ~options_.smoother  

  if options_.nobs(1) < options_.nobs(end)  % recursive estimation
    if isempty(varlist)
      varlist = lgy_;
    end
    for i = 1:size(varlist,1)
      eval(['oo_.RecursiveForecast.Mean.' varlist(i,:) '(' int2str(iter1) ',:) =' ...
            ' oo_.Forecast.Mean.' varlist(i,:) ';']); 
      eval(['oo_.RecursiveForecast.Median.' varlist(i,:) '(' int2str(iter1) ',:) =' ...
            ' oo_.Forecast.Median.' varlist(i,:) ';']); 
      eval(['oo_.RecursiveForecast.Std.' varlist(i,:) '(' int2str(iter1) ',:) =' ...
            ' oo_.Forecast.Std.' varlist(i,:) ';']); 
      eval(['oo_.RecursiveForecast.HPDinf.' varlist(i,:) '(' int2str(iter1) ',:) =' ...
            ' oo_.Forecast.HPDinf.' varlist(i,:) ';']); 
      eval(['oo_.RecursiveForecast.HPDsup.' varlist(i,:) '(' int2str(iter1) ',:) =' ...
            ' oo_.Forecast.HPDsup.' varlist(i,:) ';']); 
      eval(['oo_.RecursiveForecast.HPDTotalinf.' varlist(i,:) '(' int2str(iter1) ',:) =' ...
            ' oo_.Forecast.HPDTotalinf.' varlist(i,:) ';']); 
      eval(['oo_.RecursiveForecast.HPDTotalsup.' varlist(i,:) '(' int2str(iter1) ',:) =' ...
            ' oo_.Forecast.HPDTotalsup.' varlist(i,:) ';']); 
    end
    iter1 = iter1+1;    
  end
end        % end for gend=options_.nobs(1):options_.nobs(end)

if options_.nobs(1) < options_.nobs(end) & options_.forecast > 1
  if isempty(varlist)
    varlist = lgy_;
    nvar    = size(lgy_,1);
    SelecVariables = transpose(1:nvar);
  else
    nvar = size(varlist,1);
    SelecVariables = [];
    for i=1:nvar
      if ~isempty(strmatch(varlist(i,:),lgy_,'exact'))
    SelecVariables = [SelecVariables;strmatch(varlist(i,:),lgy_, ...
                          'exact')];
      else
    error(['Estimation: ' varlist(i,:) ' isn''t an endogenous' ...
           'variable'])
      end   
    end
    IdObs    = zeros(n_varobs,1);
    for j=1:n_varobs
      for i=1:nvar
    iobs = strmatch(options_.varobs(j,:),varlist,'exact');
      end
      if ~isempty(iobs)
    IdObs(j,1) = iobs;
      end    
    end 
  end
  k = 3+min(options_.nobs(end)-options_.nobs(1)+options_.forecast, ...
      size(rawdata,1)-options_.nobs(1));
  data2 = rawdata(end-k+1:end,:);
  [nbplt,nr,nc,lr,lc,nstar] = pltorg(nvar);
  for i = 1:size(varlist,1)
    if mod(i,nstar) == 1
      hfig = figure('Name','Out of sample forecasts');
    end
    subplot(nr,nc,i)
    hold on
    if any(i==IdObs)
      k2 = find(i==IdObs);
      if options_.loglinear == 1
    plot(1:k,exp(data2(end-k+1:end,k2))','-k','linewidth',2);
      else
    plot(1:k,data2(k2,end-k+1:end,k2)','-k','linewidth',2);
      end
      offsetx = 3;
    else
      offsetx = 0;
    end      
    for j=1:options_.nobs(end)-options_.nobs(1)+1
      y = eval(['oo_.RecursiveForecast.Mean.' varlist(i,:) '(j,:)']);
      y1 = eval(['oo_.RecursiveForecast.HPDinf.' varlist(i,:) '(j,:)']);
      y2 = eval(['oo_.RecursiveForecast.HPDsup.' varlist(i,:) '(j,:)']);
      y3 = eval(['oo_.RecursiveForecast.HPDTotalinf.' varlist(i,:) '(j,:)']);
      y4 = eval(['oo_.RecursiveForecast.HPDTotalsup.' varlist(i,:) '(j,:)']);
      plot(offsetx+[j:j+options_.forecast-1],y,'-b','linewidth',2)
      plot(offsetx+[j:j+options_.forecast-1],y1,'--g', ...
       'linewidth',1.5)
      plot(offsetx+[j:j+options_.forecast-1],y2,'--g', ...
       'linewidth',1.5)
      plot(offsetx+[j:j+options_.forecast-1],y3,'--r', ...
       'linewidth',1.5)
      plot(offsetx+[j:j+options_.forecast-1],y4,'--r','linewidth',1.5)
    end
%    set(gca,'XTick',offsetx+[1 10 20 30 40 50 60 70 80 90]);
%    set(gca,'XTickLabel',{'1';'10';'20';'30';'40';'50';'60';'70';'80';'90'});
    %   xlim([1 options_.forecast+10]);
    if any(k==IdObs)
      plot([offsetx+1 offsetx+1],ylim,'-c')
    end
    box on
    title(deblank(varlist(i,:)),'Interpreter','none')
    hold off
  end
end







      
% SA 07-31-2004     * Added TeX output.
%                   * Prior plots are done by calling plot_priors.m.
%                   * All the computations related to the metropolis-hastings are made
%                   in a new version of metropolis.m.
%                   * Corrected a bug related to prior's bounds.
%                   * ...
%                   * If you do not want to see all the figures generated by dynare, you can use the option
%                   nograph. The figures will be done and saved in formats eps, pdf and fig (so that you
%                   should be able to modify the plots within matlab) but each figure will be erased from the
%                   workspace when completed.
% SA 08-04-2004     Corrected a bug related to the display of the Smooth shocks and variables plots,
%                   for ML and posterior mode estimation. 
% SA 09-03-2004     Compilation of TeX appendix moved to dynare.m.                                                                                                                                                                         h
