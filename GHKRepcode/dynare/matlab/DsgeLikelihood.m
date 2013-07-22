function [fval,cost_flag,ys,trend_coeff,info] = DsgeLikelihood(xparam1,gend,data)
% stephane.adjemian@cepremap.cnrs.fr [09-07-2004]
%
% Adapted from mj_optmumlik.m
global bayestopt_ exo_nbr dr_ estim_params_ Sigma_e_ options_ xparam1_test
  global trend_coeff_

  fval = [];
  cost_flag = [];
  ys = [];
  trend_coeff = [];

  xparam1_test = xparam1;
  cost_flag   = 1;
  nobs 	= size(options_.varobs,1);

  if options_.mode_compute ~= 1 & any(xparam1 < bayestopt_.lb)
    k = find(xparam1 < bayestopt_.lb);
    fval = bayestopt_.penalty+sum((bayestopt_.lb(k)-xparam1(k)).^2);
    cost_flag = 0;
    return;
  end
  if options_.mode_compute ~= 1 & any(xparam1 > bayestopt_.ub)
    k = find(xparam1 > bayestopt_.ub);
    fval = bayestopt_.penalty+sum((xparam1(k)-bayestopt_.ub(k)).^2);
    cost_flag = 0;
    return;
  end
  Q = Sigma_e_;
  for i=1:estim_params_.nvx
    k =estim_params_.var_exo(i,1);
    Q(k,k) = xparam1(i)*xparam1(i);
  end
  offset = estim_params_.nvx;
  if estim_params_.nvn
    H = zeros(nobs,nobs);
    for i=1:estim_params_.nvn
      k = estim_params_.var_endo(i,1);
      H(k,k) = xparam1(i+offset)*xparam1(i+offset);
    end
    offset = offset+estim_params_.nvn;
  end	
  if estim_params_.ncx
    for i=1:estim_params_.ncx
      k1 =estim_params_.corrx(i,1);
      k2 =estim_params_.corrx(i,2);
      Q(k1,k2) = xparam1(i+offset)*sqrt(Q(k1,k1)*Q(k2,k2));
      Q(k2,k1) = Q(k1,k2);
    end
    [CholQ,testQ] = chol(Q);
    if testQ 	%% The variance-covariance matrix of the structural innovations is not definite positive.
		%% We have to compute the eigenvalues of this matrix in order to build the penalty.
		a = eig(Q);
		k = a<0;
		if k > 0
		  fval = bayestopt_.penalty+sum(-a(k));
		  cost_flag = 0;
		  return
		end
    end
    offset = offset+estim_params_.ncx;
  end
  if estim_params_.nvn & estim_params_.ncn 
    for i=1:estim_params_.ncn
      k1 = options_.lgyidx2varobs(estim_params_.corrn(i,1));
      k2 = options_.lgyidx2varobs(estim_params_.corrn(i,2));
      H(k1,k2) = xparam1(i+offset)*sqrt(H(k1,k1)*H(k2,k2));
      H(k2,k1) = H(k1,k2);
    end
    [CholH,testH] = chol(H);
    if testH
      a = eig(H);
      k = a<0;
      if k > 0
	fval = bayestopt_.penalty+sum(-a(k));
	cost_flag = 0;
	return
      end
    end
    offset = offset+estim_params_.ncn;
  end
  for i=1:estim_params_.np
    assignin('base',deblank(estim_params_.param_names(i,:)),xparam1(i+offset));
  end
  Sigma_e_ = Q;
  %------------------------------------------------------------------------------
  % 2. call model setup & reduction program
  %------------------------------------------------------------------------------
  [T,R,SteadyState,info] = dynare_resolve;
  if info(1) == 1 | info(1) == 2 | info(1) == 5
    fval = bayestopt_.penalty+1;
    cost_flag = 0;
    return
  elseif info(1) == 3 | info(1) == 4 | info(1) == 20
    fval = bayestopt_.penalty+info(2)^2;
    cost_flag = 0;
    return
  end
  if options_.loglinear == 1
    constant = log(SteadyState(bayestopt_.mfys));
  else
    constant = SteadyState(bayestopt_.mfys);
  end
  if bayestopt_.with_trend == 1
    trend_coeff = zeros(nobs,1);
    for i=1:nobs
      trend_coeff(i) = evalin('base',bayestopt_.trend_coeff{i});
    end
    trend = constant*ones(1,gend)+trend_coeff*(1:gend);
  else
    trend = constant*ones(1,gend);
  end
  start = options_.presample+1;
  np    = size(T,1);
  mf    = bayestopt_.mf;
  %------------------------------------------------------------------------------
  % 3. Initial condition of the Kalman filter
  %------------------------------------------------------------------------------
  if options_.lik_init == 1		% Kalman filter
    Pstar = lyapunov_symm(T,R*Q*transpose(R));
    Pinf	= [];
  elseif options_.lik_init == 2	% Old Diffuse Kalman filter
    Pstar = 10*eye(np);
    Pinf	= [];
  elseif options_.lik_init == 3	% Diffuse Kalman filter
    Pstar = zeros(np,np);
    ivs = bayestopt_.i_T_var_stable;
    Pstar(ivs,ivs) = lyapunov_symm(T(ivs,ivs),R(ivs,:)*Q* ...
				   transpose(R(ivs,:)));
    Pinf  = bayestopt_.Pinf;
  end
  %------------------------------------------------------------------------------
  % 4. Likelihood evaluation
  %------------------------------------------------------------------------------
  if estim_params_.nvn
    if options_.kalman_algo == 1
      LIK = DiffuseLikelihoodH1(T,R,Q,H,Pinf,Pstar,data,trend,start);
      if isinf(LIK) & ~estim_params_.ncn %% The univariate approach considered here doesn't 
					 %%	apply when H has some off-diagonal elements.
					 LIK = DiffuseLikelihoodH3(T,R,Q,H,Pinf,Pstar,data,trend,start);
      elseif isinf(LIK) & estim_params_.ncn
	LIK = DiffuseLikelihoodH3corr(T,R,Q,H,Pinf,Pstar,data,trend,start);
      end
    elseif options_.kalman_algo == 3
      if ~estim_params_.ncn %% The univariate approach considered here doesn't 
			    %%	apply when H has some off-diagonal elements.
			    LIK = DiffuseLikelihoodH3(T,R,Q,H,Pinf,Pstar,data,trend,start);
      else
	LIK = DiffuseLikelihoodH3corr(T,R,Q,H,Pinf,Pstar,data,trend,start);
      end	
    end	  
  else
    if options_.kalman_algo == 1
      LIK = DiffuseLikelihood1(T,R,Q,Pinf,Pstar,data,trend,start);
      if isinf(LIK)
	LIK = DiffuseLikelihood3(T,R,Q,Pinf,Pstar,data,trend,start);
      end
    elseif options_.kalman_algo == 3
      LIK = DiffuseLikelihood3(T,R,Q,Pinf,Pstar,data,trend,start);
    end 	
  end
  if imag(LIK) ~= 0
    likelihood = bayestopt_.penalty;
  else
    likelihood = LIK;
  end
  % ------------------------------------------------------------------------------
  % Adds prior if necessary
  % ------------------------------------------------------------------------------
  lnprior = priordens(xparam1,bayestopt_.pshape,bayestopt_.p1,bayestopt_.p2,bayestopt_.p3,bayestopt_.p4);
  fval    = (likelihood-lnprior);
  
  