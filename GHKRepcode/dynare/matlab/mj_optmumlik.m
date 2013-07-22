% Program calculating the posterior density 
% 1. define xparam
% 2. call model setup & reduction program
% 3. prepare state space variables and kalman-filter setup
% 4. evaluate likelihood with kalman filter
% 5. evaluate prior
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
function [fval,cost_flag,atT,innov,ys,trend_coeff] = ...
      mj_optmumlik(xparam1,gend,rawdata,algo);

% algo = 1: computes filter + likelihood
% alog = 2: computes filter + likelihood + smoother

global bayestopt_ exo_nbr dr_ estim_params_ Sigma_e_ options_ xparam_test ...
    dr1_test_ trend_coeff_

xparam_test = xparam1;
cost_flag = 1;
if options_.mode_compute ~= 1 & any(xparam1 < bayestopt_.lb)
  k = find(xparam1 < bayestopt_.lb);
  fval = bayestopt_.penalty+sum(bayestopt_.lb(k)-xparam1(k));
  cost_flag = 0;
  return;
end
if options_.mode_compute ~= 1 & any(xparam1 > bayestopt_.ub)
  k = find(xparam1 > bayestopt_.ub);
  fval = bayestopt_.penalty+sum(xparam1(k)-bayestopt_.ub(k));
  cost_flag = 0;
  return;
end

nobs = size(options_.varobs,1);

q = Sigma_e_;
for i=1:estim_params_.nvx
  k =estim_params_.var_exo(i,1);
  q(k,k) = xparam1(i)*xparam1(i);
end

offset = estim_params_.nvx;
h = zeros(nobs,nobs);
for i=1:estim_params_.nvn
  k =estim_params_.var_endo(i,1);
  h(k,k) = xparam1(i+offset)*xparam1(i+offset);
end

offset = offset+estim_params_.nvn;
for i=1:estim_params_.ncx
  k1 =estim_params_.corrx(i,1);
  k2 =estim_params_.corrx(i,2);
  q(k1,k2) = xparam1(i+offset)*sqrt(q(k1,k1)*q(k2,k2));
  q(k2,k1) = q(k1,k2);
end

offset = offset+estim_params_.ncx;
for i=1:estim_params_.ncn
  k1 =estim_params_.corrn(i,1);
  k2 =estim_params_.corrn(i,2);
  h(k1,k2) = xparam1(i+offset)*sqrt(h(k1,k1)*h(k2,k2));
  h(k2,k1) = h(k1,k2);
end

offset = offset+estim_params_.ncn;
for i=1:estim_params_.np
   assignin('base',deblank(estim_params_.param_names(i,:)),xparam1(i+offset));
end

%------------------------------------------------------------------------------
% 2. call model setup & reduction program
%------------------------------------------------------------------------------

[A,B,ys] = dynare_resolve;

if dr1_test_(1) == 1
    fval = bayestopt_.penalty*exp(dr1_test_(2));
    cost_flag = 0;
    return;
elseif dr1_test_(1) == 2;
    fval = bayestopt_.penalty*exp(dr1_test_(2));
    cost_flag = 0;
    return;
elseif dr1_test_(1) == 3;
    fval = bayestopt_.penalty*exp(dr1_test_(2));
    cost_flag = 0;
    return;    
end

if options_.loglinear == 1
  ys1 = log(ys(bayestopt_.mfys));
else
  ys1 = ys(bayestopt_.mfys);
end

aa = A;

np = size(A,1);
mf = eye(np);
mf = bayestopt_.mf;

% Set initial values                                             @


at = zeros(np,gend+1);            
gconst = log(2*pi);
lik = zeros(gend,1);

if options_.lik_init == 1
  p0 = lyapunov_symm(aa,B*q*B');
elseif options_.lik_init == 2
  p0=eye(np)*10.0; 
end
pt = p0;        
BqB = B*q*B';

trend_coeff = zeros(nobs,1);
if bayestopt_.with_trend == 1
  nx1 = estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+ ...
	estim_params_.ncn;
  for i=1:nobs
    trend_coeff(i) = eval(bayestopt_.trend_coeff{i});
  end
end

not_steady = 1;
ldetf_old = NaN;
warning_state = warning;
warning off;
if algo == 1
  for t = 1:gend
    if not_steady
      ptt1 = aa*pt*aa'+BqB;               
      f = ptt1(mf,mf)+h;    
      ldetf = log(det(f));
      finv = inv(f);                    
      if any(isinf(finv)) | ~isreal(ldetf)
	disp('singularity in Kalman filter');
        fval = bayestopt_.penalty;
	warning(warning_state);
	cost_flag = 0;
	dr1_test_(1) = 4;
	return
      end
      pt = ptt1-ptt1(:,mf)*finv*ptt1(mf,:);
      if abs(ldetf-ldetf_old) < 1e-12
	not_steady = 0;
      end
      ldetf_old = ldetf;
    end
    att1 = aa*at(:,t);                
    v  = rawdata(t,:)'-att1(mf,:)-ys1;
    if bayestopt_.with_trend == 1
      v = v - trend_coeff*t;
    end
    at(:,t+1) = att1+ptt1(:,mf)*finv*v;   
    if t > options_.presample
      lik(t,1) = ldetf+v'*finv*v;
    end
  end 
elseif algo == 2
  PPt = zeros(np^2,gend);
  for t = 1:gend
    if not_steady
      ptt1 = aa*pt*aa'+BqB;               
      f = ptt1(mf,mf)+h;    
      ldetf = log(det(f));
      finv = inv(f);                    
      if any(isinf(finv)) | ~isreal(ldetf)
%        disp('singularity in Kalman filter');
        fval = bayestopt_.penalty;
	warning(warning_state);
	cost_flag = 0;
	dr1_test_(1) = 4;
	return
      end
      pt = ptt1-ptt1(:,mf)*finv*ptt1(mf,:);
      if abs(ldetf-ldetf_old) < 1e-12
	not_steady = 0;
      end
      ldetf_old = ldetf;
    end
    att1 = aa*at(:,t);                
    v  = rawdata(t,:)' - att1(mf,:) - ys1;
    if bayestopt_.with_trend == 1
      v = v - trend_coeff*t;
    end
    at(:,t+1) = att1+ptt1(:,mf)*finv*v;   
    PPt(:,t) = pt(:);
    if t > options_.presample
      lik(t,1) = ldetf+v'*finv*v;
    end
  end 
  atT =zeros(np,gend+1);
  atT(:,gend+1) = at(:,gend+1);
  innov = zeros(exo_nbr,gend);
  for t = gend:-1:1
    pt = reshape(PPt(:,t),np,np);
    ptt1 = aa*pt*aa'+BqB;
    Pstar = pt*aa'*pinv(ptt1);
    atT(:,t) = at(:,t)+Pstar*(atT(:,t+1)-aa*at(:,t));
    shocks1 = atT(:,t+1)-aa*atT(:,t);
    innov(:,t) = B\shocks1;
  end
end

warning(warning_state);

likelihood = 0.5*sum(nobs*gconst+lik(options_.presample+1:end));

if imag(likelihood) ~= 0
  
   likelihood = 10000000;
end

% ------------------------------------------------------------------------------
% PRIOR SPECIFICATION
% ------------------------------------------------------------------------------

lnprior = priordens(xparam1, bayestopt_.pshape, bayestopt_.p1, bayestopt_.p2, bayestopt_.p3, bayestopt_.p4 );

fval = (likelihood-lnprior);

% 11/18/03 MJ changed input parameters for priordens()