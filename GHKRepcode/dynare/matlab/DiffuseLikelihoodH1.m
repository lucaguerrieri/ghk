function LIK = DiffuseLikelihoodH1(T,R,Q,H,Pinf,Pstar,Y,trend,start)
% stephane.adjemian@cepremap.cnrs.fr [07-19-2004]
%
%   See "Filtering and Smoothing of State Vector for Diffuse State Space
%   Models", S.J. Koopman and J. Durbin (2003, in Journal of Time Series 
%   Analysis, vol. 24(1), pp. 85-98).  
%
%   THE PROBLEM:
%
%   y_t =   Z_t * \alpha_t + \varepsilon_t
%   \alpha_{t+1} = T_t  * \alpha_t + R_t * \eta_t
%
%   with:
%
%   \alpha_1 = a + A*\delta + R_0*\eta_0
%
%   m*q matrix A and m*(m-q) matrix R_0 are selection matrices (their
%   columns constitue all the columns of the m*m identity matrix) so that 
%
%       A'*R_0 = 0 and A'*\alpha_1 = \delta
%
%   We assume that the vector \delta is distributed as a N(0,\kappa*I_q)
%   for a given  \kappa > 0. So that the expectation of \alpha_1 is a and
%   its variance is P, with
%
%       P = \kappa*P_{\infty} + P_{\star}
%
%           P_{\infty} = A*A'
%           P_{\star}  = R_0*Q_0*R_0'
%
%   P_{\infty} is a m*m diagonal matrix with q ones and m-q zeros. 
%
%
%   and where:
%
%   y_t             is a pp*1 vector
%   \alpha_t        is a mm*1 vector
%   \varepsilon_t   is a pp*1 multivariate random variable (iid N(0,H_t))
%   \eta_t          is a rr*1 multivariate random variable (iid N(0,Q_t))
%   a_1             is a mm*1 vector
%
%   Z_t     is a pp*mm matrix
%   T_t     is a mm*mm matrix
%   H_t     is a pp*pp matrix
%   R_t     is a mm*rr matrix
%   Q_t     is a rr*rr matrix
%   P_1     is a mm*mm matrix
%
%
%   FILTERING EQUATIONS:
%
%   v_t = y_t - Z_t* a_t
%   F_t = Z_t * P_t * Z_t' + H_t
%   K_t = T_t * P_t * Z_t' * F_t^{-1}
%   L_t = T_t - K_t * Z_t
%   a_{t+1} = T_t * a_t + K_t * v_t
%   P_{t+1} = T_t * P_t * L_t' + R_t*Q_t*R_t'
%
%
%   DIFFUSE FILTERING EQUATIONS:
%
%   a_{t+1} = T_t*a_t + K_{\infty,t}v_t
%   P_{\infty,t+1} = T_t*P_{\infty,t}*L_{\infty,t}'
%   P_{\ast,t+1}  = T_t*P_{\ast,t}*L_{\infty,t}' - K_{\infty,t}*F_{\infty,t}*K_{\ast,t}' + R_t*Q_t*R_t'
%   K_{\infty,t}   = T_t*P_{\infty,t}*Z_t'*F_{\infty,t}^{-1}
%   v_t = y_t - Z_t*a_t
%   L_{\infty,t} = T_t - K_{\infty,t}*Z_t
%   F_{\infty,t} = Z_t*P_{\infty,t}*Z_t'
%   K_{\ast,t}  = (T_t*P_{\ast,t}*Z_t' + K_{\infty,t}*F_{\ast,t})*F_{\infty,t}^{-1}
%   F_{\ast,t}  = Z_t*P_{\ast,t}*Z_t' + H_t
%
%	Matrix Finf is assumed to be non singular. If this is not the case we have
%   to switch to another algorithm (NewAlg=3).
%
%	start = options_.presample
  global bayestopt_ options_
  
  mf = bayestopt_.mf;
  smpl = size(Y,2);
  mm   = size(T,2);
  pp   = size(Y,1);
  a    = zeros(mm,1);
  dF = 1;
  QQ   = R*Q*transpose(R);
  t    = 0;
  lik  = zeros(smpl+1,1);
  LIK  = Inf;
  lik(smpl+1) = smpl*pp*log(2*pi);
  notsteady   = 1;
  crit        = options_.kalman_tol;
  reste       = 0;
  while rank(Pinf,crit) & t < smpl
    t     = t+1;
    v  	  = Y(:,t)-a(mf)-trend(:,t);
    Finf  = Pinf(mf,mf);
    if rcond(Finf) < crit 
      if ~all(abs(Finf(:))<crit)
	return
      else
	iFstar	= inv(Pstar(mf,mf)+H);
	dFstar	= det(Pstar(mf,mf)+H);
	Kstar	= Pstar(:,mf)*iFstar;
	lik(t)	= log(dFstar) + transpose(v)*iFstar*v;
	Pinf	= T*Pinf*transpose(T);
	Pstar	= T*(Pstar-Pstar(:,mf)*transpose(Kstar))*transpose(T)+QQ;
	a		= T*(a+Kstar*v);
      end
    else
      lik(t)	= log(det(Finf));
      iFinf	= inv(Finf);
      Kinf	= Pinf(:,mf)*iFinf;					%%	premultiplication by the transition matrix T is removed (stephane) 
      Fstar	= Pstar(mf,mf)+H;
      Kstar	= (Pstar(:,mf)-Kinf*Fstar)*iFinf; 	%%	premultiplication by the transition matrix T is removed (stephane)
      Pstar	= T*(Pstar-Pstar(:,mf)*transpose(Kinf)-Pinf(:,mf)*transpose(Kstar))*transpose(T)+QQ;
      Pinf	= T*(Pinf-Pinf(:,mf)*transpose(Kinf))*transpose(T);
      a		= T*(a+Kinf*v);					
    end  
  end
  if t == smpl                                                           
    error(['There isn''t enough information to estimate the initial' ... 
	   ' conditions of the nonstationary variables']);                   
  end                                                                    
  F_singular = 1;
  while notsteady & t < smpl
    t  = t+1;
    v  = Y(:,t)-a(mf)-trend(:,t);
    F  = Pstar(mf,mf)+H;
    oldPstar  = Pstar;
    dF = det(F);
    if rcond(F) < crit 
      if ~all(abs(F(:))<crit)
	return
      else
	a         = T*a;
	Pstar     = T*Pstar*transpose(T)+QQ;
      end
    else  
      F_singular = 0;
      iF		  = inv(F);
      lik(t)    = log(dF)+transpose(v)*iF*v;
      K         = Pstar(:,mf)*iF; %% premultiplication by the transition matrix T is removed (stephane)
      a         = T*(a+K*v);		%% --> factorization of the transition matrix...
      Pstar     = T*(Pstar-Pstar(:,mf)*iF*Pstar(mf,:))*transpose(T)+QQ;	%% ... idem (stephane)
    end
    notsteady = ~(max(max(abs(Pstar-oldPstar)))<crit);
  end
  if F_singular == 1
    error(['The variance of the forecast error remains singular until the' ...
	   'end of the sample'])
  end
  reste = smpl-t;
  while t < smpl
    t = t+1;
    v = Y(:,t)-a(mf)-trend(:,t);
    a = T*(a+K*v);
    lik(t) = transpose(v)*iF*v;
  end
  lik(t) = lik(t) + reste*log(dF);
  LIK    = .5*(sum(lik(start:end))-(start-1)*lik(smpl+1)/smpl);% Minus the
							       % log-likelihood.
							       