function LIK = DiffuseLikelihoodH3(T,R,Q,H,Pinf,Pstar,Y,trend,start)
% changes by M. Ratto
% introduced new global variable id_ for termination of DKF
% introduced a persistent fmax, in order to keep track the max order of
% magnitude of the 'zero' values in Pinf at DKF termination
% new icc counter for Finf steps in DKF
% new termination for DKF
% likelihood terms for Fstar must be cumulated in DKF also when Pinf is non
% zero. this bug is fixed.
%
% stephane.adjemian@cepremap.cnrs.fr [07-19-2004]
% 
%   See "Filtering and Smoothing of State Vector for Diffuse State Space
%   Models", S.J. Koopman and J. Durbin (2003, in Journal of Time Series 
%   Analysis, vol. 24(1), pp. 85-98).  
%
%	Case where F_{\infty,t} is singular ==> Univariate treatment of multivariate
%	time series.
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
%   a_{t+1} = T_t*a_t + K_{\ast,t}v_t
%   P_{\infty,t+1} = T_t*P_{\infty,t}*T_t'
%   P_{\ast,t+1}  = T_t*P_{\ast,t}*L_{\ast,t}' + R_t*Q_t*R_t'
%   K_{\ast,t}   = T_t*P_{\ast,t}*Z_t'*F_{\ast,t}^{-1}
%   v_t = y_t - Z_t*a_t
%   L_{\ast,t} = T_t - K_{\ast,t}*Z_t
%   F_{\ast,t}  = Z_t*P_{\ast,t}*Z_t' + H_t
global bayestopt_ options_
  
mf = bayestopt_.mf;
pp     = size(Y,1);
mm     = size(T,1);
smpl   = size(Y,2);
a      = zeros(mm,1);
QQ     = R*Q*transpose(R);
t      = 0;
lik    = zeros(smpl+1,1);
lik(smpl+1) = smpl*pp*log(2*pi); %% the constant of minus two times the log-likelihood 
notsteady 	= 1;
crit      	= options_.kalman_tol;
crit1      	= 1.e-6;
newRank	  	= rank(Pinf,crit1);
icc = 0;
while newRank & t < smpl %% Matrix Finf is assumed to be zero
  t = t+1;
  for i=1:pp
    v(i) 	= Y(i,t)-a(mf(i))-trend(i,t);
    Fstar 	= Pstar(mf(i),mf(i))+H(i,i);
    Finf	= Pinf(mf(i),mf(i));
    Kstar 	= Pstar(:,mf(i));
    if Finf > crit & newRank
      icc = icc + 1;
      Kinf	= Pinf(:,mf(i));
      a		= a + Kinf*v(i)/Finf;
      Pstar	= Pstar + Kinf*transpose(Kinf)*Fstar/(Finf*Finf) - ...
	  (Kstar*transpose(Kinf)+Kinf*transpose(Kstar))/Finf;
      Pinf	= Pinf - Kinf*transpose(Kinf)/Finf;
      lik(t) 	= lik(t) + log(Finf);
      % start new termination criterion for DKF
      if ~isempty(options_.diffuse_d),  
	newRank = (icc<options_.diffuse_d);  
	%if newRank & any(diag(Pinf(mf,mf))>crit)==0; %  M. Ratto this line is BUGGY
	if newRank & (any(diag(Pinf(mf,mf))>crit)==0 & rank(Pinf,crit1)==0); 
	  options_.diffuse_d = icc;
	  newRank=0;
	  disp('WARNING: Change in OPTIONS_.DIFFUSE_D in univariate DKF')
	  disp(['new OPTIONS_.DIFFUSE_D = ',int2str(icc)])
	  disp('You may have to reset the optimisation')
	end
      else
	%newRank = any(diag(Pinf(mf,mf))>crit);     % M. Ratto this line is BUGGY 
	newRank = (any(diag(Pinf(mf,mf))>crit) | rank(Pinf,crit1));                 
	if newRank==0, 
	  P0=	T*Pinf*transpose(T);
	  %newRank = any(diag(P0(mf,mf))>crit);   % M. Ratto this line is BUGGY
	  newRank = (any(diag(Pinf(mf,mf))>crit) | rank(P0,crit1));   
	  if newRank==0, 
	    options_.diffuse_d = icc;
	  end
	end                    
      end,
      % end new termination and checks for DKF and fmax
    elseif Finf > crit 
      %% Note that : (1) rank(Pinf)=0 implies that Finf = 0, (2) outside this loop (when for some i and t the condition
      %% rank(Pinf)=0 is satisfied we have P = Pstar and F = Fstar and (3) Finf = 0 does not imply that
      %% rank(Pinf)=0. [stéphane,11-03-2004].	  
      %if rank(Pinf) == 0
      % the likelihood terms should alwasy be cumulated, not only
      % when Pinf=0, otherwise the lik would depend on the ordering
      % of observed variables
      lik(t)	= lik(t) + log(Fstar) + v(i)*v(i)/Fstar;
      %end
      a 	= a + Kstar*v(i)/Fstar;
      Pstar	= Pstar - Kstar*transpose(Kstar)/Fstar;					
    else
      % disp(['zero F term in DKF for observed ',int2str(i),' ',num2str(Fi)])
    end
  end
  if newRank
    oldRank = rank(Pinf,crit1);
  else
    oldRank = 0;
  end
  a 		= T*a;
  Pstar 	= T*Pstar*transpose(T)+QQ;
  Pinf	= T*Pinf*transpose(T);
  if newRank
    newRank = rank(Pinf,crit1);
  end
  if oldRank ~= newRank
    disp('DiffuseLiklihoodH3 :: T does influence the rank of Pinf!')	
  end		 		
end
if t == smpl                                                           
  error(['There isn''t enough information to estimate the initial' ... 
	 ' conditions of the nonstationary variables']);                   
end                                                                    
while notsteady & t < smpl
  t = t+1;
  for i=1:pp
    v(i) = Y(i,t) - a(mf(i)) - trend(i,t);
    Fi   = Pstar(mf(i),mf(i))+H(i,i);
    if Fi > crit
      Ki	= Pstar(:,mf(i));
      a		= a + Ki*v(i)/Fi;
      Pstar 	= Pstar - Ki*transpose(Ki)/Fi;
      lik(t) 	= lik(t) + log(Fi) + v(i)*v(i)/Fi;
    end
  end	
  oldP 		= Pstar;
  a 		= T*a;
  Pstar 	= T*Pstar*transpose(T) + QQ;
  notsteady 	= ~(max(max(abs(Pstar-oldP)))<crit);
end
while t < smpl
  t = t+1;
  for i=1:pp
    v(i) = Y(i,t) - a(mf(i)) - trend(i,t);
    Fi   = Pstar(mf(i),mf(i))+H(i,i);
    if Fi > crit
      Ki 		= Pstar(:,mf(i));
      a 		= a + Ki*v(i)/Fi;
      Pstar 	= Pstar - Ki*transpose(Ki)/Fi;
      lik(t) 	= lik(t) + log(Fi) + v(i)*v(i)/Fi;
    end
  end	
  a = T*a;
end
LIK = .5*(sum(lik(start:end))-(start-1)*lik(smpl+1)/smpl);

