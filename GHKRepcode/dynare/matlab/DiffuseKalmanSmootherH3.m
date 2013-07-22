function [alphahat,epsilonhat,etahat,a1, aK] = DiffuseKalmanSmootherH3(T,R,Q,H,Pinf1,Pstar1,Y,trend,pp,mm,smpl,mf)
% Modified by M. Ratto
% New output argument aK: 1-step to nk-stpe ahed predictions)
% New input argument nk: max order of predictions in aK
% New global variable id_ where the DKF stops (common with
% diffuselikelihood3)
% New icc variable to count number of iterations for Finf steps
% Pstar % Pinf simmetric
% New termination of DKF iterations based on id_
% Li also stored during DKF iterations !!
% some bugs corrected in the DKF part of the smoother (Z matrix and
% alphahat)
%
% stephane.adjemian@cepremap.cnrs.fr [09-16-2004]
% 
%   See "Fast Filtering and Smoothing for Multivariate State Space
%   Models", S.J. Koopman and J. Durbin (2000, in Journal of Time Series 
%   Analysis, vol. 21(3), pp. 281-296).  

global options_

nk = options_.nk;
spinf   	= size(Pinf1);
spstar  	= size(Pstar1);
v       	= zeros(pp,smpl);
a       	= zeros(mm,smpl+1);
a1			= a;
aK          = zeros(nk,mm,smpl+nk);
Fstar   	= zeros(pp,smpl);
Finf		= zeros(pp,smpl);
Ki       	= zeros(mm,pp,smpl);
Li      	= zeros(mm,mm,pp,smpl);
Linf    	= zeros(mm,mm,pp,smpl);
L0      	= zeros(mm,mm,pp,smpl);
Kstar   	= zeros(mm,pp,smpl);
P       	= zeros(mm,mm,smpl+1);
P1			= P;
Pstar   	= zeros(spstar(1),spstar(2),smpl+1); Pstar(:,:,1) = Pstar1;
Pinf    	= zeros(spinf(1),spinf(2),smpl+1); Pinf(:,:,1) = Pinf1;
Pstar1 		= Pstar;
Pinf1  		= Pinf;
crit   	 	= options_.kalman_tol;
crit1       = 1.e-6;
steady  	= smpl;
rr      	= size(Q,1);
QQ      	= R*Q*transpose(R);
QRt			= Q*transpose(R);
alphahat   	= zeros(mm,smpl);
etahat	   	= zeros(rr,smpl);
epsilonhat      = zeros(size(Y));
r 		   	= zeros(mm,smpl);

Z = zeros(pp,mm);
for i=1:pp;
	Z(i,mf(i)) = 1;
end

t = 0;
icc=0;
newRank	  = rank(Pinf(:,:,1),crit1);
while newRank & t < smpl
  t = t+1;
  a1(:,t) = a(:,t);
  Pstar(:,:,t)=tril(Pstar(:,:,t))+transpose(tril(Pstar(:,:,t),-1));
  Pinf(:,:,t)=tril(Pinf(:,:,t))+transpose(tril(Pinf(:,:,t),-1));
  Pstar1(:,:,t) = Pstar(:,:,t);
  Pinf1(:,:,t) = Pinf(:,:,t);
  for i=1:pp
    v(i,t) 	= Y(i,t)-a(mf(i),t)-trend(i,t);
    Fstar(i,t) 	= Pstar(mf(i),mf(i),t) + H(i,i);
    Finf(i,t)	= Pinf(mf(i),mf(i),t);
    Kstar(:,i,t) 	= Pstar(:,mf(i),t);
    if Finf(i,t) > crit & newRank,  % original MJ: if Finf(i,t) > crit
      icc=icc+1;
      Kinf(:,i,t)	= Pinf(:,mf(i),t);
      Linf(:,:,i,t)  	= eye(mm) - Kinf(:,i,t)*Z(i,:)/Finf(i,t);
      L0(:,:,i,t)  	= (Kinf(:,i,t)*Fstar(i,t)/Finf(i,t) - Kstar(:,i,t))*Z(i,:)/Finf(i,t);
      a(:,t)		= a(:,t) + Kinf(:,i,t)*v(i,t)/Finf(i,t);
      Pstar(:,:,t)	= Pstar(:,:,t) + ...
	  Kinf(:,i,t)*transpose(Kinf(:,i,t))*Fstar(i,t)/(Finf(i,t)*Finf(i,t)) - ...
	  (Kstar(:,i,t)*transpose(Kinf(:,i,t)) +...
	   Kinf(:,i,t)*transpose(Kstar(:,i,t)))/Finf(i,t);
      Pinf(:,:,t)	= Pinf(:,:,t) - Kinf(:,i,t)*transpose(Kinf(:,i,t))/Finf(i,t);
      Pstar(:,:,t)=tril(Pstar(:,:,t))+transpose(tril(Pstar(:,:,t),-1));
      Pinf(:,:,t)=tril(Pinf(:,:,t))+transpose(tril(Pinf(:,:,t),-1));
      % new terminiation criteria by M. Ratto
      P0=Pinf(:,:,t);
      %             newRank = any(diag(P0(mf,mf))>crit);
      %             if newRank==0, options_.diffuse_d = i; end,
      if ~isempty(options_.diffuse_d),  
	newRank = (icc<options_.diffuse_d);  
	%if newRank & any(diag(P0(mf,mf))>crit)==0; 
	if newRank & (any(diag(P0(mf,mf))>crit)==0 & rank(P0,crit1)==0); 
	  disp('WARNING!! Change in OPTIONS_.DIFFUSE_D in univariate DKF')
	  options_.diffuse_d = icc;
	  newRank=0;
	end
      else
	%newRank = any(diag(P0(mf,mf))>crit);                 
	newRank = (any(diag(P0(mf,mf))>crit) | rank(P0,crit1));                 
	if newRank==0, 
	  options_.diffuse_d = icc;
	end                    
      end,
      if newRank==0, 
	options_.diffuse_d=i;
      end                    
      % end new terminiation criteria by M. Ratto
    else 
      %% Note that : (1) rank(Pinf)=0 implies that Finf = 0, (2) outside this loop (when for some i and t the condition
      %% rank(Pinf)=0 is satisfied we have P = Pstar and F = Fstar and (3) Finf = 0 does not imply that
      %% rank(Pinf)=0. [stéphane,11-03-2004].	  
      Li(:,:,i,t)    = eye(mm)-Kstar(:,i,t)*Z(i,:)/Fstar(i,t);  % we need to store Li for DKF smoother
      a(:,t) 		= a(:,t) + Kstar(:,i,t)*v(i,t)/Fstar(i,t);
      Pstar(:,:,t)	= Pstar(:,:,t) - Kstar(:,i,t)*transpose(Kstar(:,i,t))/Fstar(i,t);
      Pstar(:,:,t)=tril(Pstar(:,:,t))+transpose(tril(Pstar(:,:,t),-1));
    end
  end
  a(:,t+1) 	 	= T*a(:,t);
  for jnk=1:nk,
    aK(jnk,:,t+jnk) 	 	= T^jnk*a(:,t);
  end
  Pstar(:,:,t+1)	= T*Pstar(:,:,t)*transpose(T)+ QQ;
  Pinf(:,:,t+1)	= T*Pinf(:,:,t)*transpose(T);
  P0=Pinf(:,:,t+1);
  if newRank,
    %newRank = any(diag(P0(mf,mf))>crit);
    newRank	  = rank(P0,crit1);
  end
end


d = t;
P(:,:,d+1) = Pstar(:,:,d+1);
Linf  = Linf(:,:,:,1:d);
L0  = L0(:,:,:,1:d);
Fstar = Fstar(:,1:d);
Finf = Finf(:,1:d);
Kstar = Kstar(:,:,1:d);
Pstar = Pstar(:,:,1:d);
Pinf  = Pinf(:,:,1:d);
Pstar1 = Pstar1(:,:,1:d);
Pinf1  = Pinf1(:,:,1:d);
notsteady = 1;
while notsteady & t<smpl
  t = t+1;
  a1(:,t) = a(:,t);
  P(:,:,t)=tril(P(:,:,t))+transpose(tril(P(:,:,t),-1));
  P1(:,:,t) = P(:,:,t);
  for i=1:pp
    v(i,t)  = Y(i,t) - a(mf(i),t) - trend(i,t);
    Fi(i,t) = P(mf(i),mf(i),t);
    Ki(:,i,t) = P(:,mf(i),t) + H(i,i);
    if Fi(i,t) > crit
      Li(:,:,i,t)    = eye(mm)-Ki(:,i,t)*Z(i,:)/Fi(i,t);
      a(:,t) = a(:,t) + Ki(:,i,t)*v(i,t)/Fi(i,t);
      P(:,:,t) = P(:,:,t) - Ki(:,i,t)*transpose(Ki(:,i,t))/Fi(i,t);
      P(:,:,t)=tril(P(:,:,t))+transpose(tril(P(:,:,t),-1));
    end
  end
  a(:,t+1) = T*a(:,t);
  for jnk=1:nk,
    aK(jnk,:,t+jnk) 	 	= T^jnk*a(:,t);
  end
  P(:,:,t+1) = T*P(:,:,t)*transpose(T) + QQ;
  notsteady   = ~(max(max(abs(P(:,:,t+1)-P(:,:,t))))<crit);
end
P_s=tril(P(:,:,t))+transpose(tril(P(:,:,t),-1));
Fi_s = Fi(:,t);
Ki_s = Ki(:,:,t);
L_s  =Li(:,:,:,t);
if t<smpl
  t_steady = t+1;
  P  = cat(3,P(:,:,1:t),repmat(P(:,:,t),[1 1 smpl-t_steady+1]));
  Fi = cat(2,Fi(:,1:t),repmat(Fi_s,[1 1 smpl-t_steady+1]));
  Li  = cat(4,Li(:,:,:,1:t),repmat(L_s,[1 1 smpl-t_steady+1]));
  Ki  = cat(3,Ki(:,:,1:t),repmat(Ki_s,[1 1 smpl-t_steady+1]));
end
while t<smpl
  t=t+1;
  a1(:,t) = a(:,t);
  for i=1:pp
    v(i,t)      = Y(i,t) - a(mf(i),t) - trend(i,t);
    if Fi_s(i) > crit
      a(:,t) = a(:,t) + Ki_s(:,i)*v(i,t)/Fi_s(i);
    end
  end
  a(:,t+1) = T*a(:,t);
  for jnk=1:nk,
    aK(jnk,:,t+jnk)	= T^jnk*a(:,t);
  end
end
a1(:,t+1) = a(:,t+1);
ri=r;
t = smpl+1;
while t>d+1 & t>2,
  t = t-1;
  for i=pp:-1:1
    if Fi(i,t) > crit
      ri(:,t)=transpose(Z(i,:))/Fi(i,t)*v(i,t)+transpose(Li(:,:,i,t))*ri(:,t);
    end
  end
  r(:,t-1) = ri(:,t);
  alphahat(:,t)	= a1(:,t) + P1(:,:,t)*r(:,t-1);
  etahat(:,t)		= QRt*r(:,t);
  ri(:,t-1) = transpose(T)*ri(:,t);
end
if d
  r0 = zeros(mm,d); r0(:,d) = ri(:,d);
  r1 = zeros(mm,d);
  for t = d:-1:2
    for i=pp:-1:1
      if Finf(i,t) > crit & ~(t==d & i>options_.diffuse_d),  % use of options_.diffuse_d to be sure of DKF termination
					     %r1(:,t) = transpose(Z)*v(:,t)/Finf(i,t) + ... BUG HERE in transpose(Z)
					     r1(:,t) = transpose(Z(i,:))*v(i,t)/Finf(i,t) + ...
						       transpose(L0(:,:,i,t))*r0(:,t) + transpose(Linf(:,:,i,t))*r1(:,t);
					     r0(:,t) = transpose(Linf(:,:,i,t))*r0(:,t);
      elseif Fstar(i,t) > crit % step needed whe Finf == 0
	r0(:,t)=transpose(Z(i,:))/Fstar(i,t)*v(i,t)+Li(:,:,i,t)'*r0(:,t);
      end
    end
    alphahat(:,t)	= a1(:,t) + Pstar1(:,:,t)*r0(:,t) + Pinf1(:,:,t)*r1(:,t);
    r(:,t-1)		= r0(:,t);
    etahat(:,t)		= QRt*r(:,t);
    r0(:,t-1) = transpose(T)*r0(:,t);
    r1(:,t-1) = transpose(T)*r1(:,t);
  end
  r0_0 = r0(:,1);
  r1_0 = r1(:,1);
  for i=pp:-1:1
    if Finf(i,1) > crit,
      %r1_0 = transpose(Z)*v(:,1)/Finf(i,1) + ... %bug with Z here
      r1_0 = transpose(Z(i,:))*v(i,1)/Finf(i,1) + ...
	     transpose(L0(:,:,i,1))*r0_0 + transpose(Linf(:,:,i,1))*r1_0;
      r0_0 = transpose(Linf(:,:,i,1))*r0_0;
    elseif Fstar(i,1) > crit, % step needed when Finf=0
      r0_0=transpose(Z(i,:))/Fstar(i,1)*v(i,1)+Li(:,:,i,1)'*r0_0;
    end
  end
  %alphahat(:,1)  	= a(:,1) + Pstar(:,:,1)*r0_0 + Pinf(:,:,1)*r1_0; %this line is buggy
  alphahat(:,1)  	= a1(:,1) + Pstar1(:,:,1)*r0_0 + Pinf1(:,:,1)*r1_0;
  etahat(:,1)		= QRt*r(:,1);
else
  r0 = ri(:,1);
  for i=pp:-1:1
    if Fi(i,1) > crit
      r0=transpose(Z(i,:))/Fi(i,1)*v(i,1)+transpose(Li(:,:,i,1))*r0;
    end
  end 
  %alphahat(:,1)	= a(:,1) + P(:,:,1)*r0;  % this line is buggy
  alphahat(:,1)	= a1(:,1) + P1(:,:,1)*r0;
  etahat(:,1)	= QRt*r(:,1);
end
epsilonhat = Y-alphahat(mf,:)-trend;

