function [alphahat,epsilonhat,etahat,a, aK] = DiffuseKalmanSmootherH1(T,R,Q,H,Pinf1,Pstar1,Y,trend,pp,mm,smpl,mf)
% modified by M. Ratto:
% new output argument aK (1-step to k-step predictions)
% new options_.nk: the max step ahed prediction in aK (default is 4)
% new crit1 value for rank of Pinf
% it is assured that P is symmetric 
%
% stephane.adjemian@cepremap.cnrs.fr [09-16-2004]
% 
%   See "Filtering and Smoothing of State Vector for Diffuse State Space
%   Models", S.J. Koopman and J. Durbin (2003, in Journal of Time Series 
%   Analysis, vol. 24(1), pp. 85-98).  

global options_

nk = options_.nk;
spinf   	= size(Pinf1);
spstar  	= size(Pstar1);
v       	= zeros(pp,smpl);
a       	= zeros(mm,smpl+1);
iF      	= zeros(pp,pp,smpl);
Fstar   	= zeros(pp,pp,smpl);
iFinf   	= zeros(pp,pp,smpl);
K       	= zeros(mm,pp,smpl);
L       	= zeros(mm,mm,smpl);
Linf    	= zeros(mm,mm,smpl);
Kstar   	= zeros(mm,pp,smpl);
P       	= zeros(mm,mm,smpl+1);
Pstar   	= zeros(spstar(1),spstar(2),smpl+1); Pstar(:,:,1) = Pstar1;
Pinf    	= zeros(spinf(1),spinf(2),smpl+1); Pinf(:,:,1) = Pinf1;
crit    	= options_.kalman_tol;
crit1       = 1.e-8;
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
while rank(Pinf(:,:,t+1),crit1) & t<smpl
    t = t+1;
    v(:,t) 		 	= Y(:,t) - a(mf,t) - trend(:,t);
    if rcond(Pinf(mf,mf,t)) < crit
    	return		
    end
    iFinf(:,:,t) 	= inv(Pinf(mf,mf,t));
    Kinf(:,:,t)	 	= T*Pinf(:,mf,t)*iFinf(:,:,t);
    a(:,t+1) 	 	= T*a(:,t) + Kinf(:,:,t)*v(:,t);
    aK(1,:,t+1) 	 	= a(:,t+1);
    for jnk=2:nk,
        aK(jnk,:,t+jnk) 	 	= T^(jnk-1)*a(:,t+1);
    end
    Linf(:,:,t)  	= T - Kinf(:,:,t)*Z;
    Fstar(:,:,t) 	= Pstar(mf,mf,t) + H;
    Kstar(:,:,t) 	= (T*Pstar(:,mf,t)-Kinf(:,:,t)*Fstar(:,:,t))*iFinf(:,:,t);
    Pstar(:,:,t+1)	= T*Pstar(:,:,t)*transpose(T)-T*Pstar(:,mf,t)*transpose(Kinf(:,:,t))-Kinf(:,:,t)*Pinf(mf,mf,t)*transpose(Kstar(:,:,t)) + QQ;
    Pinf(:,:,t+1)	= T*Pinf(:,:,t)*transpose(T)-T*Pinf(:,mf,t)*transpose(Kinf(:,:,t));
end
d = t;
P(:,:,d+1) = Pstar(:,:,d+1);
iFinf = iFinf(:,:,1:d);
Linf  = Linf(:,:,1:d);
Fstar = Fstar(:,:,1:d);
Kstar = Kstar(:,:,1:d);
Pstar = Pstar(:,:,1:d);
Pinf  = Pinf(:,:,1:d);
notsteady = 1;
while notsteady & t<smpl
    t = t+1;
    v(:,t)      = Y(:,t) - a(mf,t) - trend(:,t);
    P(:,:,t)=tril(P(:,:,t))+transpose(tril(P(:,:,t),-1));
    if rcond(P(mf,mf,t) + H) < crit
    	return		
    end    
    iF(:,:,t)   = inv(P(mf,mf,t) + H);
    K(:,:,t)    = T*P(:,mf,t)*iF(:,:,t);
    L(:,:,t)    = T-K(:,:,t)*Z;
    a(:,t+1)    = T*a(:,t) + K(:,:,t)*v(:,t);    
    aK(1,:,t+1) 	 	= a(:,t+1);
    for jnk=2:nk,
        aK(jnk,:,t+jnk) 	 	= T^(jnk-1)*a(:,t+1);
    end
    P(:,:,t+1)  = T*P(:,:,t)*transpose(T)-T*P(:,mf,t)*transpose(K(:,:,t)) + QQ;
    notsteady   = ~(max(max(abs(P(:,:,t+1)-P(:,:,t))))<crit);
end
K_s = K(:,:,t);
iF_s = iF(:,:,t);
P_s = P(:,:,t+1);
if t<smpl
	t_steady = t+1;
	P  = cat(3,P(:,:,1:t),repmat(P(:,:,t),[1 1 smpl-t_steady+1]));
	iF = cat(3,iF(:,:,1:t),repmat(inv(P_s(mf,mf)+H),[1 1 smpl-t_steady+1]));
	L  = cat(3,L(:,:,1:t),repmat(T-K_s*Z,[1 1 smpl-t_steady+1]));
	K  = cat(3,K(:,:,1:t),repmat(T*P_s(:,mf)*iF_s,[1 1 smpl-t_steady+1]));
end
while t<smpl
    t=t+1;
    v(:,t) = Y(:,t) - a(mf,t) - trend(:,t);
    a(:,t+1) = T*a(:,t) + K_s*v(:,t);
    aK(1,:,t+1) 	 	= a(:,t+1);
    for jnk=2:nk,
        aK(jnk,:,t+jnk) 	 	= T^(jnk-1)*a(:,t+1);
    end
end
t = smpl+1;
while t>d+1 & t>2
	t = t-1;
    r(:,t-1) = transpose(Z)*iF(:,:,t)*v(:,t) + transpose(L(:,:,t))*r(:,t);
    alphahat(:,t)	= a(:,t) + P(:,:,t)*r(:,t-1);
	etahat(:,t)		= QRt*r(:,t);
end
if d
	r0 = zeros(mm,d); r0(:,d) = r(:,d);
	r1 = zeros(mm,d);
	for t = d:-1:2
    	r0(:,t-1) = transpose(Linf(:,:,t))*r0(:,t);
		r1(:,t-1) = transpose(Z)*(iFinf(:,:,t)*v(:,t)-transpose(Kstar(:,:,t))*r0(:,t)) + transpose(Linf(:,:,t))*r1(:,t);
		alphahat(:,t)	= a(:,t) + Pstar(:,:,t)*r0(:,t-1) + Pinf(:,:,t)*r1(:,t-1);
		etahat(:,t)		= QRt*r0(:,t);
	end
	r0_0 = transpose(Linf(:,:,1))*r0(:,1);
	r1_0 = transpose(Z)*(iFinf(:,:,1)*v(:,1)-transpose(Kstar(:,:,1))*r0(:,1)) + transpose(Linf(:,:,1))*r1(:,1);
	alphahat(:,1)  	= a(:,1) + Pstar(:,:,1)*r0_0 + Pinf(:,:,1)*r1_0;
	etahat(:,1)		= QRt*r0(:,1);
else
    r0 = transpose(Z)*iF(:,:,1)*v(:,1) + transpose(L(:,:,1))*r(:,1);
    alphahat(:,1)	= a(:,1) + P(:,:,1)*r0;
    etahat(:,1)	= QRt*r(:,1);
end
epsilonhat = Y-alphahat(mf,:)-trend;
