function w=matrictint(S,XXi,T)
%function w=matrictint(S,XXi,T)
%  S:  usually sample cross product matrix of LS residuals
% XXi:  inv(X'X) matrix for rhs variables
%  T:  number of observations
%  w:  log of integrated posterior for SUR or RF VAR with det(Sigma)^(-(m+1)/2) Jeffreys-like prior
%  To get the log of the integral of the likelihood for a VAR with T observations, 
%   k rhs variables in each equation, and m equations, set T=T-m-1 and subtract .5*m*(m+1)*log(2*pi).
% We are integrating the exponential of -.5*T*m*log(2*pi)-.5*(T+m+1)*log(det(Sigma))-.5*trace(Sigma\S(beta)).
k=size(XXi,1);
m=size(S,1);
[cx,p]=chol(XXi);
[cs,q]=chol(S);
%cx=cschol(XXi);
%cs=cschol(S);
if any(diag(cx)<100*eps)
    error('singular XXi')
end
if any(diag(cs<100*eps))
    error('singular S')
end
w=(-T+k+(m-1)/2)*m*.5*log(pi)-(T-k)*sum(log(diag(cs)))+m*sum(log(diag(cx)))+ggammaln(m,(T-k)/2);

function lgg=ggammaln(m,ndf)
%function gg=ggamma(m,ndf)
% From 8.2.22 on p.427 of Box and Tiao, this is the log of generalized
% gamma divided by gamma(.5)^(.5*m*(m-1))
if ndf<=(m-1)/2
    error('too few df in ggammaln')
else
    %lgg=.5*m*(m-1)*gammaln(.5); % normalizing factor not used in Wishart integral
    garg=ndf+.5*(0:-1:1-m);
    lgg=sum(gammaln(garg));
end
    
