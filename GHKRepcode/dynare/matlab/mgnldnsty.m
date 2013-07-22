function w=mgnldnsty(ydata,lags,xdata,breaks,lambda,mu,mnprior,vprior,train,flat)
%function w=mgnldnsty(ydata,lags,xdata,breaks,lambda,mu,mnprior,vprior,train,flat)
% ydata:        endogenous variable data matrix, including initial condition dates.
% xdata:        exogenous variable data matrix, including initial condition dates.
% breaks:       breaks in the data.  The first lags data points after a break are used
%               as new initial conditions, not data points for the fit.
% lambda:       weight on the co-persistence prior dummy observation.  (5 is reasonable)
%               lambda>0 => x variables included; lambda<0 => x variables excluded;
% mu:           weight on variable-by-variable sum of coeffs dummy obs. (1 is reasonable)
% mnprior.tight:weight on the Minnesota prior dummies.  Prior std dev on first lag is
%               1/mnprior.tight
% mnprior.decay:prior std dev on lag j is 1/j^decay
% vprior.sig:   vector of nv prior std dev's of equation shocks
% vprior.w:     weight on vcv dummies.  (1 is reasonable; higher values tighten up.)
% train:        If present and non-zero, this is the point in the sample at which the
%               "training sample" ends.  Prior x likelihood to this point is weighted to
%               integrate to 1, and therefore is treated as if it were itself the prior.
%               To do a pure training sample prior, set lambda=mu=0, mnprior=vprior=[],
%               train>lags.  
%
%flat:          Even with lambda=mu=0, vprior=mnprior=[], det(Sigma)^(-(nv+1)/2) is used
%               as a "prior", unless flat=1.  flat, if present, must be 1 or 0.
%               flat=1 is likely not to work unless train is reasonably large.
if nargin<10,flat=0;end
[T,nv]=size(ydata);
[Tx,nx]=size(xdata);
if Tx ~= T, error('ydata and xdata length mismatch'),end
[ydum,xdum,pbreaks]=varprior(nv,nx,lags,mnprior,vprior);
var=rfvar3([ydata;ydum],lags,[xdata;xdum],[breaks;T;T+pbreaks],lambda,mu);
Tu=size(var.u,1);
w=matrictint(var.u'*var.u,var.xxi,Tu-flat*(nv+1))-flat*.5*nv*(nv+1)*log(2*pi);
if nargin>8
    if ~isempty(train) & train>0
        if train <= lags
            error('end of training sample <= # of lags')
        end
        Tp=train;
        tbreaks=breaks(find(breaks<train));
    else
        Tp=lags;
        tbreaks=[];
    end
else
    Tp=lags;
    tbreaks=[];
end
varp=rfvar3([ydata(1:Tp,:);ydum],lags,[xdata(1:Tp);xdum],[tbreaks;Tp;Tp+pbreaks],lambda,mu);
Tup=size(varp.u,1);
wp=matrictint(varp.u'*varp.u,varp.xxi,Tup-flat*(nv+1)/2)-flat*.5*nv*(nv+1)*log(2*pi);
w=w-wp;