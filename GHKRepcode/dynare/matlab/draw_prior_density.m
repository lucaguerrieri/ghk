
function [x,f,abscissa,dens,binf,bsup] = draw_prior_density(indx);
% stephane.adjemian@cepremap.cnrs.fr [07-15-2004]

global bayestopt_

pmean   = bayestopt_.pmean;
pshape  = bayestopt_.pshape; 
p1      = bayestopt_.p1;
p2      = bayestopt_.p2;
p3      = bayestopt_.p3;
p4      = bayestopt_.p4;

truncprior = 10^(-3);

if pshape(indx) == 1     %/* BETA Prior */
    density = inline('((bb-x).^(b-1)).*(x-aa).^(a-1)./(beta(a,b)*(bb-aa)^(a+b-1))','x','a','b','aa','bb');
    mu = (p1(indx)-p3(indx))/(p4(indx)-p3(indx));
    stdd = p2(indx)/(p4(indx)-p3(indx));
    a = (1-mu)*mu^2/stdd^2 - mu;
    b = a*(1/mu-1);
    aa = p3(indx);
    bb = p4(indx);
    infbound = qbeta(truncprior,a,b)*(bb-aa)+aa;
    supbound = qbeta(1-truncprior,a,b)*(bb-aa)+aa;
    stepsize = (supbound-infbound)/200;
    abscissa = infbound:stepsize:supbound;
    dens = density(abscissa,a,b,aa,bb);
elseif pshape(indx) == 2  %/* GAMMA PRIOR */
    mu = p1(indx)-p3(indx);
    b  = p2(indx)^2/mu;
    a  = mu/b;
    infbound = mj_qgamma(truncprior,a)*b; 
    supbound = mj_qgamma(1-truncprior,a)*b;
    stepsize = (supbound-infbound)/200;
    abscissa = infbound:stepsize:supbound;
    dens = exp(lpdfgam(abscissa,a,b));
    abscissa = abscissa + p3(indx);
elseif pshape(indx) == 3  %/* GAUSSIAN PRIOR */
    density = inline('inv(sqrt(2*pi)*b)*exp(-0.5*((x-a)/b).^2)','x','a','b');
    a = p1(indx);
    b = p2(indx);
    infbound = qnorm(truncprior,a,b); 
    supbound = qnorm(1-truncprior,a,b);
    stepsize = (supbound-infbound)/200;
    abscissa = infbound:stepsize:supbound;
    dens = density(abscissa,a,b);  
elseif pshape(indx) == 4  %/* INVGAMMA PRIOR type 1 */
    density = inline('2*inv(gamma(nu/2))*(x.^(-nu-1))*((s/2)^(nu/2)).*exp(-s./(2*x.^2))','x','s','nu');
    nu = p2(indx);
    s  = p1(indx);
    a  = nu/2;
    b  = 2/s;
    infbound = 1/sqrt(mj_qgamma(1-10*truncprior,a)*b); 
    supbound = 1/sqrt(mj_qgamma(10*truncprior,a)*b);
    stepsize = (supbound-infbound)/200;
    abscissa = infbound:stepsize:supbound;
    dens = density(abscissa,s,nu);  
elseif pshape(indx) == 5  %/* UNIFORM PRIOR */
    density = inline('(x.^0)/(b-a)','x','a','b');
    a  = p1(indx);
    b  = p2(indx);
    infbound = a; 
    supbound = b;
    stepsize = (supbound-infbound)/200;
    abscissa = infbound:stepsize:supbound;
    dens = density(abscissa,a,b);  
elseif pshape(indx) == 6  %/*  INVGAMMA PRIOR type 2 */        
    density = inline('inv(gamma(nu/2))*(x.^(-.5*(nu+2)))*((s/2)^(nu/2)).*exp(-s./(2*x))','x','s','nu');
    nu = p2(indx);
    s  = p1(indx);
    a  = nu/2;
    b  = 2/s;
    infbound = 1/(qgamma(1-truncprior,a)*b); 
    supbound = 1/(qgamma(truncprior,a)*b);
    stepsize = (supbound-infbound)/200;
    abscissa = infbound:stepsize:supbound;
    dens = density(abscissa,s,nu);  
end 

k = [1:length(dens)];
if pshape(indx) ~= 5 
    [junk,k1] = max(dens);
    if k1 == 1 | k1 == length(dens)
        k = find(dens < 10);
    end
end
binf = abscissa(k(1));
bsup = abscissa(k(length(k)));
x = abscissa(k);
f = dens(k);