function compdist(xparam1, x2, pltopt, figurename)
global bayestopt_  estim_params_ lgx_ lgy_ options_

% NOTE: If pltopt ~= 'All' compdist.m just draws prior densities.  

%% Set density estimation parameters:
number_of_grid_points = 2^9;      % 2^9 = 512 !... Must be a power of two.
bandwidth = 0;                    % Rule of thumb optimal bandwidth parameter.
kernel_function = 'gaussian';     % You can switch to: 'epanechnikov', 'quartic', 'triangle', 
% 'triweight', 'uniform' or 'cosinus' kernels (iff bandwidth=0, see posterior_density_estimate.m). 
truncprior = 10^(-3);             
  
  
npar=length(xparam1);
nruns=length(x2);
icol=ceil(sqrt(npar));
iraw=icol;
if (icol-1)*(icol-2)>=npar
    iraw = icol-2;
    icol=icol-1;
elseif (icol)*(icol-2)>=npar
    iraw = icol-2;
elseif icol*(icol-1)>=npar
    iraw=icol-1;
end

pmean=bayestopt_.pmean;
pshape=bayestopt_.pshape; 
p1 = bayestopt_.p1;
p2 = bayestopt_.p2;
p3 = bayestopt_.p3;
p4 = bayestopt_.p4;
  
figure('Name',figurename)
for i=1:npar;
  if i<=estim_params_.nvx
    vname = deblank(lgx_(estim_params_.var_exo(i,1),:));
    nam=['SE_{',vname,'}'];
  elseif  i<=(estim_params_.nvx+estim_params_.nvn)
    deblank(options_.varobs(estim_params_.var_endo(i-estim_params_.nvx,1),:));
    nam=['SE_{EOBS_',vname,'}'];
  elseif  i<=(estim_params_.nvx+estim_params_.nvn+estim_params_.ncx)
    j = i - (estim_params_.nvx+estim_params_.nvn);
    k1 = estim_params_.corrx(j,1);
    k2 = estim_params_.corrx(j,2);
    vname = [deblank(lgx_(k1,:)) ',' deblank(lgx_(k2,:))];
    nam=['CC_{',vname,'}'];
  elseif  i<=(estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+ ...
              estim_params_.ncn)
    j = i - (estim_params_.nvx+estim_params_.nvn+estim_params_.ncx);
    k1 = estim_params_.corrn(j,1);
    k2 = estim_params_.corrn(j,2);
    vname = [deblank(lgx_(k1,:)) ',' deblank(lgx_(k2,:))];
    nam=['CC_{EOBS_',vname,'}'];
  else
    j = i - (estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+ ...
      estim_params_.ncn);
    nam = deblank(estim_params_.param_names(j,:));
  end
  subplot(iraw, icol, i);   
  if strcmpi(pltopt,'all'); % Estimation of the density...
      [abscissa,ff,h] = posterior_density_estimate(x2(round(options_.mh_drop*nruns):end,i),...
                      number_of_grid_points,bandwidth,kernel_function);
      plot(abscissa,ff,'-k','linewidth',2); 
  end;
  a = 0;
  b = 0;
  if pshape(i) == 1;     %/* BETA Prior */
      density = inline('((bb-x).^(b-1)).*(x-aa).^(a-1)./(beta(a,b)*(bb-aa)^(a+b-1))','x','a','b','aa','bb');
      mu = (p1(i)-p3(i))/(p4(i)-p3(i));
      stdd = p2(i)/(p4(i)-p3(i));
      a = (1-mu)*mu^2/stdd^2 - mu;
      b = a*(1/mu - 1);
      aa = p3(i);
      bb = p4(i);
      infbound = qbeta(truncprior,a,b)*(bb-aa)+aa;
      supbound = qbeta(1-truncprior,a,b)*(bb-aa)+aa;
      stepsize = (supbound-infbound)/200;
      abscissa = infbound:stepsize:supbound;
      f = density(abscissa,a,b,aa,bb);
      if strcmpi(pltopt,'all');
          top = max([max(ff);max(f)]);
      end;
  elseif pshape(i) == 2; %/* GAMMA PRIOR */
%      density = inline('((x/b).^(a-1)).*exp(-x/b)*inv(b*gamma(a))','x','a','b');
      mu = p1(i)-p3(i);
      b  = p2(i)^2/mu;
      a  = mu/b;
      infbound = mj_qgamma(truncprior,a)*b; 
      supbound = mj_qgamma(1-truncprior,a)*b;
      stepsize = (supbound-infbound)/200;
      abscissa = infbound:stepsize:supbound;
%      f = density(abscissa,a,b);
      f = exp(lpdfgam(abscissa,a,b));
      abscissa = abscissa + p3(i);
      if strcmpi(pltopt,'all');
          top = max([max(ff);max(f)]);
      end;
  elseif pshape(i) == 3; %/* GAUSSIAN PRIOR */
      density = inline('inv(sqrt(2*pi)*b)*exp(-0.5*((x-a)/b).^2)','x','a','b');
      a = p1(i);
      b = p2(i);
      infbound = qnorm(truncprior,a,b); 
      supbound = qnorm(1-truncprior,a,b);
      stepsize = (supbound-infbound)/200;
      abscissa = infbound:stepsize:supbound;
      f = density(abscissa,a,b);  
      if strcmpi(pltopt,'all');
          top = max([max(ff);max(f)]);
      end;
  elseif pshape(i) == 4; %/* INVGAMMA PRIOR type 1 */
      density = inline('2*inv(gamma(nu/2))*(x.^(-nu-1))*((s/2)^(nu/2)).*exp(-s./(2*x.^2))','x','s','nu');
      nu = p2(i);
      s  = p1(i);
      a  = nu/2;
      b  = 2/s;
      infbound = 1/sqrt(mj_qgamma(1-10*truncprior,a)*b); 
      supbound = 1/sqrt(mj_qgamma(10*truncprior,a)*b);
      stepsize = (supbound-infbound)/200;
      abscissa = infbound:stepsize:supbound;
      f = density(abscissa,s,nu);  
      if strcmpi(pltopt,'all');
          top = max([max(ff);max(f)]);
      end;
  elseif pshape(i) == 5; %/* UNIFORM PRIOR */
      density = inline('(x.^0)/(b-a)','x','a','b');
      a  = p1(i);
      b  = p2(i);
      infbound = a; 
      supbound = b;
      stepsize = (supbound-infbound)/200;
      abscissa = infbound:stepsize:supbound;
      f = density(abscissa,a,b);  
      if strcmpi(pltopt,'all');
          top = max([max(ff);max(f)]);
      end;
  elseif pshape(i) == 6; %/*  INVGAMMA PRIOR type 2 */        
      density = inline('inv(gamma(nu/2))*(x.^(-.5(nu+2)))*((s/2)^(nu/2)).*exp(-s./(2*x))','x','s','nu');
      nu = p2(i);
      s  = p1(i);
      a  = nu/2;
      b  = 2/s;
      infbound = 1/(qgamma(1-truncprior,a)*b); 
      supbound = 1/(qgamma(truncprior,a)*b);
      stepsize = (supbound-infbound)/200;
      abscissa = infbound:stepsize:supbound;
      f = density(abscissa,s,nu);  
      if strcmpi(pltopt,'all');
          top = max([max(ff);max(f)]);
      end;
  end;
  hold on;
  k = [1:length(f)];
  if pshape(i) ~= 5 
    [junk,k1] = max(f);
    if k1 == 1 | k1 == length(f)
      k = find(f < 10);
    end
  end
  hh = plot(abscissa(k),f(k),'-k','linewidth',2);
  set(hh,'color',[0.7 0.7 0.7]);
  if strcmpi(pltopt,'all');
      plot( [xparam1(i) xparam1(i)], [0,top], '--g', 'linewidth', 2);
  end;
  title(nam,'Interpreter','none');
  hold off;
end;
drawnow

% 12/01/03 MJ adapted from M. Ratto's version
% 02/16/04 SA correction to the generalized beta distibution over interval [aa,bb].
