function EvalWelfare(vName,pltOpt,NumberOfDraws)
% stephane.adjemian@ens.fr
%
% vName is a string designing the name of the welfare variable in the mod file. 
%   
global options_ oo_ lgy_ ys_ fname_ dr_ info

options_.order = 2;
options_ = set_default_option(options_,'dr_algo',0);
options_ = set_default_option(options_,'simul_algo',0);
subindx = subset();

if isempty(NumberOfDraws)
  B = 500;
else
  B = NumberOfDraws; 
end

% [1] How many mh files and simulations ?
if options_.mh_nblck == 1
  files = eval(['dir(''' fname_ '_mh*.mat'');']);
else% More than one block
  files = eval(['dir(''' fname_ '_mh*_blck1.mat'');']);
end
nn = size(files,1);
METRO = zeros(nn,3);
for i = 1:nn
  file = files(i);
  eval(['load ' file.name ';'])
  METRO(i,1) = i-1;
  METRO(i,2) = size(x2,1);
  if i>1
    METRO(i,3) = METRO(i,2)+METRO(i-1,3);
  else
    METRO(i,3) = METRO(i,2);
  end
end
TotalNumberOfDraws = METRO(end,3);
% [2] Take some draws from the metropolis and compute the welfare.
% [2.1] Set some parameters:
B = min(B,TotalNumberOfDraws);
WelfDistribution = zeros(B,1);
kk = dr_.order_var;
windx = strmatch(vName,lgy_(kk,:),'exact');
% [2.2] Get the posterior mean:
deep = get_posterior_parameters('posterior_mean');
% [2.3] Compute the posterior distribution of Welfare:
hfid = waitbar(0,'Posterior welfare distribution...');
compt = 0;
ys = ys_;
for i=1:B
  linea = 1+floor(rand*TotalNumberOfDraws);
  tmp = find(METRO(:,3)<linea);
  if isempty(tmp)
    FileNumber = 1;
    linee = linea;
  else
    FileNumber = tmp(end)+1;
    linee = linea-METRO(tmp(end),3);
  end
  file = files(FileNumber);
  eval(['load ' file.name ';'])
  DEEP = x2(linee,:);
  deep(subindx) = DEEP(subindx);
  set_parameters(deep);
  [dr,info] = resol(ys,0);
  if ~info(1)
    ys = dr.ys;
    WelfDistribution(i) = ys(kk(windx))+.5*dr.ghs2(windx);
  else
    WelfDistribution(i) = Inf;
    compt = compt+1;
  end
  waitbar(i/B,hfid);
end
close(hfid);
indx = find(~isinf(WelfDistribution));
WelfareDistribution = WelfDistribution(indx);
oo_.Welfare.PosteriorDraws = WelfareDistribution;
oo_.Welfare.PosteriorMean = mean(WelfareDistribution);
oo_.Welfare.posteriorStd = std(WelfareDistribution);
% [3] Non parametric estimation of the posterior density.
% [3.1] Set some parameters:
number_of_grid_points = 2^9;      % 2^9 = 512 !... Must be a power of two.
bandwidth = 0;                    % Rule of thumb optimal bandwidth parameter.
kernel_function = 'gaussian';     % Gaussian kernel for Fast Fourrier Transform approximaton.
% [3.2] Compute the optimal bandwidth parameter:
optimal_bandwidth = mh_optimal_bandwidth(WelfareDistribution,length(indx),bandwidth,kernel_function); 
% [3.3] Estimation of the posterior density:
[x1,f1] = kernel_density_estimate(WelfareDistribution,number_of_grid_points,...
    optimal_bandwidth,kernel_function);
oo_.Welfare.PosteriorDensity.abscissa = x1;
oo_.Welfare.PosteriorDensity.density = f1;
disp(['Percentage of mh-draws violating the B&K conditions: ' num2str(compt/B*100)  '%.'])
disp(' ')
% [4] Plot.
if pltOpt
  figure('Name','Posterior distribution of the Welfare.')
  plot(x1,f1,'-k','linewidth',2)
  axis tight
  box on
  xlabel('Welfare')
  ylabel('Posterior density')
end
options_.order = 1;