function initial_estimation_checks(xparam1,gend,data)
global bayestopt_ estim_params_ exo_nbr dsge_prior_weight
  
nv = size(data,1);
  
if nv > exo_nbr + estim_params_.nvn
  error(['Estimation can''t take place because there are less shocks than' ...
     'observed variables'])
end
  
r = rank(data);
if r < nv
  error(['Estimation can''t take place because the data are perfectly' ...
     ' correlated']);
end

if ~isempty(strmatch('dsge_prior_weight',estim_params_.param_names)) | ~isempty(dsge_prior_weight)
  [fval,cost_flag,ys,trend_coeff,info] = DsgeVarLikelihood(xparam1,gend);
else
  [fval,cost_flag,ys,trend_coeff,info] = DsgeLikelihood(xparam1,gend,data);
end

disp(['Initial value of the posterior (or likelihood): ' num2str(fval)]);

if info(1) > 0
  disp('Error in computing likelihood for initial parameter values')
  print_info(info)
end
