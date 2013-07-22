function xparam=get_posterior_parameters(type)
  global estim_params_ oo_ Sigma_e_ lgx_ options_
  
  nvx = estim_params_.nvx;
  nvn = estim_params_.nvn;
  ncx = estim_params_.ncx;
  ncn = estim_params_.ncn;
  np  = estim_params_.np;
  
  xparam = zeros(nvx+nvn+ncx+ncn+np,1);
  
  m = 1;
  for i=1:nvx
    k1 = estim_params_.var_exo(i,1);
    name1 = deblank(lgx_(k1,:));
    xparam(m) = eval(['oo_.' type '.shocks_std.' name1]);
    Sigma_e_(k1,k1) = xparam(m)^2;
    m = m+1;
  end
  
  for i=1:nvn
    k1 = estim_params_.var_endo(i,1);
    name1 = deblank(options_.varobs(k1,:));
    xparam(m) = eval(['oo_.' type '.measurement_errors_std.' name1]);
    m = m+1;
  end
  
  for i=1:ncx
    k1 = estim_params_.corrx(i,1);
    k2 = estim_params_.corrx(i,2);
    name1 = deblank(lgx_(k1,:));
    name2 = deblank(lgx_(k2,:));
    xparam(m) = eval(['oo_.' type '.shocks_corr.' name1 '_' name2]);
    Sigma_e_(k1,k2) = xparam(m);
    Sigma_e_(k2,k1) = xparam(m);
    m = m+1;
  end
  
  for i=1:ncn
    k1 = estim_params_.corrn(i,1);
    k2 = estim_params_.corrn(i,2);
    name1 = deblank(options_.varobs(k1,:));
    name2 = deblank(options_.varobs(k2,:));
    xparam(m) = eval(['oo_.' type '.measurement_errors_corr.' name1 '_' name2]);
    m = m+1;
  end

  for i=1:np
    name1 = deblank(estim_params_.param_names(i,:));
    xparam(m) = eval(['oo_.' type '.parameters.' name1]);
    assignin('base',name1,xparam(m));
    m = m+1;
  end
  

  