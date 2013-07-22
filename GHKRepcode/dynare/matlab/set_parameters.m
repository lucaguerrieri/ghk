function set_parameters(xparam1)
  global estim_params_ Sigma_e_
  
  for i=1:estim_params_.nvx
    k =estim_params_.var_exo(i,1);
    Sigma_e_(k,k) = xparam1(i)*xparam1(i);
  end
  offset = estim_params_.nvx+estim_params_.nvn;
  for i=1:estim_params_.ncx
    k1 =estim_params_.corrx(i,1);
    k2 =estim_params_.corrx(i,2);
    Sigma_e_(k1,k2) = xparam1(i+offset)*sqrt(Sigma_e_(k1,k1)*Sigma_e_(k2,k2));
    Sigma_e_(k2,k1) = Sigma_e_(k1,k2);
  end
  offset = offset+estim_params_.ncx+estim_params_.ncn;
  for i=1:estim_params_.np
    assignin('base',deblank(estim_params_.param_names(i,:)),xparam1(i+offset));
  end
