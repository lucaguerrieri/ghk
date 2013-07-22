% the beginning and the end of this function may be adapted by the userx
function [loss,vx,info]=osr_obj(x,params,weights);
  global ys_ Sigma_e_ endo_nbr exo_nbr optimal_Q_ it_ ykmin_ options_
  
  vx = [];
  % set parameters of the policiy rule
  np = size(params,1);
  for i=1:np
    assignin('base',deblank(params(i,:)),x(i))
  end
  
  % don't change below until the part where the loss function is computed
  it_ = ykmin_+1;
  [dr_,info] = resol(ys_,0);
  
  switch info(1)
   case 1
    loss = 1e8;
    return
   case 2
    loss = 1e8*min(1e3,info(2));
    return
   case 3
    loss = 1e8*min(1e3,info(2));
    return
   case 4
    loss = 1e8*min(1e3,info(2));
    return
   case 5
    loss = 1e8;
    return
   case 20
    loss = 1e8*min(1e3,info(2));
    return
   otherwise
  end
  
  [A,B] = kalman_transition_matrix(dr_);
  [vx, info] = lyapunov_symm(A,B*Sigma_e_*B');
  if info > 0
    loss = 1e8;
    return
  end
  weights = weights(dr_.order_var,dr_.order_var);
  vx = vx(1:endo_nbr,1:endo_nbr);
  loss = weights(:)'*vx(:);
  














