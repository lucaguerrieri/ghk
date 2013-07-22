function [A,B,ys,info] = dynare_resolve()
  global ys_ dr_
  
  [dr_,info] = resol(ys_,0);
  
  if info(1) > 0
    A = [];
    B = [];
    ys = [];
    return
  end
  
  [A,B] = kalman_transition_matrix(dr_);
  ys = dr_.ys;
