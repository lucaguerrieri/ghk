% Copyright (C) 2001 Michel Juillard
%
function steady_()

  global it_ ex_ exe_ iy_ fname_ ex_det_ exe_det_
  global ykmin_ ykmax_ ys_ xkmin_ xkmax_ exo_det_nbr

  x = ys_ ;
  xlen = xkmin_ + xkmax_ + 1 ;
  nn = size(iy_,2) ;
  it_ = ykmin_+1 ;
  temp = ex_ ;
  ex_ = ones(xlen,1)*exe_' ;

  if exo_det_nbr > 0
    tempdet = ex_det_ ;
    ex_det_ = ones(ykmin_+1,1)*exe_det_' ;
  end

  if exist([fname_ '_steadystate'])
    [ys_,check] = feval([fname_ '_steadystate'],x);
  else
    [ys_,check] = dynare_solve([fname_ '_fff'],x);
  end

  if check ~= 0
    error('STEADY: convergence problems')
  end

  if exo_det_nbr > 0
    ex_det_ = tempdet;
  end
  ex_ = temp ;

% 06/24/01 MJ: steady_ no results printer; steady with printed results
% 07/31/03 MJ: in case of convergence problem steady stops with an error







