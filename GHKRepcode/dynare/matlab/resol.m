% Copyright (C) 2001 Michel Juillard
%
function [dr,info]=resol(ys,check_flag)
% info: same as dr1 
% plus: 
% 11 .... same as dr1 for dr_algo = 2
% 20: can't find steady state info(2) contains sum of sqare residuals
  
global jacobia_ iy_ ykmin_ ykmax_ gstep_ exo_nbr exo_det_nbr endo_nbr
global ex_ ex_det_ valf_ it_ exe_ exe_det_ xkmin_ xkmax_ 
global fname_ means_ stderrs_ lgy_ maxit_
global dynatol_ options_

options_ = set_default_option(options_,'olr',0);
info = 0;

it_ = ykmin_ + 1 ;

if exo_nbr == 0
  exe_ = [] ;
end

% check if ys is steady state
tempex = ex_;
tempexdet = ex_det_;
ex_ = repmat(exe_',xkmin_+xkmax_+1,1);
if exo_det_nbr > 0 
  ex_det_ = ones(ykmin_+1,1)*exe_det_';
end
fh = str2func([fname_ '_fff']);
if max(abs(feval(fh,ys))) > dynatol_ & options_.olr == 0
  if exist([fname_ '_steadystate'])
    [dr.ys,check1] = feval([fname_ '_steadystate'],ys);
  else
    [dr.ys,check1] = dynare_solve([fname_ '_fff'],ys);
  end
  if check1
    info(1) = 20;
    resid = feval(fh,ys);
    info(2) = resid'*resid; % penalty...
    return
  end
else 
  dr.ys = ys;
end

dr.fbias = zeros(endo_nbr,1);
[dr,info] = dr1(dr,check_flag);

if info(1)
  return
end

if options_.dr_algo == 1 & options_.order > 1
  dr.ys = dynare_solve('dr2',ys,dr);
  dr.fbias = 2*feval([fname_ '_fff'],dr.ys);
  [dr, info1] = dr1(dr,check_flag);
  if info1(1)
    info(1) = info(1)+10;
    return
  end
end
ex_det_ = tempexdet;
ex_ = tempex;
tempex = [];

% 01/01/2003 MJ added dr_algo == 1
% 08/24/2001 MJ uses Schmitt-Grohe and Uribe (2001) constant correction
%               in dr.ghs2 
% 05/26/2003 MJ added temporary values for ex_