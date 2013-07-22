% Copyright (C) 2003 Michel Juillard
%
% set
function set_shocks(flag,k,ivar,values)
  global ex_ exe_ ex_det_ exe_det_
  
  n1 = size(ex_,1);
  n2 = size(ex_det_,1);
  if k(end) > n1 
    if flag <= 1
      ex_ = [ex_; ones(k(end)-n1,1)*exe_'];
    else
      ex_det_ = [ex_det_; ones(k(end)-n2,1)*exe_det_'];
    end
  end  
  
  if flag == 0
    ex_(k,ivar) = ones(length(k),1).*values;
  elseif flag == 1
    ex_(k,ivar) = ex_(k,ivar).*values;
  elseif flag == 2
    ex_det_(k,ivar) = ones(length(k),1).*values;
  elseif flag == 3
    ex_det_(k,ivar) = ex_det_(k,ivar).*values;
  end

  % 05/29/03 MJ