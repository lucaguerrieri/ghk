% Copyright (C) 2001 Michel Juillard
%
function make_ex_
  global exe_ ex_ ex0_ xkmin_ xkmax_ exo_nbr options_ exe_det_ ex_det_ ...
      exo_det_nbr ex_det0_ ykmin_ ykmax_
  
  options_ = set_default_option(options_,'periods',0);
  
  if isempty(exe_)
    exe_ = zeros(exo_nbr,1);
  end
  
  if exo_det_nbr > 1 & isempty(exe_det_)
    exe_det_ = zeros(exo_det_nbr,1);
  end
  
  if isempty(ex_)
    if isempty(ex0_)
      ex_ = [ones(xkmin_+options_.periods+xkmax_,1)*exe_'];
    else
      ex_ = [ones(xkmin_,1)*ex0_';ones(options_.periods+xkmax_,1)*exe_'];
    end
  elseif size(ex_,2) < length(exe_)
    k = size(ex_,2)+1:length(exe_)
    if isempty(ex0_)
      ex_ = [ex_ ones(xkmin_+size(ex_,1)+xkmax_,1)*exe_(k)'];
    else
      ex_ = [ex_ [ones(xkmin_,1)*ex0_(k)'; ones(size(ex_,1)-xkmin_+xkmax_, ...
						1)*exe_(k)']];
    end
  elseif size(ex_,1) < xkmin_+xkmax_+options_.periods
    if isempty(ex0_)
      ex_ = [ex_; ones(xkmin_+options_.periods+xkmax_-size(ex_,1),1)*exe_'];
    else
      ex_ = [ones(xkmin_,1)*ex0_'; ex_; ones(options_.periods+xkmax_-size(ex_, ...
						  1),1)*exe_'];
    end
  end
  
  if exo_det_nbr > 0
    if isempty(ex_det_)
      if isempty(ex_det0_)
	ex_det_ = [ones(ykmin_+options_.periods+ykmax_,1)*exe_det_'];
      else
	ex_det_ = [ones(ykmin_,1)*ex_det0_';ones(options_.periods+ykmax_,1)*exe_det_'];
      end
    elseif size(ex_det_,2) < length(exe_det_)
      k = size(ex_det_,2)+1:length(exe_det_)
      if isempty(ex_det0_)
	ex_det_ = [ex_det_ ones(ykmin_+size(ex_det_,1)+ykmax_,1)*exe_det_(k)'];
      else
	ex_det_ = [ex_det_ [ones(ykmin_,1)*ex_det0_(k)'; ones(size(ex_det_,1)-ykmin_+ykmax_, ...
						  1)*exe_det_(k)']];
      end
    elseif size(ex_det_,1) < ykmin_+ykmax_+options_.periods
      if isempty(ex_det0_)
	ex_det_ = [ex_det_; ones(ykmin_+options_.periods+ykmax_-size(ex_det_,1),1)*exe_det_'];
      else
	ex_det_ = [ones(ykmin_,1)*ex_det0_'; ex_det_; ones(options_.periods+ykmax_-size(ex_det_, ...
						  1),1)*exe_det_'];
      end
    end
  end
	     