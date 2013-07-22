% Copyright (C) 2001 Michel Juillard
%
function make_y_
  global ys_ y_ ys0_ ykmin_ ykmax_ endo_nbr options_
  
  options_ = set_default_option(options_,'periods',0);
  
  if isempty(ys_)
    ys_ = ones(endo_nbr,1);
  end
  
  
  if isempty(y_)
    if isempty(ys0_)
      y_ = [ys_*ones(1,ykmin_+options_.periods+ykmax_)];
    else
      y_ = [ys0_*ones(1,ykmin_);ys_*ones(1,options_.periods+ykmax_)];
    end
  elseif size(y_,1) < length(ys_)
    k = size(y_,1)+1:length(ys_)
    if isempty(ys0_)
      y_ = [y_; ys_(k)*ones(1,ykmin_+size(y_,1)+ykmax_)];
    else
      y_ = [y_; [ys0_(k)*ones(1,ykmin_); ys_(k)*ones(1,size(y_,2)-ykmin_+ ...
						     ykmax_)]];
    end
  elseif size(y_,2) < ykmin_+ykmax_+options_.periods
    if isempty(ys0_)
      y_ = [y_ ys_*ones(1,ykmin_+options_.periods+ykmax_-size(y_,2),1)];
    else
      y_ = [ys0_*ones(1,ykmin_) y_  ys_*ones(1,options_.periods+ykmax_-size(y_, ...
						  2))];
    end
  end
    
	     