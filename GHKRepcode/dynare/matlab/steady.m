% Copyright (C) 2001 Michel Juillard
%
function steady(linear)

  global ys_ lgy_ ys0_ ykmin_ ykmax_ y_ options_

  steady_;
  
  disp(' ')
  disp('STEADY-STATE RESULTS:')
  disp(' ')
  for i=1:size(ys_,1)
    disp(sprintf('%s \t\t %g',lgy_(i,:),ys_(i)));
  end
  
% overwrites the initialization of y_ in case it was
% set by initval
  if isempty(ys0_)
    y_(:,1:ykmin_) = ys_ * ones(1,ykmin_);
  else
    options_ =set_default_option(options_,'periods',1);
    y_(:,ykmin_+1:ykmin_+options_.periods+ykmax_) = ys_ * ones(1,options_.periods+ykmax_);
  end
  
% 06/24/01 MJ steady print results; steady_ doesn't
% 09/22/01 FC corrected lgy(i,:)
% 05/29/03 MJ sets initial values of y_

