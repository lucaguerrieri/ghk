% Copyright (C) 2001 Michel Juillard
%
function resid0
  global iter_ valf_ ex_ y_ it_ exe_ ys_ iy_ ykmin_ ykmax_ endval_ z
  global fname_
  
  n = size(iy_,2);
%  if ~ valf_ | size(y_,2) ~= iter_+ykmin_+ykmax_
  if ~ valf_ 
    if size(ys_,1) == 1 & ys_ == 0
      ys_ = zeros(size(ys_,1),1) ;
    end
    y_ = ys_*ones(1,iter_+ykmin_+ykmax_) ;
    if endval_ == 1
      y_(:,1:ykmin_) = ys0_*ones(1,ykmin_) ;
    end
  end

  i = iy_';
  iyr0 = find(i(:));

  y =y_(:);
  it_ = ykmin_+1;
  fh = str2func([fname_ '_fff']);
  z = feval(fh,ys_);


  disp(z'); 






