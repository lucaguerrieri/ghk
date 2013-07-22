% Copyright (C) 2001 Michel Juillard
%
function resid(iter_)
  global valf_ ex_ y_ it_ exe_ ys_ ys0_ iy_ ykmin_ ykmax_ endval_ z
  global fname_ xkmin_ xkmax_
  
  ex_ = ones(xkmin_+xkmax_+iter_,1)*exe_';
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
  z = zeros(n,iter_);
  fh = str2func([fname_ '_ff']);
  for it_=ykmin_+1:iter_+ykmin_
    z(:,it_-ykmin_) = feval(fh,y(iyr0));
    iyr0 = iyr0 + n;
  end

  disp([[1:iter_]' z']); 






