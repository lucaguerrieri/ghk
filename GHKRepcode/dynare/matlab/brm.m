% Copyright (C) 2001 Michel Juillard
%
function [e_1, e_2, e_inf]=brm(v,nbr,drop)
  global y_ iy_ iter_ ykmin_ ykmax_ it_ endo_nbr ex_ fname_
  
  i=(iy_ > 0)' ;
  iyr0 = find(i(:)) ;

  ex_old = ex_;
  ex_ = ex_(drop+1:end,:);
  y = y_(:,drop+1:end);
  y = y(:);
  
  iter = iter_ - drop - ykmin_ - ykmax_;

  z = zeros(endo_nbr,iter);
  fh = str2func([fname_ '_ff']);
  for it_ = ykmin_+1 : iter+ykmin_
    z(:,it_) = feval(fh,y(iyr0));
    iyr0 = iyr0 + endo_nbr;
  end
  
    t1 = z(nbr,:)';
    t2 = v(drop+1:end-2);
    t = t1./t2;
  e_1 = log10(mean(abs(t)));
  e_2 = log10(var(t));
  e_inf = log10(max(abs(t)));

  ex_ = ex_old;
