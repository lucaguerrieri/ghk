% Copyright (C) 2001 Michel Juillard
%
function [x,f,fvec,check]=lnsrch(xold,fold,g,p,stpmax,func,j1,j2,varargin)

  x = [];
  f = [];
  fvec = [];
  alf = 1e-4 ;
  tolx = 3.7e-11 ;
  alam = 1;
  
  x = xold;
  nn = length(j2);
  summ = sqrt(sum(p.*p)) ;
  if ~isfinite(summ)
    warning(['Some element of Newton direction isn''t finite. Jacobian maybe' ...
	   ' singular or there is a problem with initial values'])
    check = 1;
    return
  end
  
  if summ > stpmax
    p=p.*stpmax/summ ;
  end

  slope = g'*p ;
  
  test = max(abs(p)'./max([abs(xold(j2))';ones(1,nn)])) ;
  alamin = tolx/test ;

  if alamin > 0.1
    alamin = 0.1;
  end
  
  while 1
    if alam < alamin
      check = 1 ;
      return
    end
    
    x(j2) = xold(j2) + (alam*p) ;
    fvec = feval(func,x,varargin{:}) ;
    fvec = fvec(j1);
    f = 0.5*fvec'*fvec ;

    if any(isnan(fvec))
      alam = alam/2 ;
      alam2 = alam ;
      f2 = f ;
      fold2 = fold ;
    else

      if f <= fold+alf*alam*slope
	check = 0;
	break ;
      else
	if alam == 1
	  tmplam = -slope/(2*(f-fold-slope)) ;
	else
	  rhs1 = f-fold-alam*slope ;
	  rhs2 = f2-fold2-alam2*slope ;
	  a = (rhs1/(alam^2)-rhs2/(alam2^2))/(alam-alam2) ;
	  b = (-alam2*rhs1/(alam^2)+alam*rhs2/(alam2^2))/(alam-alam2) ;
	  if a == 0
	    tmplam = -slope/(2*b) ;
	  else
	    disc = (b^2)-3*a*slope ;

	    if disc < 0
	      warning ('Roundoff problem in nlsearch') ;
	      check = 1;
	      return
	    else
	      tmplam = (-b+sqrt(disc))/(3*a) ;
	    end

	  end

	  if tmplam > 0.5*alam
	    tmplam = 0.5*alam;
	  end

	end

	alam2 = alam ;
	f2 = f ;
	fold2 = fold ;
	alam = max([tmplam;(0.1*alam)]) ;
      end
    end
  end

% 01/14/01 MJ lnsearch is now a separate function
% 01/12/03 MJ check for finite summ to avoid infinite loop when Jacobian
%             is singular or model is denormalized