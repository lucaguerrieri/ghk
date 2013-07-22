% Copyright (C) 2001 Michel Juillard
%
function [x,check] = solve1(func,x,j1,j2,varargin)

  global gstep_ fjac debug_ unit_root_ lgy_

  nn = length(j1);
  
  fjac = zeros(nn,nn) ;
  g = zeros(nn,1) ;

  tolf = 100*eps^(2/3) ;
  tolmin = 3.7e-11 ;
  tolx = 3.7e-11 ;

  stpmx = 100 ;
  maxit = 2000 ;

  check = 0 ;

  fvec = feval(func,x,varargin{:});
  fvec = fvec(j1);
  
  i = find(~isfinite(fvec));
  
  if ~isempty(i)
    disp(['STEADY:  numerical initial values incompatible with the following' ...
	  ' equations'])
    disp(j1(i)')
  end
  
  f = 0.5*fvec'*fvec ;

  if max(abs(fvec)) < 0.01*tolf
    return ;
  end

  stpmax = stpmx*max([sqrt(x'*x);nn]) ;
  first_time = 1;
  for its = 1:maxit

    dh = max(abs(x(j2)),gstep_*ones(nn,1))*eps^(1/3);
    for j = 1:nn
      xdh = x ;
      xdh(j2(j)) = xdh(j2(j))+dh(j) ;
      t = feval(func,xdh,varargin{:});
      fjac(:,j) = (t(j1) - fvec)./dh(j) ;
      g(j) = fvec'*fjac(:,j) ;
    end

    if debug_
      disp(['cond(fjac) ' num2str(cond(fjac))])
    end
    
    if unit_root_
      if first_time
	first_time = 0;
	[q,r,e]=qr(fjac);
	n = sum(abs(diag(r)) < 1e-12);
	fvec = q'*fvec;
	p = e*[-r(1:end-n,1:end-n)\fvec(1:end-n);zeros(n,1)];
	disp(' ')
	disp('STEADY with unit roots:')
	disp(' ')
	if n > 0
	  disp(['   The following variable(s) kept their value given in INITVAL' ...
		' or ENDVAL'])
	  disp(char(e(:,end-n+1:end)'*lgy_))
	else
	  disp('   STEADY can''t find any unit root!')
	end
      else
	[q,r]=qr(fjac*e);
	fvec = q'*fvec;
	p = e*[-r(1:end-n,1:end-n)\fvec(1:end-n);zeros(n,1)];
      end	
%    elseif cond(fjac) > 10*sqrt(eps)
    elseif cond(fjac) > 1/sqrt(eps)
	fjac2=fjac'*fjac;
	p=-(fjac2+sqrt(nn*eps)*max(sum(abs(fjac2)))*eye(nn))\(fjac'*fvec);
    else
      p = -fjac\fvec ;
    end
    xold = x ;
    fold = f ;

    [x,f,fvec,check]=lnsrch1(xold,fold,g,p,stpmax,func,j1,j2,varargin{:});

    if debug_
      disp([its f])
      disp([xold x])
    end
      
    if check > 0
      den = max([f;0.5*nn]) ;
      if max(abs(g).*max([abs(x(j2)');ones(1,nn)])')/den < tolmin
	return
      else
	disp (' ')
	disp (['SOLVE: Iteration ' num2str(its)])
	disp (['Spurious convergence.'])
	disp (x)
	return
      end

      if max(abs(x-xold)./max([abs(x);ones(1,nn)])') < tolx
	disp (' ')
	disp (['SOLVE: Iteration ' num2str(its)])
	disp (['Convergence on dX.'])
	disp (x)
	return
      end
    elseif f < tolf
      return
    end
  end
  
  check = 1;
  disp(' ')
  disp('SOLVE: maxit has been reached')

% 01/14/01 MJ lnsearch is now a separate function
% 01/16/01 MJ added varargin to function evaluation
% 04/13/01 MJ added test  f < tolf !!
% 05/11/01 MJ changed tests for 'check' so as to remove 'continue' which is
%             an instruction which appears only in version 6










