% Copyright (C) 2001 Michel Juillard
%
function [x,check] = dynare_solve(func,x,varargin)
  global gstep_ options_ debug_
  

% unfinished
  jacobian_flag = 0;   

  options_ = set_default_option(options_,'solve_algo',2);
  check = 0;
  func = str2func(func);
  if options_.solve_algo == 0
    if ~isempty(which('fsolve')) & sscanf(version('-release'),'%d') >= 13;
      options=optimset('fsolve');
      options.MaxFunEvals = 20000;
      options.TolFun=1e-8;
      options.Display = 'off';
      [x,fval,exitval,output] = fsolve(func,x,options,varargin{:});
      if exitval > 0
	check = 0;
      else
	check = 1;
      end
      return
    else 
      options_.solve_algo = 1;
    end
  end

  if options_.solve_algo == 1
    nn = size(x,1) ;
    [x,check]=solve1(func,x,1:nn,1:nn,jacobian_flag,varargin{:});
  elseif options_.solve_algo == 2
    nn = size(x,1) ;
    %    tolf = eps^(2/3) ;
    tolf = 1e-9;

    fjac = zeros(nn,nn) ;

    fvec = feval(func,x,varargin{:});

    i = find(~isfinite(fvec));
    
    if ~isempty(i)
      if debug_
	disp(['STEADY:  numerical initial values incompatible with the following' ...
	      ' equations'])
	disp(i')
	error('exiting ...')
      else
	check = 1;
	return
      end
    end
    
    f = 0.5*fvec'*fvec ;

    if max(abs(fvec)) < 0.01*tolf
      return ;
    end

    dh = max(abs(x),gstep_*ones(nn,1))*eps^(1/3);
    for j = 1:nn
      xdh = x ;
      xdh(j) = xdh(j)+dh(j) ;
      fjac(:,j) = (feval(func,xdh,varargin{:}) - fvec)./dh(j) ;
    end

    [j1,j2,r,s] = dmperm(fjac);
    
    for i=length(r)-1:-1:1
      [x,check]=solve1(func,x,j1(r(i):r(i+1)-1),j2(r(i):r(i+1)-1),varargin{:});
      if check
	if debug_
	  error(sprintf('Solve block = %d check = %d\n',i,check));
	else
	  return
	end
      end
    end
    [x,check]=solve1(func,x,1:nn,1:nn,varargin{:});
      
  end
%    fvec1 = feval(func,x,varargin{:})

  % 08/28/03 MJ add a final call to solve1 for solve_algo == 1 in case
  %             initvals generates 'false' zeros in the Jacobian
  
