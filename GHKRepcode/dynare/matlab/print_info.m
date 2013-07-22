% Copyright (C) 2005 Michel Juillard
%
function print_info(info)
  global options_

  options_ = set_default_option(options_,'noprint',0);
  
  if ~options_.noprint
    if info(1) > 10 & info(1) < 20
      disp('Failure in dr_algo=2')
      info(1) = info(1) - 10;
    end
    switch info(1)
     case 1
      error(['The model doesn''t determine the current variables' ...
	     ' uniquely'])
     case 2
      error(['MJDGGES returns the following error code' ...
	     int2str(info(2))])
     case 3
      error(['Blanchard Kahn conditions are not satisfied: no stable' ...
	     ' equilibrium'])
     case 4
      error(['Blanchard Kahn conditions are not satisfied:' ...
	     ' indeterminacy'])
     case 5
      error(['Blanchard Kahn conditions are not satisfied:' ...
	     ' indeterminacy due to rank failure'])
     case 20
      error(['Impossible to find the steady state. Either the model' ...
	     ' doesn''t have a unique steady state of the guess values' ...
	     ' are too far from the solution']) 
     case 30
      error('Variance can''t be computed')
     otherwise
      error('This case shouldn''t happen. Contact the authors of Dynare')
    end
  end
  