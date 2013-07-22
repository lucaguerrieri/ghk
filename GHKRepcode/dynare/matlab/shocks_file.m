% Copyright (C) 2001 Michel Juillard
%
function []=shocks_file(ms_flag)
  global options_ M_ exo_nbr exo_det_nbr lgx_ lgx_det_ ex_ exe_ ex_det_ exe_det_
  
  if isfield(options_,'shocks_file')
    if exo_det_nbr > 0
      x = read_variables(options_.shocks_file,lgx_det_,[]);
      n = size(x,1);
      M_.ex_det_length = n;
      if size(ex_det_,1) >= n
	if ms_flag == 1
	  ex_det_(1:n,:) = ex_det_(1:n,:).*x;
	else
	  ex_det_(1:n,:) = x;
	end
      else
	if ms_flag == 1
	  ex_det_ = ones(n,1)*exe_det_';
	  ex_det_ = ex_det_.*x;
	else
	  ex_det_ = x;
	end
      end
    else
      x = read_variables(options_.shocks_file,lgx_det_);
      n = size(x,1);
      if size(ex_,1) >= n
	if ms_flag == 1
	  ex_(1:n,:) = ex_(1:n,:).*x;
	else
	  ex_(1:n,:) = x;
	end
      else
	if ms_flag == 1
	  ex_ = ones(n,1)*exe_';
	  ex_ = ex_.*x;
	else
	  ex_ = x;
	end
      end
    end
  end
  