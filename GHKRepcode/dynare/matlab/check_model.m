function check_model()
  global iy_ ykmin_ ykmax_ xkmin_ xkmax_ iy_ exo_det_nbr
  
  xlen = xkmin_ + xkmax_ + 1;
  if ~ iy_(ykmin_+1,:) > 0
  error ('RESOL: Error in model specification: some variables don"t appear as current') ;
end

if xlen > 1
  error (['RESOL: stochastic exogenous variables must appear only at the' ...
	  ' current period. Use additional endogenous variables']) ;
end

if (exo_det_nbr > 0) & (ykmin_ > 1 | ykmax_ > 1)
  error(['Exogenous deterministic variables are currently only allowed in' ...
	 ' models with leads and lags on only one period'])
end

