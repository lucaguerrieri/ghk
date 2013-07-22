% targets and iy order: 1) variances 2) correlations 
% 3) constraints on Sigma_e_ itself 4) autocorrelations
function objective=calib_obj2(Sigma_e_,A,ghu1,ghx,ghu,targets,var_weights,iy,nar)
  global vx fold
  
  objective = cell (nar+3);
  Gamma_y = cell(nar+1,1);
  Sigma_e_=diag(Sigma_e_);
  nx = size(ghx,2);
  b=ghu1*Sigma_e_*ghu1';
  vx = lyapunov_symm(A,b);
  Gamma_y{1} = ghx*vx*ghx'+ ghu*Sigma_e_*ghu';
  if ~isempty(targets{1})
    objective{1} = sqrt(Gamma_y{1}(iy{1}));
  end

  sy = sqrt(diag(Gamma_y{1}));
  sy = sy *sy';
  if ~isempty(targets{2})
    objective{2} = Gamma_y{1}(iy{2})./(sy(iy{2})+1e-10);
  end
  
  if ~isempty(targets{3})
    objective{3} = Sigma_e_(iy{3});
  end
  
  % autocorrelations
  if nar > 0
    vxy = (A*vx*ghx'+ghu1*Sigma_e_*ghu');
    
    Gamma_y{2} = ghx*vxy./(sy+1e-10);
    if ~isempty(targets{4})
      objective{4} = Gamma_y{2}(iy{4});
    end
    
    for i=2:nar
      vxy = A*vxy;
      Gamma_y{i+1} = ghx*vxy./(sy+1e-10);
      if ~isempty(targets{i+3})
	objecitve{i+3} = Gamma_y{i+1}(iy{i+3});
      end
    end
  end

  % 11/04/02 MJ generalized for correlations, autocorrelations and
  %             constraints on Sigma_e_