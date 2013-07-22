% targets and iy order: 1) variances 2) correlations 
% 3) constraints on Sigma_e_ itself 4) autocorrelations
function f=calib_obj(Sigma_e_,A,ghu1,ghx,ghu,targets,var_weights,iy,nar)
  global vx fold
  
  Gamma_y = cell(nar+1,1);
%  Sigma_e_ = Sigma_e_'*Sigma_e_;
  Sigma_e_=diag(Sigma_e_);
  nx = size(ghx,2);
  b=ghu1*Sigma_e_*ghu1';
  vx = [];
  if isempty(vx)
    vx = lyapunov_symm(A,b);
  else
    [vx,status] = bicgstab(@f_var,b(:),vx(:),1e-8,50,A,nx);
    if status
      vx = lyapunov_symm(A,b);
    else
      vx=reshape(vx,nx,nx);
    end
  end
  Gamma_y{1} = ghx*vx*ghx'+ ghu*Sigma_e_*ghu';
  f = 0;
  if ~isempty(targets{1})
    e = targets{1}-sqrt(Gamma_y{1}(iy{1}));
    f = e'*(var_weights{1}.*e);
  end

  sy = sqrt(diag(Gamma_y{1}));
  sy = sy *sy';
  if ~isempty(targets{2})
    e = targets{2}-Gamma_y{1}(iy{2})./(sy(iy{2})+1e-10);
    f = f+e'*(var_weights{2}.*e);
  end
  
  if ~isempty(targets{3})
    e = targets{3}-sqrt(Sigma_e_(iy{3}));
    f = f+e'*(var_weights{3}.*e);
  end
  
  % autocorrelations
  if nar > 0
    vxy = (A*vx*ghx'+ghu1*Sigma_e_*ghu');
    
    Gamma_y{2} = ghx*vxy./(sy+1e-10);
    if ~isempty(targets{4})
      e = targets{4}-Gamma_y{2}(iy{4});
      f = f+e'*(var_weights{4}.*e);
    end
    
    for i=2:nar
      vxy = A*vxy;
      Gamma_y{i+1} = ghx*vxy./(sy+1e-10);
      if ~isempty(targets{i+3})
	e = targets{i+3}-Gamma_y{i+1}(iy{i+3});
	f = f+e'*(var_weights{i+3}.*e);
      end
    end
  end
  if isempty(fold) | f < 2*fold
    fold = f;
    vxold = vx;
  end
  % 11/04/02 MJ generalized for correlations, autocorrelations and
  %             constraints on Sigma_e_
  % 01/25/03 MJ targets std. dev. instead of variances
  