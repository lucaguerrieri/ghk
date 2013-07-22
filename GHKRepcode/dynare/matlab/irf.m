function y_=irf(dr, e1, long_, drop_, replic, iorder)
  global lgy_ ykmin_ endo_nbr exo_nbr exe_ ykmin_ xkmin_ Sigma_e_ iter_ lgx_

  old_iter = iter_;
  iter_ = long_;
  
  temps = repmat(dr.ys,1,ykmin_);

  y_	= 0;
  if iorder == 1
    iter_ = long_;
    y1_ = repmat(dr.ys,1,iter_);
    ex2_ = zeros(iter_,exo_nbr);
    ex2_(1,:) = e1';
    y2_ = simult_(repmat(dr.ys,1,ykmin_),dr,ex2_,iorder);
    y_ = y2_(:,ykmin_+1:end)-y1_;
  else
    % eliminate shocks with 0 variance
    i_exo_var = setdiff([1:exo_nbr],find(diag(Sigma_e_) == 0));
    nxs = length(i_exo_var);
    ex1_ = zeros(long_+drop_+xkmin_,exo_nbr);
    ex2_ = ex1_;
    chol_S = chol(Sigma_e_(i_exo_var,i_exo_var));

    for j = 1: replic
      randn('seed',j);
      ex1_(:,i_exo_var) = randn(long_+drop_+xkmin_,nxs)*chol_S;
      ex2_ = ex1_;
      ex2_(drop_+1,:) = ex2_(drop_+1,:)+e1';   
      y1_ = simult_(repmat(dr.ys,1,ykmin_),dr,ex1_,iorder);
      y2_ = simult_(repmat(dr.ys,1,ykmin_),dr,ex2_,iorder);
      y_ = y_+(y2_(:,ykmin_+drop_+1:end)-y1_(:,ykmin_+drop_+1:end));
    end
    y_=y_/replic;
  end
  iter_ = old_iter;

% 01/18/02 MJ corrected for many lags
% 03/11/22 MJ input is now entire shock vector e1 (for orthogonalized
%              IRFs)


