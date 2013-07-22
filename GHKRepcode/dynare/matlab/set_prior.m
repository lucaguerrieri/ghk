function [xparam1,estim_params_,bayestopt_,lb,ub]=set_prior(estim_params_)
  global lgx_ lgy_ options_
  
  nvx = size(estim_params_.var_exo,1);
  nvn = size(estim_params_.var_endo,1);
  ncx = size(estim_params_.corrx,1);
  ncn = size(estim_params_.corrn,1);
  np = size(estim_params_.param_vals,1);
  
  estim_params_.nvx = nvx;
  estim_params_.nvn = nvn;
  estim_params_.ncx = ncx;
  estim_params_.ncn = ncn;
  estim_params_.np = np;
  
  xparam1 = [];
  ub = [];
  lb = [];
  bayestopt_.pshape = [];
  bayestopt_.pmean = [];
  bayestopt_.pstdev = [];
  bayestopt_.p1 = [];
  bayestopt_.p2 = [];
  bayestopt_.p3 = [];
  bayestopt_.p4 = [];
  bayestopt_.jscale = [];
  bayestopt_.name = [];
  if nvx
    xparam1 = estim_params_.var_exo(:,2);
    ub = estim_params_.var_exo(:,4); 
    lb = estim_params_.var_exo(:,3); 
    bayestopt_.pshape =  estim_params_.var_exo(:,5);
    bayestopt_.pmean =  estim_params_.var_exo(:,6);
    bayestopt_.pstdev =  estim_params_.var_exo(:,7);
    bayestopt_.p3 =  estim_params_.var_exo(:,8);
    bayestopt_.p4 =  estim_params_.var_exo(:,9);
    bayestopt_.jscale =  estim_params_.var_exo(:,10);
    bayestopt_.name = cellstr(lgx_(estim_params_.var_exo(:,1),:));
  end
  if nvn
    for i=1:nvn
      estim_params_.var_endo(i,1) = strmatch(deblank(lgy_(estim_params_.var_endo(i,1),:)),deblank(options_.varobs),'exact');
    end
    xparam1 = [xparam1; estim_params_.var_endo(:,2)];
    ub = [ub; estim_params_.var_endo(:,4)]; 
    lb = [lb; estim_params_.var_endo(:,3)]; 
    bayestopt_.pshape = [ bayestopt_.pshape; estim_params_.var_endo(:,5)];
    bayestopt_.pmean = [ bayestopt_.pmean; estim_params_.var_endo(:,6)];
    bayestopt_.pstdev = [ bayestopt_.pstdev; estim_params_.var_endo(:,7)];
    bayestopt_.p3 = [ bayestopt_.p3; estim_params_.var_endo(:,8)];
    bayestopt_.p4 = [ bayestopt_.p4; estim_params_.var_endo(:,9)];
    bayestopt_.jscale = [ bayestopt_.jscale; estim_params_.var_endo(:,10)];
    bayestopt_.name = cellstr(strvcat(char(bayestopt_.name),...
				      lgy_(estim_params_.var_endo(:,1),:)));
  end
  if ncx
    xparam1 = [xparam1; estim_params_.corrx(:,3)];
    ub = [ub; max(min(estim_params_.corrx(:,5),1),-1)];
    lb = [lb; max(min(estim_params_.corrx(:,4),1),-1)];
    bayestopt_.pshape = [ bayestopt_.pshape; estim_params_.corrx(:,6)];
    bayestopt_.pmean = [ bayestopt_.pmean; estim_params_.corrx(:,7)];
    bayestopt_.pstdev = [ bayestopt_.pstdev; estim_params_.corrx(:,8)];
    bayestopt_.p3 = [ bayestopt_.p3; estim_params_.corrx(:,9)];
    bayestopt_.p4 = [ bayestopt_.p4; estim_params_.corrx(:,10)];
    bayestopt_.jscale = [ bayestopt_.jscale; estim_params_.corrx(:,11)];
    bayestopt_.name = cellstr(strvcat(char(bayestopt_.name),...
				      char(strcat(cellstr(lgx_(estim_params_.corrx(:,1),:)),...
						  ',',...
						  cellstr(lgx_(estim_params_.corrx(:,2),:))))));
  end
  if ncn
    xparam1 = [xparam1; estim_params_.corrn(:,3)];
    ub = [ub; max(min(estim_params_.corrn(:,5),1),-1)];
    lb = [lb; max(min(estim_params_.corrn(:,4),1),-1)];
    bayestopt_.pshape = [ bayestopt_.pshape; estim_params_.corrn(:,6)];
    bayestopt_.pmean = [ bayestopt_.pmean; estim_params_.corrn(:,7)];
    bayestopt_.pstdev = [ bayestopt_.pstdev; estim_params_.corrn(:,8)];
    bayestopt_.p3 = [ bayestopt_.p3; estim_params_.corrn(:,9)];
    bayestopt_.p4 = [ bayestopt_.p4; estim_params_.corrn(:,10)];
    bayestopt_.jscale = [ bayestopt_.jscale; estim_params_.corrn(:,11)];
    bayestopt_.name = cellstr(strvcat(char(bayestopt_.name),...
				      char(strcat(cellstr(lgy_(estim_params_.corrn(:,1),:)),...
						  ',',...
						  cellstr(lgy_(estim_params_.corrn(:,2),:))))));
  end
  if np
    xparam1 = [xparam1; estim_params_.param_vals(:,1)];
    ub = [ub; estim_params_.param_vals(:,3)];
    lb = [lb; estim_params_.param_vals(:,2)];
    bayestopt_.pshape = [ bayestopt_.pshape; estim_params_.param_vals(:,4)];
    bayestopt_.pmean = [ bayestopt_.pmean; estim_params_.param_vals(:,5)];
    bayestopt_.pstdev = [ bayestopt_.pstdev; estim_params_.param_vals(:,6)];
    bayestopt_.p3 = [ bayestopt_.p3; estim_params_.param_vals(:,7)];
    bayestopt_.p4 = [ bayestopt_.p4; estim_params_.param_vals(:,8)];
    bayestopt_.jscale = [ bayestopt_.jscale; estim_params_.param_vals(:, ...
						  9)];
    bayestopt_.name = cellstr(strvcat(char(bayestopt_.name),estim_params_.param_names));
  end

  bayestopt_.ub = ub;
  bayestopt_.lb = lb;
  
  bayestopt_.p1 = bayestopt_.pmean;
  bayestopt_.p2 = bayestopt_.pstdev;
  
  % generalized location parameters by default for beta distribution
  k = find(bayestopt_.pshape == 1);
  k1 = find(isnan(bayestopt_.p3(k)));
  bayestopt_.p3(k(k1)) = zeros(length(k1),1);
  k1 = find(isnan(bayestopt_.p4(k)));
  bayestopt_.p4(k(k1)) = ones(length(k1),1);
  
  % generalized location parameter by default for gamma distribution
  k = find(bayestopt_.pshape == 2);
  k1 = find(isnan(bayestopt_.p3(k)));
  bayestopt_.p3(k(k1)) = zeros(length(k1),1);
  
  % truncation parameters by default for normal distribution
  k = find(bayestopt_.pshape == 3);
  k1 = find(isnan(bayestopt_.p3(k)));
  bayestopt_.p3(k(k1)) = -Inf*ones(length(k1),1);
  k1 = find(isnan(bayestopt_.p4(k)));
  bayestopt_.p4(k(k1)) = Inf*ones(length(k1),1);

  k = find(bayestopt_.pshape == 4);
  for i=1:length(k)
    [bayestopt_.p1(k(i)),bayestopt_.p2(k(i))] = ...
	inverse_gamma_specification(bayestopt_.pmean(k(i)),bayestopt_.pstdev(k(i)),1);
  end
  
  k = find(bayestopt_.pshape == 5);
  for i=1:length(k)
    [bayestopt_.pmean(k(i)),bayestopt_.pstdev(k(i)),bayestopt_.p1(k(i)),bayestopt_.p2(k(i))] = ...
	uniform_specification(bayestopt_.pmean(k(i)),bayestopt_.pstdev(k(i)),bayestopt_.p3(k(i)),bayestopt_.p4(k(i)));
  end
  
  k = find(bayestopt_.pshape == 6);
  for i=1:length(k)
    [bayestopt_.p1(k(i)),bayestopt_.p2(k(i))] = ...
	inverse_gamma_specification(bayestopt_.pmean(k(i)),bayestopt_.pstdev(k(i)),2);
  end
  
  k = find(isnan(xparam1));
  xparam1(k) = bayestopt_.pmean(k);
  
  