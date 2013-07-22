function [Sigma_e_,info] = calib(var_indices,targets,var_weights,nar,cova,Sigma_e_)
  global ys_ endo_nbr exo_nbr lgy_ lgx_ ykmin_ vx options_
  
  ncstr = 0;
  ni = size(var_indices,1);
  for i=1:nar+3
    ncstr = ncstr + size(var_indices{i},1);
  end
  if cova
    if ncstr < exo_nbr*(exo_nbr+1)/2
      error(['number of preset variances is smaller than number of shock' ...
	     ' variances and covariances to be estimated !'])
    end
  else
    if ncstr < exo_nbr
      error(['number of preset variances is smaller than number of shock' ...
	     ' variances to be estimated !'])
    end
  end
  
  check_model;
  
  % computes approximate solution at order 1
  [dr_, info] = resol(ys_,0);

  if info(1)
    print_info(info);
  end  

  ghx = dr.ghx;
  ghu = dr.ghu;
  npred = dr.npred;
  nstatic = dr.nstatic;
  kstate = dr.kstate;
  order = dr.order_var;
  iv(order) = [1:endo_nbr];
  iv = iv';
  nx = size(ghx,2);

  ikx = [nstatic+1:nstatic+npred];
  
  A = zeros(nx,nx);
  A(1:npred,:)=ghx(ikx,:);
  offset_r = npred;
  offset_c = 0;
  i0 = find(kstate(:,2) == ykmin_+1);
  n0 = size(i0,1);
  for i=ykmin_:-1:2
    i1 = find(kstate(:,2) == i);
    n1 = size(i1,1);
    j = zeros(n1,1);
    for j1 = 1:n1
      j(j1) = find(kstate(i0,1)==kstate(i1(j1),1));
    end
    A(offset_r+1:offset_r+n1,offset_c+j)=eye(n1);
    offset_r = offset_r + n1;
    offset_c = offset_c + n0;
    i0 = i1;
    n0 = n1;
  end
  ghu1 = [ghu(ikx,:);zeros(nx-npred,exo_nbr)];
%  IA = speye(nx*nx)-kron(A,A);
%  kron_ghu = kron(ghu1,ghu1);
  
  % vx1 such that vec(sigma_x) = vx1 * vec(Sigma_e_) (predetermined vars) 
vx1 = [];
  % vx1 = IA\kron_ghu;
  IA = [];
  kron_ghu = [];
  
  % computes required variables and indices among required variables
  iiy = [];
  for i=1:nar+3
    if i ~= 3 & ~isempty(var_indices{i})
      iiy = union(iiy, iv(var_indices{i}(:,1)));
    end
  end
  if ~isempty(var_indices{2})
    iiy = union(iiy, iv(var_indices{2}(:,2)));
  end
  ny = size(iiy,1);

  for i=1:nar+3
    if i ~= 3 & ~isempty(var_indices{i})
      var_indices{i}(:,1) = indnv(iv(var_indices{i}(:,1)),iiy);
    end
    if i ~= 2 & i ~= 3 & ~isempty(var_indices{i})
      var_indices{i} = sub2ind([ny ny],var_indices{i},var_indices{i});
    end
  end
  if ~isempty(var_indices{2})
    var_indices{2}(:,2) = indnv(iv(var_indices{2}(:,2)),iiy);
    var_indices{2} = sub2ind([ny ny],var_indices{2}(:,1),var_indices{2}(:,2));
  end
  if ~isempty(var_indices{3})
    var_indices{3} = sub2ind([exo_nbr exo_nbr],var_indices{3}(:,1),var_indices{3}(:,2));
  end
  if isempty(Sigma_e_)
    Sigma_e_ = 0.01*eye(exo_nbr);
    b = 0.1*ghu1*ghu1';
  else
    b = ghu1*Sigma_e_*ghu1';
    Sigma_e_ = chol(Sigma_e_+1e-14*eye(exo_nbr));
  end
  options=optimset('LargeScale','on','MaxFunEvals',20000*ny,'TolX',1e-4, ...
		   'TolFun',1e-4,'Display','Iter','MaxIter',10000);
%  [Sigma_e_,f]=fminunc(@calib_obj,Sigma_e_,options,A,ghu1,ghx(iiy,:),ghu(iiy,:),targets,var_weights,var_indices,nar);
  [Sigma_e_,f]=fmincon(@calib_obj,diag(Sigma_e_).^2,-eye(exo_nbr),zeros(exo_nbr,1),[],[],[],[],[],options,A,ghu1,ghx(iiy,:),ghu(iiy,:),targets,var_weights,var_indices,nar);
  Sigma_e_ = diag(Sigma_e_);
  
  objective = calib_obj2(diag(Sigma_e_),A,ghu1,ghx(iiy,:),ghu(iiy,:), ...
			 targets,var_weights,var_indices,nar);
  if ~options_.noprint
    disp('CALIBRATION')
    disp('')
    if ~isempty(var_indices{1})
      disp(sprintf('%16s %14s %14s %14s %14s','Std. Dev','Target','Obtained','Diff'));
      str = '   ';
      for i=1:size(var_indices{1},1)
	[i1,i2] = ind2sub([ny ny],var_indices{1}(i));
	str = sprintf('%16s: %14.2f %14.2f %14.2f',lgy_(order(iiy(i1)),:),targets{1}(i),objective{1}(i),objective{1}(i)-targets{1}(i));
	disp(str);
      end
    end
    if ~isempty(var_indices{2})
      disp(sprintf('%32s %14s %14s','Correlations','Target','Obtained','Diff'));
      str = '   ';
      for i=1:size(var_indices{2},1)
	[i1,i2]=ind2sub([ny ny],var_indices{2}(i));
	str = sprintf('%16s,%16s: %14.2f %14.2f %14.2f',lgy_(order(iiy(i1)),:), ...
		      lgy_(order(iiy(i2)),:),targets{2}(i),objective{2}(i),objective{2}(i)-targets{2}(i));
	disp(str);
      end
    end
    if ~isempty(var_indices{3})
      disp(sprintf('%32s %16s %16s','Constrained shocks (co)variances','Target','Obtained'));
      str = '   ';
      for i=1:size(var_indices{3},1)
	[i1,i2]=ind2sub([exo_nbr exo_nbr],var_indices{3}(i));
	if i1 == i2
	  str = sprintf('%32s: %16.4f %16.4f',lgx_(order(i1),:), ...
			targets{3}(i),objective{3}(i));
	else
	  str = sprintf('%16s,%16s: %16.4f %16.4f',lgx_(order(i1),:), ...
			lgx_(order(i2), :),targets{3}(i),objective{3}(i));
	end
	disp(str);
      end
    end
    flag = 1;
    for j=4:nar+3
      if ~isempty(var_indices{j})
	if flag
	  disp(sprintf('%16s %16s %16s','Autocorrelations','Target','Obtained'));
	  str = '   ';
	  flag = 0;
	end
	for i=1:size(var_indices{j},1)
	  [i1,i2] = ind2sub([ny ny],var_indices{j}(i));
	  str = sprintf('%16s(%d): %16.4f %16.4f',lgy_(order(iiy(i1)),:), ...
			j-3,targets{j}(i),objective{j}(i));
	  disp(str);
	end
      end
    end    
    
    disp('');
    disp('Calibrated variances')
    str = '   ';
    for i=1:exo_nbr
      str = [str sprintf('%16s',lgx_(i,:))];
    end
    disp(str);
    disp('');
    str = '      ';
    for i=1:exo_nbr
      str = [str sprintf('%16f',Sigma_e_(i,i))];
    end
    disp(str);
  end % end if options_.noprint

  
  % 10/9/02 MJ