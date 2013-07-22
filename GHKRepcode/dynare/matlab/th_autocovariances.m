% Copyright (C) 2001 Michel Juillard
%
% computes the theoretical auto-covariances, Gamma_y, for an AR(p) process 
% with coefficients dr.ghx and dr.ghu and shock variances Sigma_e_
% for a subset of variables ivar (indices in lgy_)
% Theoretical HP filtering is available as an option

function [Gamma_y,ivar]=th_autocovariances(dr,ivar)
  global lgy_ endo_nbr exo_nbr Sigma_e_ ykmin_ lgx_orig_ord_ options_

  if sscanf(version('-release'),'%d') < 13
    warning off
  else
    eval('warning off MATLAB:dividebyzero')
  end
  nar = options_.ar;
  Gamma_y = cell(nar+1,1);
  if isempty(ivar)
    ivar = [1:endo_nbr]';
  end
  nvar = size(ivar,1);
  
  ghx = dr.ghx;
  ghu = dr.ghu;
  npred = dr.npred;
  nstatic = dr.nstatic;
  kstate = dr.kstate;
  order = dr.order_var;
  iv(order) = [1:length(order)];
  nx = size(ghx,2);
  
  ikx = [nstatic+1:nstatic+npred];
  
  A = zeros(nx,nx);
  k0 = kstate(find(kstate(:,2) <= ykmin_+1),:);
  i0 = find(k0(:,2) == ykmin_+1);
  i00 = i0;
  n0 = length(i0);
  A(i0,:) = ghx(ikx,:);
  AS = ghx(:,i0);
  ghu1 = zeros(nx,exo_nbr);
  ghu1(i0,:) = ghu(ikx,:);
  for i=ykmin_:-1:2
    i1 = find(k0(:,2) == i);
    n1 = size(i1,1);
    j1 = zeros(n1,1);
    j2 = j1;
    for k1 = 1:n1
      j1(k1) = find(k0(i00,1)==k0(i1(k1),1));
      j2(k1) = find(k0(i0,1)==k0(i1(k1),1));
    end
%    A(i1,i0(j2))=eye(n1);
    AS(:,j1) = AS(:,j1)+ghx(:,i1);
    i0 = i1;
  end
  b = ghu1*Sigma_e_*ghu1';


  [A,B] = kalman_transition_matrix(dr);
  % index of predetermined variables in A
  i_pred = [nstatic+(1:npred) endo_nbr+1:length(A)];
  A = A(i_pred,i_pred);

  [v,d] = eig(A);
  abs_root = abs(diag(d));
  i_ur = find(abs_root > 2-options_.qz_criterium & ...
	      abs_root < options_.qz_criterium);

  if ~isempty(i_ur)
    % non-stationary variables have non-zero entries in eigenvectors
    % associated with unit roots
    ns_var = find(any(abs(v(1:npred,i_ur))>1e-7,2));

    %right eigenvectors
    v1 = inv(v);

    % remove zero frequency from A
    A = A-real(v(:,i_ur)*d(i_ur,i_ur)*v1(i_ur,:));
    
    % return only variance of stationary variables
    i_ivar = find(~ismember(ivar,dr.order_var(ns_var+nstatic)));
    ivar = ivar(i_ivar);
  end
  % order of variables with preset variances in ghx and ghu
  iky = iv(ivar);
  
  aa = ghx(iky,:);
  bb = ghu(iky,:);

  if options_.order == 2
    vx =  lyapunov_symm(A,b); 
    Ex = (dr.ghs2(ikx)+dr.ghxx(ikx,:)*vx(:)+dr.ghuu(ikx,:)*Sigma_e_(:))/2;
    Ex = (eye(n0)-AS(ikx,:))\Ex;
    Gamma_y{nar+3} = AS(iky,:)*Ex+(dr.ghs2(iky)+dr.ghxx(iky,:)*vx(:)+dr.ghuu(iky,:)*Sigma_e_(:))/2;
  end
  if options_.hp_filter == 0
    if options_.order < 2
      vx =  lyapunov_symm(A,b);
    end
    Gamma_y{1} = aa*vx*aa'+ bb*Sigma_e_*bb';
    k = find(abs(Gamma_y{1}) < 1e-12);
    Gamma_y{1}(k) = 0;
    
    % autocorrelations
    if nar > 0
      vxy = (A*vx*aa'+ghu1*Sigma_e_*bb');

      sy = sqrt(diag(Gamma_y{1}));
      sy = sy *sy';
      Gamma_y{2} = aa*vxy./sy;
      
      for i=2:nar
	vxy = A*vxy;
	Gamma_y{i+1} = aa*vxy./sy;
      end
    end
    
    % variance decomposition
    if exo_nbr > 1
      Gamma_y{nar+2} = zeros(length(ivar),exo_nbr);
      SS(lgx_orig_ord_,lgx_orig_ord_)=Sigma_e_+1e-14*eye(exo_nbr);
      cs = chol(SS)';
      b1(:,lgx_orig_ord_) = ghu1;
      b1 = b1*cs;
      b2(:,lgx_orig_ord_) = ghu(iky,:);
      b2 = b2*cs;
      vx  = lyapunov_symm(A,b1*b1');
      vv = diag(aa*vx*aa'+b2*b2');
      for i=1:exo_nbr
	vx1 = lyapunov_symm(A,b1(:,i)*b1(:,i)');
	Gamma_y{nar+2}(:,i) = abs(diag(aa*vx1*aa'+b2(:,i)*b2(:,i)'))./vv;
      end
    end
  else
    lambda = options_.hp_filter;
    ngrid = options_.hp_ngrid;
    freqs = 0 : ((2*pi)/ngrid) : (2*pi*(1 - .5/ngrid)); 
    tpos  = exp( sqrt(-1)*freqs);
    tneg  =  exp(-sqrt(-1)*freqs);
    hp1 = 4*lambda*(1 - cos(freqs)).^2 ./ (1 + 4*lambda*(1 - cos(freqs)).^2);
    
    mathp_col = [];
    IA = eye(size(A,1));
    IE = eye(exo_nbr);
    for ig = 1:ngrid
      f_omega  =(1/(2*pi))*( [inv(IA-A*tneg(ig))*ghu1;IE]...
			     *Sigma_e_*[ghu1'*inv(IA-A'*tpos(ig)) ...
		    IE]); % state variables
      g_omega = [aa*tneg(ig) bb]*f_omega*[aa'*tpos(ig); bb']; % selected variables
      f_hp = hp1(ig)^2*g_omega; % spectral density of selected filtered series
      mathp_col = [mathp_col ; (f_hp(:))'];    % store as matrix row
                                               % for ifft
    end;

    % covariance of filtered series
    imathp_col = real(ifft(mathp_col))*(2*pi);

    Gamma_y{1} = reshape(imathp_col(1,:),nvar,nvar);
    
    % autocorrelations
    if nar > 0
      sy = sqrt(diag(Gamma_y{1}));
      sy = sy *sy';
      for i=1:nar
	Gamma_y{i+1} = reshape(imathp_col(i+1,:),nvar,nvar)./sy;
      end
    end
    
    %variance decomposition
    if exo_nbr > 1 
      Gamma_y{nar+2} = zeros(nvar,exo_nbr);
      SS(lgx_orig_ord_,lgx_orig_ord_)=Sigma_e_+1e-14*eye(exo_nbr);
      cs = chol(SS)';
      SS = cs*cs';
      b1(:,lgx_orig_ord_) = ghu1;
      b2(:,lgx_orig_ord_) = ghu(iky,:);
      mathp_col = [];
      IA = eye(size(A,1));
      IE = eye(exo_nbr);
      for ig = 1:ngrid
	f_omega  =(1/(2*pi))*( [inv(IA-A*tneg(ig))*b1;IE]...
			       *SS*[b1'*inv(IA-A'*tpos(ig)) ...
		    IE]); % state variables
	g_omega = [aa*tneg(ig) b2]*f_omega*[aa'*tpos(ig); b2']; % selected variables
	f_hp = hp1(ig)^2*g_omega; % spectral density of selected filtered series
	mathp_col = [mathp_col ; (f_hp(:))'];    % store as matrix row
						 % for ifft
      end;

      imathp_col = real(ifft(mathp_col))*(2*pi);
      vv = diag(reshape(imathp_col(1,:),nvar,nvar));
      for i=1:exo_nbr
	mathp_col = [];
	SSi = cs(:,i)*cs(:,i)';
	for ig = 1:ngrid
	  f_omega  =(1/(2*pi))*( [inv(IA-A*tneg(ig))*b1;IE]...
				 *SSi*[b1'*inv(IA-A'*tpos(ig)) ...
		    IE]); % state variables
	  g_omega = [aa*tneg(ig) b2]*f_omega*[aa'*tpos(ig); b2']; % selected variables
	  f_hp = hp1(ig)^2*g_omega; % spectral density of selected filtered series
	  mathp_col = [mathp_col ; (f_hp(:))'];    % store as matrix row
						   % for ifft
	end;

	imathp_col = real(ifft(mathp_col))*(2*pi);
	Gamma_y{nar+2}(:,i) = abs(diag(reshape(imathp_col(1,:),nvar,nvar)))./vv;
      end
    end
  end
  if sscanf(version('-release'),'%d') < 13
    warning on
  else
    eval('warning on MATLAB:dividebyzero')
  end
  
  % 10/18/2002 MJ
  % 10/30/2002 added autocorrelations, Gamma_y is now a cell array
  % 01/20/2003 MJ added variance decomposition
  % 02/18/2003 MJ added HP filtering (Thanks to Jean Chateau for the code)
  % 04/28/2003 MJ changed handling of options
  % 05/19/2003 MJ don't assume lags are in increasing order in building A
  % 05/21/2003 MJ added global i_exo_var_,
  %               variance decomposition: test exo_nbr > 1
  % 05/29/2003 MJ removed global i_exo_var_
  % 06/10/2003 MJ test release for warning syntax