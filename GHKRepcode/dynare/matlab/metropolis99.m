function [PostMod,PostVar,Scale,PostMean] = ...
    samaxwell(xparam1,gend,data,mh_bounds,NumberOfIterationsForInitialization,init_scale,info,MeanPar,VarCov)
% stephane.adjemian@cepremap.cnrs.fr [11-22-2005]
global bayestopt_ exo_nbr dr_ estim_params_ Sigma_e_ options_ xparam_test
global lgy_ lgx_ fname_ ys_ xkmin_ xkmax_ ykmin_ ykmax_ endo_nbr mean_varobs
global oo_ lgx_orig_ord_ lgy_TeX_ lgx_TeX_ dsge_prior_weight
nvx   	= estim_params_.nvx;
nvn   	= estim_params_.nvn;
ncx   	= estim_params_.ncx;
ncn   	= estim_params_.ncn;
np    	= estim_params_.np ;
ns      = nvx+nvn+ncx+ncn; 
nx    	= nvx+nvn+ncx+ncn+np;
npar  	= length(xparam1);
nvobs 	= size(options_.varobs,1);
options_.lik_algo = 1;
MaxNumberOfTuningSimulations = 100000;
MaxNumberOfClimbingSimulations = 20000;
AcceptanceTarget = 0.33;

%% [1] I build a covariance matrix for the jumping distribution.
stdev = bayestopt_.pstdev;
indx = find(isinf(stdev));
stdev(indx) = 10*ones(length(indx),1);
CovJump = diag(stdev).^2;
if nargin == 9
  CovJump = VarCov;
end
ModePar = xparam1;
%% [2] For this covariance matrix I tune the scale parameter.
hh = waitbar(0,'Tuning of the scale parameter...');
set(hh,'Name','Tuning of the scale parameter (1)')
j = 1; jj  = 1;
isux = 0;
jsux = 0;
test = 0;
% TEST = 0;
ix2 = ModePar;% initial condition! 
if isempty(strmatch('dsge_prior_weight',estim_params_.param_names)) & isempty(dsge_prior_weight)
  ilogpo2 = -DsgeLikelihood(ix2,gend,data);% initial posterior density
else
  ilogpo2 = -DsgeVarLikelihood(ix2,gend);% initial posterior density
end
mlogpo2 = ilogpo2;
dd = transpose(chol(CovJump));
while j<=MaxNumberOfTuningSimulations
  proposal = init_scale*dd*randn(npar,1) + ix2;
  if all(proposal > mh_bounds(:,1)) & all(proposal < mh_bounds(:,2))
    if isempty(strmatch('dsge_prior_weight',estim_params_.param_names)) & isempty(dsge_prior_weight)
      logpo2 = -DsgeLikelihood(proposal,gend,data);
    else
      logpo2 = -DsgeVarLikelihood(proposal,gend);
    end
  else
    logpo2 = -inf;
  end
  % I move if the proposal is enough likely...
  if logpo2 > -inf & log(rand) < logpo2 - ilogpo2
    ix2 = proposal; 
    if logpo2 > mlogpo2
      ModePar = proposal;
      mlogpo2 = logpo2;
    end
    ilogpo2 = logpo2;
    isux = isux + 1;
    jsux = jsux + 1;
  else
    % ... otherwise I don't move.
  end	
  prtfrc = j/MaxNumberOfTuningSimulations;
  waitbar(prtfrc,hh,sprintf('%f done, acceptation rate %f',prtfrc,isux/j));
  if  j/500 == round(j/500)
    test1 = jsux/jj;
    test2 = isux/j;
    cfactor = (test1/AcceptanceTarget)^1.0*(test2/AcceptanceTarget)^0.0;
    if cfactor>0
      init_scale = init_scale*cfactor;
    end
    jsux = 0; jj = 0;
    if test1>AcceptanceTarget*0.90 & test1<AcceptanceTarget*1.10
      test = test+1;
    end
    if test>6
      break
    end
  end
  j = j+1;
  jj = jj + 1;
end
close(hh);
%% [3] One block metropolis, I update the covariance matrix of the jumping distribution
hh = waitbar(0,'Metropolis-Hastings...');
set(hh,'Name','Looking for the posterior covariance...')
j = 1;
isux = 0;
% ix2 = xparam1;% initial condition! 
if isempty(strmatch('dsge_prior_weight',estim_params_.param_names)) & Isempty(dsge_prior_weight)
  ilogpo2 = -DsgeLikelihood(ix2,gend,data);
else
  ilogpo2 = -DsgeVarLikelihood(ix2,gend);
end
dd = transpose(chol(CovJump));
while j<=NumberOfIterationsForInitialization
  proposal = init_scale*dd*randn(npar,1) + ix2;
  if all(proposal > mh_bounds(:,1)) & all(proposal < mh_bounds(:,2))
    if isempty(strmatch('dsge_prior_weight',estim_params_.param_names)) & isempty(dsge_prior_weight)
      logpo2 = -DsgeLikelihood(proposal,gend,data);
    else
      logpo2 = -DsgeVarLikelihood(proposal,gend);
    end
  else
    logpo2 = -inf;
  end
  % I move if the proposal is enough likely...
  if logpo2 > -inf & log(rand) < logpo2 - ilogpo2
    ix2 = proposal; 
    if logpo2 > mlogpo2
      ModePar = proposal;
      mlogpo2 = logpo2;
    end
    ilogpo2 = logpo2;
    isux = isux + 1;
    jsux = jsux + 1;
  else
    % ... otherwise I don't move.
  end	
  prtfrc = j/NumberOfIterationsForInitialization;
  waitbar(prtfrc,hh,sprintf('%f done, acceptation rate %f',prtfrc,isux/j));
  % I update the covariance matrix and the mean:
  oldMeanPar = MeanPar;
  MeanPar = oldMeanPar + (1/j)*(ix2-oldMeanPar);
  CovJump = CovJump + oldMeanPar*oldMeanPar' - MeanPar*MeanPar' + ...
                 (1/j)*(ix2*ix2' - CovJump - oldMeanPar*oldMeanPar');
  j = j+1;
end
close(hh);
%% [4 & 5] I tune the scale parameter (with the new covariance matrix) if
%% this is the last call to the routine, and I climb the hill (without
%% updating the covariance matrix)...
if strcmpi(info,'LastCall')
  % MeanPar = xparam1;  
  % ModePar = xparam1; I'm stupid!
  hh = waitbar(0,'Tuning of the scale parameter...');
  set(hh,'Name','Tuning of the scale parameter (2)')
  j = 1; jj  = 1;
  isux = 0;
  jsux = 0;
  test = 0;
  % ix2 = ModePar;% initial condition! Again stupid 
  if isempty(strmatch('dsge_prior_weight',estim_params_.param_names)) & isempty(dsge_prior_weight)
    ilogpo2 = -DsgeLikelihood(ix2,gend,data);% initial posterior density
  else
    ilogpo2 = -DsgeVarLikelihood(ix2,gend);% initial posterior density
  end
  dd = transpose(chol(CovJump));
  while j<=MaxNumberOfTuningSimulations
    proposal = init_scale*dd*randn(npar,1) + ix2;
    if all(proposal > mh_bounds(:,1)) & all(proposal < mh_bounds(:,2))
      if isempty(strmatch('dsge_prior_weight',estim_params_.param_names)) & isempty(dsge_prior_weight)
	logpo2 = -DsgeLikelihood(proposal,gend,data);
      else
	logpo2 = -DsgeVarLikelihood(proposal,gend);
      end
    else
      logpo2 = -inf;
    end
    % I move if the proposal is enough likely...
    if logpo2 > -inf & log(rand) < logpo2 - ilogpo2
      ix2 = proposal;
      if logpo2 > mlogpo2
	ModePar = proposal;
	mlogpo2 = logpo2;
      end
      ilogpo2 = logpo2;
      isux = isux + 1;
      jsux = jsux + 1;
    else
      % ... otherwise I don't move.
    end
    prtfrc = j/MaxNumberOfTuningSimulations;
    waitbar(prtfrc,hh,sprintf('%f done, acceptation rate %f',prtfrc,isux/j));  
    if  j/500 == round(j/500) 
      test1 = jsux/jj;
      test2 = isux/j;  
      cfactor = (test1/AcceptanceTarget)^1.0*(test2/AcceptanceTarget)^0.0;
      init_scale = init_scale*cfactor;
      jsux = 0; jj = 0;
      if test1>AcceptanceTarget*0.90 & test1<AcceptanceTarget*1.10
	test = test+1;
      end
      if test>6
	break
      end
    end
    j = j+1;
    jj = jj + 1;
  end
  close(hh);
  PostVar = CovJump;
  PostMean = MeanPar;
  Scale = init_scale;
  %%
  %% Now I climb the hill
  %%
  hh = waitbar(0,' ');
  set(hh,'Name','Now I am climbing the hill...')
  j = 1; jj  = 1;
  jsux = 0;
  test = 0;
  dd = transpose(chol(CovJump));
  init_scale = eye(length(dd))*init_scale;
  indx = strmatch('dsge_prior_weight',estim_params_.param_names,'exact');
  while j<=MaxNumberOfClimbingSimulations
    proposal = init_scale*dd*randn(npar,1) + ModePar;
    if all(proposal > mh_bounds(:,1)) & all(proposal < mh_bounds(:,2))
      if isempty(strmatch('dsge_prior_weight',estim_params_.param_names)) & isempty(dsge_prior_weight)
	logpo2 = -DsgeLikelihood(proposal,gend,data);
      else
	logpo2 = -DsgeVarLikelihood(proposal,gend);
      end
    else
      logpo2 = -inf;
    end
    if logpo2 > mlogpo2% I move if the proposal is higher...
      ModePar = proposal;
      mlogpo2 = logpo2;
      jsux = jsux + 1;
    else
      % otherwise I don't move...
    end
    prtfrc = j/MaxNumberOfClimbingSimulations;
    waitbar(prtfrc,hh,sprintf('%f done, acceptation rate %f',prtfrc,jsux/jj));  
    if  j/200 == round(j/200) 
      test1 = jsux/jj;
      %cfactor = (test1/AcceptanceTarget)^1.0*(test2/AcceptanceTarget)^0.0;
      %init_scale = init_scale*cfactor;
      jsux = 0; 
      jj = 0;
      if test1<0.001
	test = test+1;
      end
      if test>4% If I do not progress enough I reduce the scale parameter
               % of the jumping distribution (cooling down the system).
	if isempty(indx)
          init_scale = init_scale/1.10;
        else
          init_scale = init_scale/1.10;
          % USEFUL IF THE MODEL IS TRUE:
          % init_scale(indx+ns,indx+ns) = init_scale(indx+ns,indx+ns)*1.10/1.001; 
        end
      end
    end
    j = j+1;
    jj = jj + 1;
  end
  close(hh);
else
  PostVar = CovJump;
  PostMean = MeanPar;
  Scale = init_scale;
end
PostMod = ModePar;