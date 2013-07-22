function metropolis(xparam1,vv,gend,data,rawdata,mh_bounds)
% stephane.adjemian@cepremap.cnrs.fr [07-31-2004]
% Adapted from an older version of metropolis.m 

  global bayestopt_ exo_nbr dr_ estim_params_ Sigma_e_ options_ xparam_test
  global lgy_ lgx_ fname_ ys_ xkmin_ xkmax_ ykmin_ ykmax_ endo_nbr
  global oo_ lgx_orig_ord_ lgy_TeX_ lgx_TeX_
  global dsge_prior_weight

  TeX       = options_.TeX;
  nruns     = options_.mh_replic;
  truns     = options_.mh_replic*options_.mh_nblck;
  nblck     = options_.mh_nblck;
  nvx       = estim_params_.nvx;
  nvn       = estim_params_.nvn;
  ncx       = estim_params_.ncx;
  ncn       = estim_params_.ncn;
  np        = estim_params_.np ;
  nx        = nvx+nvn+ncx+ncn+np;
  npar      = length(xparam1);
  nvobs     = size(options_.varobs,1);
  horizon = options_.forecast;
  bayestopt_.penalty = 1e8;
  
  % options_.load_mh_file = -1;
  

  %% Determine the value of MAX_nruns, MAX_nforc, MAX_nsmoo and MAX_ninno values
  MaxNumberOfBytes = 1000000;%% This value should be adjusted
  MAX_nruns = ceil(MaxNumberOfBytes/(npar+2)/8);
  MAX_nforc = ceil(MaxNumberOfBytes/((options_.forecast+ykmin_)*length(ys_))/8);
  MAX_nsmoo = ceil(MaxNumberOfBytes/((endo_nbr)*gend)/8);
  MAX_ninno = ceil(MaxNumberOfBytes/(exo_nbr*gend)/8);
  MAX_nerro = ceil(MaxNumberOfBytes/(size(options_.varobs,1)*gend)/8);
  MAX_nfilt = ceil(MaxNumberOfBytes/((endo_nbr)*gend)/8);
  if options_.bayesian_irf
    MAX_nirfs_dsge = ceil(MaxNumberOfBytes/(options_.irf*length(ys_)*exo_nbr)/8);
    if ~isempty(dsge_prior_weight)
      MAX_nirfs_dsgevar = ceil(MaxNumberOfBytes/(options_.irf*nvobs*exo_nbr)/8)
    end
  end
  MAX_nthm1 = ceil(MaxNumberOfBytes/(length(ys_)*8));
  MAX_nthm2 = ceil(MaxNumberOfBytes/(length(ys_)*length(ys_)*8));
  MAX_nthm3 = ceil(MaxNumberOfBytes/(length(ys_)*exo_nbr*8));
  MAX_nthm4 = ceil(MaxNumberOfBytes/(length(ys_)*options_.ar*8));

  d     = chol(vv);
  options_.lik_algo = 1;

  if nruns
    if options_.load_mh_file == 0
      % Delete old mh files...
      if nblck > 1
	disp('MH: Multiple chains mode.')
      else
	disp('MH: One Chain mode.')
      end
      files = eval(['dir(''' fname_ '_mh*.mat'');']);
      if size(files,1)
	delete([fname_ '_mh*.mat']);
	disp('MH: Old _mh files succesfully erased!')
      end   
      nops = 0;         % Number Of Past Simulations.
      lfile = -1;       % Index for the last mh file.
      if nblck > 1
	disp('MH: Searching for initial values...')
	ix2 = zeros(1,npar,nblck);
	ilogpo2 = zeros(1,nblck);
	for j=1:nblck
	  validate  = 0;
	  init_iter = 0;
	  trial     = 1;
	  while validate == 0 & trial <= 10 
	    candidate = options_.mh_init_scale*randn(1,npar)*d + transpose(xparam1);
	    if all(candidate' > mh_bounds(:,1)) & all(candidate' < mh_bounds(:,2)) 
	      ix2(1,:,j) = candidate;
	      if isempty(strmatch('dsge_prior_weight',estim_params_.param_names)) & isempty(dsge_prior_weight)
                ilogpo2(1,j) = -DsgeLikelihood(ix2(1,:,j)',gend,data);
	      else
                ilogpo2(1,j) = -DsgeVarLikelihood(ix2(1,:,j)',gend);
	      end
	      j = j+1;
	      validate = 1;
	    end
	    init_iter = init_iter + 1;
	    if init_iter > 100 & validate == 0
	      disp(['MH: I couldn''t get a valid initial value in 100 trials.'])
	      disp(['MH: You should Reduce mh_init_scale...'])
	      disp(sprintf('MH: Parameter mh_init_scale is equal to %f.',options_.mh_init_scale))
	      options_.mh_init_scale = input('MH: Enter a new value...  ');
	      trial = trial+1;
	    end
	  end
	  if trial > 10 & ~validate
	    error(['MH: I''m unable to find a starting value for block ' int2str(j)]);
	  end
	end
	disp('MH: Initial values found!')
	disp(' ')
      else
	candidate = transpose(xparam1);
	if all(candidate' > mh_bounds(:,1)) & all(candidate' < mh_bounds(:,2)) 
	  ix2 = candidate;
	  if isempty(strmatch('dsge_prior_weight',estim_params_.param_names)) & isempty(dsge_prior_weight)
	    ilogpo2 = -DsgeLikelihood(ix2',gend,data);
	  else
	    ilogpo2 = -DsgeVarLikelihood(ix2',gend);
	  end
	  disp('MH: Initialization at the posterior mode.')
	  disp(' ')
	else
	  disp('MH: Initialization failed...')
	  error('MH: The posterior mode lies outside the prior bounds.')
	end
      end
      save([fname_ '_MhInitialization'],'ix2','ilogpo2');
    elseif options_.load_mh_file == 1
      disp('MH: I''m loading past metropolis-hastings simulations...')
      if nblck>1
	files = eval(['dir(''' fname_ '_mh*_blck*.mat'');']);
	bfiles = eval(['dir(''' fname_ '_mh0_blck*.mat'');']);
	lfile = size(files,1)/size(bfiles,1)-1;
      else
	files = eval(['dir(''' fname_ '_mh*.mat'');']);
	bfiles = 1;
	ifiles = eval(['dir(''' fname_ '_MhInitialization.mat'');']);
	if isempty(ifiles)
	  lfile = size(files,1)-1;
	else
	  lfile = size(files,1)-2;
	end
      end
      if ~size(files,1)
	error('MH: FAILURE :: there is no MH file to load here!')    
      end
      past_number_of_blocks = size(bfiles,1);
      if size(bfiles,1)>0 & past_number_of_blocks ~= nblck
	disp('MH: The specified number of blocks doesn''t match with the previous number of blocks!')
	disp(['MH: You declared ' int2str(nblck) ' blocks, but the previous number of blocks was ' int2str(past_number_of_blocks) '.'])
	disp(['MH: I will run the Metropolis-Hastings with ' int2str(past_number_of_blocks) ' blocks.' ])
	nblck = past_number_of_blocks;
	options_.mh_nblck = nblck;
      end
      %lfile = size(files,1)/nblck-1;
      if nblck == 1
	instr = [fname_ '_mh' int2str(lfile)];
	eval(['load ' instr]);
	clear post2;
	nops = size(logpo2,1);  
	ix2 = x2(nops,:);   
	ilogpo2 = logpo2(nops);
	clear x2  logpo2;     
	for file = 0:lfile-1
	  instr = [fname_ '_mh' int2str(file)];
	  eval(['load ' instr]);
	  clear post2 x2;
	  nops = nops + size(logpo2,1);
	end
      else 
	for b = 1:nblck
	  instr = [fname_ '_mh' int2str(lfile) '_blck' int2str(b)];
	  eval(['load ' instr]);
	  clear post2;
	  nops = length(logpo2);
	  ix2(1,:,b) = x2(nops,:);  
	  ilogpo2(b) = logpo2(nops);
	  clear x2  logpo2;     
	end
	for file = 0:lfile-1
	  instr = [fname_ '_mh' int2str(file) '_blck1'];
	  eval(['load ' instr]);
	  clear post2 x2;
	  nops = nops + length(logpo2);
	  clear logpo2;
	end
      end
      % nops is the Number Of Past Simulations. 
      disp(['MH: ... It''s done. I''ve loaded ' int2str(nops) 'simulations.'])
      disp(' ')
    elseif options_.load_mh_file == -1%%% Not ready...
      instr = [fname_ '_MhInitialization'];
      eval(['load ' instr]);
      nblck = length(ilogpo2);
      options_.mh_nblck = nblck;
      % Count the total number of saved mh files
      AllMhFiles = eval(['dir(''' fname_ '_mh*_blck*.mat'');']);
      TotalNumberOfMhFiles = size(AllMhFiles,1);
      % Count the number of saved mh files per block
      NumberOfMhFilesPerBlock = zeros(nblck,1); 
      for i = 1:nblck
	BlckMhFiles = eval(['dir(''' fname_ '_mh*_blck' int2str(i) '.mat'');']);
	NumberOfMhFilesPerBlock(i) = size(BlckMhFiles,1);
      end
      NumberOfMhFilesPerBlock
      return
    end    
    isux = 0; 
    if nblck == 1
      hh   = waitbar(0,'Please wait... Metropolis-Hastings...');
      set(hh,'Name','Metropolis-Hastings')
      if nruns <= MAX_nruns
	x2 = zeros(nruns,npar); 
	x2(1,:) = ix2(1,:);
	logpo2 = zeros(nruns,1);    
	logpo2(1) = ilogpo2;    
      else
	x2 = zeros(MAX_nruns,npar);
	x2(1,:) = ix2(1,:);
	logpo2 = zeros(MAX_nruns,1);
	logpo2(1) = ilogpo2;
      end
      irun = ~options_.load_mh_file;    %%%% irun=0 <-- previous files are loaded
      rruns = nruns-irun;
      j=1;
      while j<=rruns
	irun = irun + 1;
	if irun <= MAX_nruns
	  par = randn(1,npar)*d;
	  par = par.*bayestopt_.jscale' + ix2;  
	  if all(transpose(par) > mh_bounds(:,1)) & all(transpose(par) < mh_bounds(:,2))
	    if isempty(strmatch('dsge_prior_weight',estim_params_.param_names)) & isempty(dsge_prior_weight)
	      logpost = -DsgeLikelihood(transpose(par),gend,data);
	    else
	      logpost = -DsgeVarLikelihood(transpose(par),gend);
	    end
	  else
	    logpost = -inf;
	  end    
	  if logpost > -inf & log(rand) < logpost - ilogpo2
	    x2(irun,:) = par; 
	    ix2 = par;
	    logpo2(irun) = logpost; 
	    ilogpo2 = logpost;
	    isux = isux + 1;
	  else    
	    x2(irun,:) = ix2;
	    logpo2(irun) = ilogpo2;
	  end   
	  prtfrc = j/nruns;
	  waitbar(prtfrc,hh,sprintf('%f done, acceptation rate %f',prtfrc,isux/j));
	else
	  post2 = exp(logpo2);
	  save([fname_ '_mh' int2str(lfile+1)],'x2','logpo2','post2');
	  clear x2 logpo2 post2;
	  x2 = zeros(MAX_nruns,npar);
	  logpo2 = zeros(MAX_nruns,1);
	  lfile = lfile + 1;
	  irun = 0;
	  j = j - 1;
	end
	j = j + 1;
      end
      if nruns <= MAX_nruns
	post2 = exp(logpo2);
	save([fname_ '_mh' int2str(lfile+1)],'x2','logpo2','post2');
	clear post2 x2 logpo2;
      elseif irun <= MAX_nruns    
	x2 = x2(1:irun,:);
	logpo2 = logpo2(1:irun,1); 
	post2 = exp(logpo2);
	save([fname_ '_mh' int2str(lfile+1)],'x2','logpo2','post2');
	clear post2 x2 logpo2;
      end
      close(hh)
      disp(sprintf('Acceptation rate : %f',isux/nruns))
    else
      disp('Acceptation rates :')
      for b=1:nblck
	hh   = waitbar(0,'Please wait... Metropolis-Hastings...');
	set(hh,'Name',['Metropolis-Hastings, Block ',int2str(b)]);
	if nruns <= MAX_nruns
	  x2 = zeros(nruns,npar);   
	  x2(1,:) = ix2(1,:,b);
	  logpo2 = zeros(nruns,1);  
	  logpo2(1) = ilogpo2(1,b); 
	else
	  x2 = zeros(MAX_nruns,npar);
	  x2(1,:) = ix2(1,:,b);
	  logpo2 = zeros(MAX_nruns,1);
	  logpo2(1) = ilogpo2(1,b);
	end 
	irun  = ~options_.load_mh_file; % Previous files are loaded <-- irun=0
	rruns = nruns-irun;
	isav = 0;
	isux = 0;
	j = 1;
	while j <= rruns
	  irun = irun + 1;
	  if irun <= MAX_nruns
	    par = randn(1,npar)*d;
	    par = par.*transpose(bayestopt_.jscale) + ix2(1,:,b);  
	    if all(transpose(par) > mh_bounds(:,1)) & all(transpose(par) < mh_bounds(:,2))
	      if isempty(strmatch('dsge_prior_weight',estim_params_.param_names)) & isempty(dsge_prior_weight)
		logpost = -DsgeLikelihood(transpose(par),gend,data);
	      else
		logpost = -DsgeVarLikelihood(transpose(par),gend);
	      end
	    else
	      logpost = -inf;
	    end    
	    if logpost > -inf & log(rand) < logpost - ilogpo2(1,b)
	      x2(irun,:) = par; 
	      ix2(1,:,b) = par;
	      logpo2(irun) = logpost; 
	      ilogpo2(1,b) = logpost;
	      isux = isux + 1;
	    else    
	      x2(irun,:) = ix2(1,:,b);
	      logpo2(irun) = ilogpo2(1,b);
	    end 
	    prtfrc = j/nruns;
	    waitbar(prtfrc,hh,sprintf('%f done, acceptation rate %f',prtfrc,isux/j));
	  else
	    post2 = exp(logpo2);
	    save([fname_ '_mh' int2str(lfile+1+isav) '_blck' int2str(b)],'x2','logpo2','post2');
	    clear post2;
	    x2 = zeros(MAX_nruns,npar);
	    logpo2 = zeros(MAX_nruns,1);
	    isav = isav + 1;
	    irun = 0;
	    j=j-1;
	  end
	  j = j+1; 
	end
	if nruns <= MAX_nruns
	  post2 = exp(logpo2);
	  save([fname_ '_mh' int2str(lfile+isav+1) '_blck' int2str(b)],'x2','logpo2','post2');
	  clear post2 x2 logpo2;
	elseif irun <= MAX_nruns    
	  x2 = x2(1:irun,:);
	  logpo2 = logpo2(1:irun,1); 
	  post2 = exp(logpo2);
	  save([fname_ '_mh' int2str(lfile+isav+1) '_blck' int2str(b)],'x2','logpo2','post2');
	  clear post2 x2 logpo2;
	end
	disp(sprintf('Block %d: %f',b,isux/nruns))
	close(hh)
      end
    end
    disp(' ')
    disp(['MH: Total number of iterations       : ' int2str(nops+nruns) '.'])
  end %end if nruns
  if nblck == 1
    files = eval(['dir(''' fname_ '_mh*.mat'');']);
    nfile = size(files,1)-2;
    number_of_simulations_per_file = zeros(nfile+1,1);
    instr = [fname_ '_mh' int2str(0)];
    eval(['load ' instr]);
    clear x2 post2;
    number_of_simulations_per_file(1) = length(logpo2);
    if nfile >= 1
      for file = 1:nfile
	instr = [fname_ '_mh' int2str(file)];
	eval(['load ' instr]);
	clear post2 x2;
	number_of_simulations_per_file(file+1) = length(logpo2);
      end
    end
    clear logpo2;
    if ~nruns
      tmp  = cumsum(number_of_simulations_per_file);
      nops = tmp(nfile+1); clear tmp;
    end
  else
    files = eval(['dir(''' fname_ '_mh*_blck1.mat'');']);   
    nfile = size(files,1)-1;
    number_of_simulations_per_file = zeros(nfile+1,1);
    instr = [fname_ '_mh' int2str(0) '_blck' int2str(1)];
    eval(['load ' instr]);
    clear x2 post2;
    number_of_simulations_per_file(1) = length(logpo2);
    if nfile >= 1
      for file = 1:nfile
	instr = [fname_ '_mh' int2str(file) '_blck1'];
	eval(['load ' instr]);
	clear post2 x2;
	number_of_simulations_per_file(file+1) = length(logpo2);
      end
    end
    clear logpo2;
    if ~nruns
      tmp  = cumsum(number_of_simulations_per_file);
      nops = tmp(nfile+1); clear tmp;
      bfiles = eval(['dir(''' fname_ '_mh0_blck*.mat'');']);
      past_number_of_blocks = size(bfiles,1);
      if past_number_of_blocks ~= nblck
	nblck = past_number_of_blocks;
	options_.mh_nblck = nblck;
      end
    end
  end
  cumulated_number_of_simulations_per_file = cumsum(number_of_simulations_per_file);
  disp(['MH: Number of mh files             : ' int2str(nfile+1) ' per block.'])
  disp(['MH: Total number of generated files    : ' int2str((nfile+1)*nblck) '.'])
  disp(['MH: Total number of iterations         : ' int2str(nops+nruns) '.'])
  disp('MH: Number of simulations per file: ')
  for i=0:nfile
    disp(sprintf('    The number of simulations in file %d is: %d.',i,number_of_simulations_per_file(i+1)))
  end
  disp(' ')
  nsim = nops+nruns;
  %
  %%
  %%%
  %%%%
  %%%%% MCMC convergence diagnostics
  %%%%
  %%%
  %%
  %
  origin = 1000;
  if ~options_.nodiagnostic & nblck > 1 & nsim > origin
    %%
    %%  Univariate diagnostic : Brooks and Gelman (1998).
    %%
    step_size   = ceil((nsim-origin)/100);  % So that the computational time does not 
    ALPHA       = 0.2;                      % increase too much with the number of simulations. 
    time = 1:nsim;
    xx = origin:step_size:nsim;
    number_of_lines = length(xx);
    tmp = zeros(nsim*nblck,3);
    UDIAG = zeros(number_of_lines,6,npar);
    if nsim < origin
      error('MH: The number of simulations is to small to compute the MCMC convergence diagnostics.')
    end
    if TeX
      fidTeX = fopen([fname_ '_UnivariateDiagnostics.TeX'],'w');
      fprintf(fidTeX,'%% TeX eps-loader file generated by metropolis.m (Dynare).\n');
      fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
      fprintf(fidTeX,' \n');
    end
    disp('MH: Univariate convergence diagnostic, Brooks and Gelman (1998):')
    for j=1:npar
      fprintf('    Parameter %d...  ',j);
      for b = 1:nblck
    startline = 0;
    for n = 0:nfile
      instr = [fname_ '_mh' int2str(n) '_blck' int2str(b)];
      eval(['load ' instr]);
      clear logpo2 post2;
      tmp((b-1)*nsim+startline+1:(b-1)*nsim+cumulated_number_of_simulations_per_file(n+1),1) = x2(:,j);
      clear x2;
      startline = startline+number_of_simulations_per_file(n+1);
    end 
      end
      tmp(:,2) = kron(transpose(1:nblck),ones(nsim,1));
      tmp(:,3) = kron(ones(nblck,1),transpose(time)); 
      tmp = sortrows(tmp,1);
      ligne   = 0;
      for iter  = origin:step_size:nsim
    ligne = ligne+1;
    linea = ceil(0.5*iter);
    n     = iter-linea+1;
    cinf  = round(n*ALPHA/2);
    csup  = round(n*(1-ALPHA/2));
    CINF  = round(nblck*n*ALPHA/2);
    CSUP  = round(nblck*n*(1-ALPHA/2));
    temp  = tmp(find((tmp(:,3)>=linea) & (tmp(:,3)<=iter)),1:2);
    UDIAG(ligne,1,j) = temp(CSUP,1)-temp(CINF,1);
    moyenne = mean(temp(:,1));%% Pooled mean.
    UDIAG(ligne,3,j) = sum((temp(:,1)-moyenne).^2)/(nblck*n-1);
    UDIAG(ligne,5,j) = sum(abs(temp(:,1)-moyenne).^3)/(nblck*n-1);
    for i=1:nblck
      pmet = temp(find(temp(:,2)==i));
      UDIAG(ligne,2,j) = UDIAG(ligne,2,j) + pmet(csup,1)-pmet(cinf,1);
      moyenne = mean(pmet,1); %% Within mean. 
      UDIAG(ligne,4,j) = UDIAG(ligne,4,j) + sum((pmet(:,1)-moyenne).^2)/(n-1);
      UDIAG(ligne,6,j) = UDIAG(ligne,6,j) + sum(abs(pmet(:,1)-moyenne).^3)/(n-1);
    end
      end
      fprintf('Done! \n');
    end
    UDIAG(:,[2 4 6],:) = UDIAG(:,[2 4 6],:)/nblck;
    disp(' ')
    clear pmet temp moyenne CSUP CINF csup cinf n linea iter tmp;    
    pages = floor(npar/3);
    k = 0;  
    for i = 1:pages
      h = figure('Name','MCMC univariate diagnostic (Brooks and Gelman,1998)');
      boxplot = 1;
      if TeX
    NAMES = [];
    TEXNAMES = [];
      end
      for j = 1:3 % Loop over parameters
    k = k+1;
    [nam,namtex] = get_the_name(k,TeX);
    for crit = 1:3% Loop over criteria
      if crit == 1
        plt1 = UDIAG(:,1,k);
        plt2 = UDIAG(:,2,k);
        namnam  = [nam , ' (Interval)']; 
      elseif crit == 2
        plt1 = UDIAG(:,3,k);
        plt2 = UDIAG(:,4,k);
        namnam  = [nam , ' (m2)'];
      elseif crit == 3    
        plt1 = UDIAG(:,5,k);
        plt2 = UDIAG(:,6,k);
        namnam  = [nam , ' (m3)'];
      end
      if TeX
        NAMES = strvcat(NAMES,deblank(namnam));
        TEXNAMES = strvcat(TEXNAMES,deblank(namtex));
      end
      subplot(3,3,boxplot);
      plot(xx,plt1,'-b');     % Pooled
      hold on;
      plot(xx,plt2,'-r');     % Within (mean)
      hold off;
      xlim([xx(1) xx(number_of_lines)])
      title(namnam,'Interpreter','none')
      boxplot = boxplot + 1;
    end
      end
      eval(['print -depsc2 ' fname_ '_udiag' int2str(i)]);
      eval(['print -dpdf ' fname_ '_udiag' int2str(i)]);
      saveas(h,[fname_ '_udiag' int2str(i) '.fig']);
      if options_.nograph, close(h), end
      if TeX
    fprintf(fidTeX,'\\begin{figure}[H]\n');
    for jj = 1:size(NAMES,1)
      fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TEXNAMES(jj,:)));
    end    
    fprintf(fidTeX,'\\centering \n');
    fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_udiag%s}\n',fname_,int2str(i));
    fprintf(fidTeX,'\\caption{Univariate convergence diagnostics for the Metropolis-Hastings.\n');
    fprintf(fidTeX,'The first, second and third columns are respectively the criteria based on\n');
    fprintf(fidTeX,'the eighty percent interval, the second and third moments.}');
    fprintf(fidTeX,'\\label{Fig:UnivariateDiagnostics:%s}\n',int2str(i));
    fprintf(fidTeX,'\\end{figure}\n');
    fprintf(fidTeX,'\n');
      end
    end
    reste = npar-k;
    if reste
      if reste == 1
    nr = 3;
    nc = 1;
      elseif reste == 2;
    nr = 2;
    nc = 3;
      end
      if TeX
    NAMES = [];
    TEXNAMES = [];
      end
      h = figure('Name','MCMC univariate diagnostic (Brooks and Gelman, 1998)');
      boxplot = 1;
      for j = 1:reste
    k = k+1;
    [nam,namtex] = get_the_name(k,TeX);
    for crit = 1:3
      if crit == 1
        plt1 = UDIAG(:,1,k);
        plt2 = UDIAG(:,2,k);
        namnam  = [nam , ' (Interval)']; 
      elseif crit == 2
        plt1 = UDIAG(:,3,k);
        plt2 = UDIAG(:,4,k);
        namnam  = [nam , ' (m2)'];
      elseif crit == 3    
        plt1 = UDIAG(:,5,k);
        plt2 = UDIAG(:,6,k);
        namnam  = [nam , ' (m3)'];
      end
      if TeX
        NAMES = strvcat(NAMES,deblank(namnam));
        TEXNAMES = strvcat(TEXNAMES,deblank(namtex));
      end
      subplot(nr,nc,boxplot);
      plot(xx,plt1,'-b');                   % Pooled
      hold on;
      plot(xx,plt2,'-r');                   % Within (mean)
      hold off;
      xlim([xx(1) xx(number_of_lines)]);
      title(namnam,'Interpreter','none');
      boxplot = boxplot + 1;
    end
      end
      eval(['print -depsc2 ' fname_ '_udiag' int2str(pages+1)]);
      eval(['print -dpdf ' fname_ '_udiag' int2str(pages+1)]);
      saveas(h,[fname_ '_udiag' int2str(pages+1) '.fig']);
      if options_.nograph, close(h), end
      if TeX
    fprintf(fidTeX,'\\begin{figure}[H]\n');
    for jj = 1:size(NAMES,1);
      fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TEXNAMES(jj,:)));
    end    
    fprintf(fidTeX,'\\centering \n');
    fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_udiag%s}\n',fname_,int2str(pages+1));
    if reste == 2
      fprintf(fidTeX,'\\caption{Univariate convergence diagnostics for the Metropolis-Hastings.\n');
      fprintf(fidTeX,'The first, second and third columns are respectively the criteria based on\n');
      fprintf(fidTeX,'the eighty percent interval, the second and third moments.}');
    elseif reste == 1
      fprintf(fidTeX,'\\caption{Univariate convergence diagnostics for the Metropolis-Hastings.\n');
      fprintf(fidTeX,'The first, second and third rows are respectively the criteria based on\n');
      fprintf(fidTeX,'the eighty percent interval, the second and third moments.}');
    end
    fprintf(fidTeX,'\\label{Fig:UnivariateDiagnostics:%s}\n',int2str(pages+1));
    fprintf(fidTeX,'\\end{figure}\n');
    fprintf(fidTeX,'\n');
    fprintf(fidTeX,'% End Of TeX file.');
    fclose(fidTeX);
      end
    end % if reste > 0
    clear UDIAG;
    %%
    %% Multivariate diagnostic.
    %%
    if TeX
      fidTeX = fopen([fname_ '_MultivariateDiagnostics.TeX'],'w');
      fprintf(fidTeX,'%% TeX eps-loader file generated by metropolis.m (Dynare).\n');
      fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
      fprintf(fidTeX,' \n');
      NAMES = [];
    end
    tmp = zeros(nsim*nblck,3);
    MDIAG = zeros(number_of_lines,6);
    for b = 1:nblck
      startline = 0;
      for n = 0:nfile
    instr = [fname_ '_mh' int2str(n) '_blck' int2str(b)];
    eval(['load ' instr]);
    clear x2 post2;
    tmp((b-1)*nsim+startline+1:(b-1)*nsim+cumulated_number_of_simulations_per_file(n+1),1) = logpo2;
    startline = startline+number_of_simulations_per_file(n+1);
      end   
    end
    clear logpo2;
    tmp(:,2) = kron(transpose(1:nblck),ones(nsim,1));
    tmp(:,3) = kron(ones(nblck,1),transpose(time)); 
    tmp = sortrows(tmp,1);
    ligne   = 0;
    for iter  = origin:step_size:nsim
      ligne = ligne+1;
      linea = ceil(0.5*iter);
      n     = iter-linea+1;
      cinf  = round(n*ALPHA/2);
      csup  = round(n*(1-ALPHA/2));
      CINF  = round(nblck*n*ALPHA/2);
      CSUP  = round(nblck*n*(1-ALPHA/2));
      temp  = tmp(find((tmp(:,3)>=linea) & (tmp(:,3)<=iter)),1:2);
      MDIAG(ligne,1) = temp(CSUP,1)-temp(CINF,1);
      moyenne = mean(temp(:,1));%% Pooled mean.
      MDIAG(ligne,3) = sum((temp(:,1)-moyenne).^2)/(nblck*n-1);
      MDIAG(ligne,5) = sum(abs(temp(:,1)-moyenne).^3)/(nblck*n-1);
      for i=1:nblck
	pmet = temp(find(temp(:,2)==i));
	MDIAG(ligne,2) = MDIAG(ligne,2) + pmet(csup,1)-pmet(cinf,1);
	moyenne = mean(pmet,1); %% Within mean. 
	MDIAG(ligne,4) = MDIAG(ligne,4) + sum((pmet(:,1)-moyenne).^2)/(n-1);
	MDIAG(ligne,6) = MDIAG(ligne,6) + sum(abs(pmet(:,1)-moyenne).^3)/(n-1);
      end
    end
    MDIAG(:,[2 4 6],:) = MDIAG(:,[2 4 6],:)/nblck;  
    h = figure('Name','Multivatiate diagnostic');
    boxplot = 1;
    for crit = 1:3
      if crit == 1
	plt1 = MDIAG(:,1);
	plt2 = MDIAG(:,2);
	namnam  = 'Interval'; 
      elseif crit == 2
	plt1 = MDIAG(:,3);
	plt2 = MDIAG(:,4);
	namnam  = 'm2';
      elseif crit == 3    
	plt1 = MDIAG(:,5);
	plt2 = MDIAG(:,6);
	namnam  = 'm3';
      end
      if TeX
	NAMES = strvcat(NAMES,namnam);
      end
      subplot(3,1,boxplot);
      plot(xx,plt1,'-b');  % Pooled
      hold on
      plot(xx,plt2,'-r');  % Within (mean)
      hold off
      xlim([xx(1) xx(number_of_lines)])
      title(namnam,'Interpreter','none');
      boxplot = boxplot + 1;
    end
    eval(['print -depsc2 ' fname_ '_mdiag']);
    eval(['print -dpdf ' fname_ '_mdiag']);
    saveas(h,[fname_ '_mdiag.fig']);
    if options_.nograph, close(h), end
    if TeX
      fprintf(fidTeX,'\\begin{figure}[H]\n');
      for jj = 1:3
	fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),' ');
      end    
      fprintf(fidTeX,'\\centering \n');
      fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_mdiag}\n',fname_);
      fprintf(fidTeX,'\\caption{Multivariate convergence diagnostics for the Metropolis-Hastings.\n');
      fprintf(fidTeX,'The first, second and third rows are respectively the criteria based on\n');
      fprintf(fidTeX,'the eighty percent interval, the second and third moments. The different \n');
      fprintf(fidTeX,'parameters are aggregated using the posterior kernel.}');
      fprintf(fidTeX,'\\label{Fig:MultivariateDiagnostics}\n');
      fprintf(fidTeX,'\\end{figure}\n');
      fprintf(fidTeX,'\n');
      fprintf(fidTeX,'% End Of TeX file.');
      fclose(fidTeX);
    end
  end % End of if ~options_.nodiagnostic
  %%
  %% Now i discard some simulations...
  %%
  trun = cumulated_number_of_simulations_per_file(nfile+1);
  irun = floor(options_.mh_drop*trun)+1;
  ffil = 0;       % The first MH file we have to read...
  ifil = irun;    % and the first line we have to read in this file.
  for ffil = 0:nfile
    if irun <= cumulated_number_of_simulations_per_file(ffil+1)
      break
    end
    ifil = ifil-number_of_simulations_per_file(ffil+1);
  end
  trun = trun-irun+1;
  fprintf('MH: I''ll use mh-files %d to %d.\n',ffil,nfile);
  fprintf('MH: In mh-file number %d i''ll start at line %d.\n',ffil,ifil);
  fprintf('MH: Finally the total number of simulations is %d.\n',trun);
  disp(' ');
  %
  %%
  %%%
  %%%%
  %%%%% Modified harmonic mean
  %%%%
  %%%
  %%
  %
  fprintf('MH: I''m computing the posterior mean... ');
  MU = zeros(1,npar);
  lpost_mode = -Inf;
  for  b = 1:nblck
    if nblck > 1
      instr = [fname_ '_mh' int2str(ffil) '_blck' int2str(b)];
    else
      instr = [fname_ '_mh' int2str(ffil)];
    end
    eval(['load ' instr]); clear post2;
    MU(1,:) = MU(1,:) + sum(x2(ifil:end,:),1);
    lpost_mode = max(lpost_mode,max(logpo2(ifil:end,1)));
  end
  for n = ffil+1:nfile
    for b = 1:nblck
      if nblck > 1
	instr = [fname_ '_mh' int2str(n) '_blck' int2str(b)];
      else
	instr = [fname_ '_mh' int2str(n)];
      end
      eval(['load ' instr]);
      clear post2;
      MU(1,:) = MU(1,:) + sum(x2,1);
      lpost_mode = max(lpost_mode,max(logpo2));
    end
  end
  clear x2 logpo2;
  MU = MU/(trun*nblck);
  fprintf(' Done!\n');
  fprintf('MH: I''m computing the posterior covariance matrix... ');
  SIGMA = zeros(npar,npar);
  for b = 1:nblck
    if nblck > 1
      instr = [fname_ '_mh' int2str(ffil) '_blck' int2str(b)];
    else
      instr = [fname_ '_mh' int2str(ffil)];
    end  
    eval(['load ' instr]);
    clear post2 logpo2;
    SIGMA = SIGMA + transpose(x2(ifil:end,:)-ones(size(x2(ifil:end,:),1),1)*MU)*...
	    (x2(ifil:end,:)-ones(size(x2(ifil:end,:),1),1)*MU);
  end               
  for n = ffil+1:nfile
    for b = 1:nblck
      if nblck > 1
	instr = [fname_ '_mh' int2str(n) '_blck' int2str(b)];
      else
	instr = [fname_ '_mh' int2str(n)];
      end
      eval(['load ' instr]);
      clear post2 logpo2;
      SIGMA = SIGMA + transpose(x2-ones(size(x2,1),1)*MU)*(x2-ones(size(x2,1),1)*MU);
    end             
  end
  clear x2;
  SIGMA =  SIGMA/(trun*nblck);%<=== Variance of the parameters (ok!)
  fprintf(' Done!\n');
  disp(' ');
  disp('MH: I''m computing the posterior log marginale density (modified harmonic mean)... ');
  detSIGMA = det(SIGMA);
  invSIGMA = inv(SIGMA);
  marginal = zeros(9,2);
  linee = 0;
  check_coverage  = 1;
  increase        = 1;
  while check_coverage
    for p = 0.1:0.1:0.9;
      critval = qchisq(p,npar);
      tmp = 0;
      for k = ffil:nfile
    inst = [fname_ '_mh' int2str(k)];
    if k == ffil
      i1 = ifil;
    else
      i1 = 1;
    end
    EndOfFile = number_of_simulations_per_file(k+1);
    for b=1:nblck
      if nblck > 1
        instr = [inst '_blck' int2str(b)];
      else
        instr = inst;
      end  
      load(instr,'x2','logpo2');
      for i = i1:EndOfFile
        deviation  = (x2(i,:)-MU)*invSIGMA*(x2(i,:)-MU)';
        if deviation <= critval
          lftheta = -log(p)-(npar*log(2*pi)+log(detSIGMA)+deviation)/2;
          tmp = tmp + exp(lftheta - logpo2(i)+lpost_mode);
        end
      end
    end 
      end
      clear x2 logpo2;
      linee = linee + 1;    
      marginal(linee,:) = [p,lpost_mode-log(tmp/(trun*nblck))];
    end
    if abs((marginal(9,2)-marginal(1,2))/marginal(9,2)) > 0.01 | isinf(marginal(1,2))
      if increase == 1
	disp('MH: The support of the weighting density function is not large enough...')
	disp('MH: I increase the variance of this distribution.')
	increase = 1.2*increase;
	invSIGMA = inv(SIGMA*increase);
	detSIGMA = det(SIGMA*increase);
	linee    = 0;   
      else
	disp('MH: Let me try again.')
	increase = 1.2*increase;
	invSIGMA = inv(SIGMA*increase);
	detSIGMA = det(SIGMA*increase);
	linee    = 0;
	if increase > 20
	  check_coverage = 0;
	  clear invSIGMA detSIGMA increase;
	  disp('MH: There''s probably a problem with the modified harmonic mean estimator.')    
	end    
      end    
    else
      check_coverage = 0;
      clear invSIGMA detSIGMA increase;
      disp('MH: Modified harmonic mean estimator, done!')
    end
  end
  %
  %%
  %%%
  %%%%
  %%%%% Highest Probability Intervals (coverage is given by options_.mh_conf_sig)
  %%%%
  %%%
  %%
  %
  disp(' ')
  fprintf('MH: I''m computing the Highest Probability Intervals... ');
  post_mean = transpose(MU);
  n = trun*nblck;
  n1    = round((1-options_.mh_conf_sig)*n);
  k = zeros(n1,1);
  tmp = zeros(n,1);
  if nblck == 1
    for i = 1:npar
      EndOfFile = number_of_simulations_per_file(ffil+1)-ifil+1;
      instr = [fname_ '_mh' int2str(ffil)];
      eval(['load ' instr]);
      clear post2 logpo2;
      tmp(1:EndOfFile) = x2(ifil:end,i);
      OldEndOfFile = EndOfFile;
      for f = ffil+1:nfile
	NewEndOfFile = number_of_simulations_per_file(f+1);
	instr = [fname_ '_mh' int2str(f)];
	eval(['load ' instr]);
	clear post2 logpo2;
	tmp(OldEndOfFile+1:OldEndOfFile+NewEndOfFile) = x2(:,i);
	OldEndOfFile = OldEndOfFile + NewEndOfFile;
      end
      clear x2;
      tmp = sort(tmp);
      j2 = n-n1;
      for j1 = 1:n1
	k(j1) = tmp(j2)-tmp(j1);
	j2 = j2 + 1;
      end
      [kmin,k1] = min(k);
      min_interval(i,:) = [tmp(k1) tmp(k1)+kmin];
    end
    clear tmp;
  else
    for i = 1:npar
      EndOfFile = number_of_simulations_per_file(ffil+1)-ifil+1;
      NewStartLine = 0;
      inst = [fname_ '_mh' int2str(ffil)];
      for b = 1:nblck
	instr = [inst '_blck' int2str(b)];
	eval(['load ' instr]);
	clear post2 logpo2;
	tmp(NewStartLine+1:NewStartLine+EndOfFile,1) = x2(ifil:end,i);
	NewStartLine = NewStartLine + EndOfFile;
      end
      for f = ffil+1:nfile
	EndOfFile = number_of_simulations_per_file(f+1);
	inst = [fname_ '_mh' int2str(f)];
	for B = 1:nblck
	  instr = [inst '_blck' int2str(b)];
	  eval(['load ' instr]);
	  clear post2 logpo2;
	  tmp(NewStartLine+1:NewStartLine+EndOfFile,1) = x2(:,i);
	  NewStartLine = NewStartLine + EndOfFile;
	end
      end
      clear x2;
      tmp = sort(tmp);
      j2 = n-n1;
      for j1 = 1:n1
	k(j1) = tmp(j2)-tmp(j1);
	j2 = j2 + 1;
      end
      [kmin,k1] = min(k);
      min_interval(i,:) = [tmp(k1) tmp(k1)+kmin];
    end
    clear tmp;
  end
  fprintf(' Done!\n');
  %
  %%
  %%%
  %%%%
  %%%%% Print results
  %%%%
  %%%
  %%
  %%
  %% [1] On screen
  %%
  disp(' ');
  disp(' ')
  marginal
  disp(' ')
  disp(' ')
  disp('ESTIMATION RESULTS')
  disp(' ')
  disp(sprintf('Log data density is %f.',mean(marginal(:,2))))
  oo_.MarginalDensity.ModifiedHarmonicMean = mean(marginal(:,2));
  pnames=['     ';'beta ';'gamm ';'norm ';'invg ';'unif ';'invg2'];
  tit2 = sprintf('%10s %7s %10s %14s %4s %6s\n',' ','prior mean', ...
         'post. mean','conf. interval','prior','pstdev');
  ip = nvx+nvn+ncx+ncn+1;
  if np
    disp(' ')
    disp('parameters')
    disp(tit2)
    for i=1:np
      disp(sprintf('%12s %7.3f %8.4f %7.4f %7.4f %4s %6.4f', ...
           deblank(estim_params_.param_names(i,:)), ...
           bayestopt_.pmean(ip),post_mean(ip),min_interval(ip,:), ...
           pnames(bayestopt_.pshape(ip)+1,:), ...
           bayestopt_.pstdev(ip)));
      eval(['oo_.posterior_mean.parameters.' deblank(estim_params_.param_names(i,:)) ' = post_mean(ip);']);
      eval(['oo_.posterior_hpdinf.parameters.' deblank(estim_params_.param_names(i,:)) ' = min_interval(ip,1);']); 
      eval(['oo_.posterior_hpdsup.parameters.' deblank(estim_params_.param_names(i,:)) ' = min_interval(ip,2);']);
      ip = ip+1;
    end
  end
  if nvx
    ip = 1;
    disp(' ')
    disp('standard deviation of shocks')
    disp(tit2)
    for i=1:nvx
      k = estim_params_.var_exo(i,1);
      disp(sprintf('%12s %7.3f %8.4f %7.4f %7.4f %4s %6.4f', ...
           deblank(lgx_(k,:)),bayestopt_.pmean(ip),post_mean(ip), ...
           min_interval(ip,:),pnames(bayestopt_.pshape(ip)+1,:), ...
           bayestopt_.pstdev(ip))); 
      Sigma_e_(k,k) = post_mean(ip)*post_mean(ip);
      eval(['oo_.posterior_mean.shocks_std.' deblank(lgx_(k,:)) ' = post_mean(ip);']);
      eval(['oo_.posterior_hpdinf.shocks_std.' deblank(lgx_(k,:)) ' = min_interval(ip,1);']); 
      eval(['oo_.posterior_hpdsup.shocks_std.' deblank(lgx_(k,:)) ' = min_interval(ip,2);']);
      ip = ip+1;
    end
  end
  if nvn
    disp(' ')
    disp('standard deviation of measurement errors')
    disp(tit2)
    ip = nvx+1;
    for i=1:nvn
      disp(sprintf('%12s %7.3f %8.4f %7.4f %7.4f %4s %6.4f', ...
           deblank(options_.varobs(estim_params_.var_endo(i,1),:)),...
           bayestopt_.pmean(ip), ...
           post_mean(ip),min_interval(ip,:), ...
           pnames(bayestopt_.pshape(ip)+1,:), ...
           bayestopt_.pstdev(ip)));
      eval(['oo_.posterior_mean.measurement_errors_std.' deblank(options_.varobs(estim_params_.var_endo(i,1),:)) ' = post_mean(ip);']);
      eval(['oo_.posterior_hpdinf.measurement_errors_std.' deblank(options_.varobs(estim_params_.var_endo(i,1),:)) ' = min_interval(ip,1);']); 
      eval(['oo_.posterior_hpdsup.measurement_errors_std.' deblank(options_.varobs(estim_params_.var_endo(i,1),:)) ' = min_interval(ip,2);']);            
      ip = ip+1;
    end
  end
  if ncx
    disp(' ')
    disp('correlation of shocks')
    disp(tit2)
    ip = nvx+nvn+1;
    for i=1:ncx
      k1 = estim_params_.corrx(i,1);
      k2 = estim_params_.corrx(i,2);
      name = [deblank(lgx_(k1,:)) ',' deblank(lgx_(k2,:))];
      disp(sprintf('%12s %7.3f %8.4f %7.4f %7.4f %4s %6.4f', name, ...
           bayestopt_.pmean(ip),post_mean(ip),min_interval(ip,:), ...
           pnames(bayestopt_.pshape(ip)+1,:), ...
           bayestopt_.pstdev(ip)));
      eval(['oo_.posterior_mean.shocks_corr.' deblank(lgx_(k1,:)) '_' deblank(lgx_(k2,:)) ' = post_mean(ip);']);
      eval(['oo_.posterior_hpdinf.shocks_corr.' deblank(lgx_(k1,:)) '_' deblank(lgx_(k2,:)) ' = min_interval(ip,1);']); 
      eval(['oo_.posterior_hpdsup.shocks_corr.' deblank(lgx_(k1,:)) '_' deblank(lgx_(k2,:)) ' = min_interval(ip,2);']);      
      Sigma_e_(k1,k2) = post_mean(ip)*sqrt(Sigma_e_(k1,k1)*Sigma_e_(k2,k2));
      Sigma_e_(k2,k1) = Sigma_e_(k1,k2);
      ip = ip+1;
    end
  end
  if ncn
    disp(' ')
    disp('correlation of measurement errors')
    disp(tit2)
    ip = nvx+nvn+ncx+1;
    for i=1:ncn
      k1 = estim_params_.corrn(i,1);
      k2 = estim_params_.corrn(i,2);
      name = [deblank(lgy_(k1,:)) ',' deblank(lgy_(k2,:))];
      disp(sprintf('%12s %7.3f %8.4f %7.4f %7.4f %4s %6.4f', name, ...
           bayestopt_.pmean(ip),post_mean(ip),min_interval(ip,:), ...
           pnames(bayestopt_.pshape(ip)+1,:), ...
           bayestopt_.pstdev(ip))); 
      eval(['oo_.posterior_mean.measurement_errors_corr.' deblank(lgy_(k1,:)) '_' deblank(lgy_(k2,:)) ' = post_mean(ip);']);
      eval(['oo_.posterior_hpdinf.measurement_errors_corr.' deblank(lgy_(k1,:)) '_' deblank(lgy_(k2,:)) ' = min_interval(ip,1);']); 
      eval(['oo_.posterior_hpdsup.measurement_errors_corr.' deblank(lgy_(k1,:)) '_' deblank(lgy_(k2,:)) ' = min_interval(ip,2);']);      
      ip = ip+1;
    end
  end
  %%
  %% [1] In a TeX file
  %%
  if TeX 
    if np
      ip = nvx+nvn+ncx+ncn+1;
      fidTeX = fopen([fname_ '_MH_Posterior_1.TeX'],'w');
      fprintf(fidTeX,'%% TeX-table generated by metropolis.m (Dynare).\n');
      fprintf(fidTeX,'%% RESULTS FROM METROPOLIS HASTINGS (parameters)\n');
      fprintf(fidTeX,['%% ' datestr(now,0)]);
      fprintf(fidTeX,' \n');
      fprintf(fidTeX,' \n');
      fprintf(fidTeX,'{\\tiny \n');
      fprintf(fidTeX,'\\begin{table}[H]\n');
      fprintf(fidTeX,'\\centering\n');
      fprintf(fidTeX,'\\begin{tabular}{l|lccccc} \n');
      fprintf(fidTeX,'\\hline\\hline \\\\ \n');
      fprintf(fidTeX,'  & Prior distribution & Prior mean  & Prior s.d. & Post. mean & HPD inf & HPD sup  \\\\ \n');
      fprintf(fidTeX,'\\hline \\\\ \n');
      for i=1:np
	fprintf(fidTeX,' $%s$ & %s & %7.3f & %6.4f & %8.4f & %7.4f & %7.4f \\\\ \n',...
		deblank(estim_params_.tex(i,:)), ...
		deblank(pnames(bayestopt_.pshape(ip)+1,:)), ...
		bayestopt_.pmean(ip), ...
		bayestopt_.pstdev(ip), ...
		post_mean(ip), ...
		min_interval(ip,1), ...
		min_interval(ip,2));
	ip = ip+1;
      end   
      fprintf(fidTeX,'\\hline\\hline \n');
      fprintf(fidTeX,'\\end{tabular}\n ');    
      fprintf(fidTeX,'\\caption{Results from Metropolis Hastings (parameters)}\n ');
      fprintf(fidTeX,'\\label{Table:MhPosterior:1}\n');
      fprintf(fidTeX,'\\end{table}\n');
      fprintf(fidTeX,'} \n');
      fprintf(fidTeX,'%% End of TeX file.\n');
      fclose(fidTeX);
    end
    if nvx
      ip = 1;
      fidTeX = fopen([fname_ '_MH_Posterior_2.TeX'],'w');
      fprintf(fidTeX,'%% TeX-table generated by metropolis.m (Dynare).\n');
      fprintf(fidTeX,'%% RESULTS FROM METROPOLIS HASTINGS (standard deviation of structural shocks)\n');
      fprintf(fidTeX,['%% ' datestr(now,0)]);
      fprintf(fidTeX,' \n');
      fprintf(fidTeX,' \n');
      fprintf(fidTeX,'{\\tiny \n');
      fprintf(fidTeX,'\\begin{table}[H]\n');
      fprintf(fidTeX,'\\centering\n');
      fprintf(fidTeX,'\\begin{tabular}{l|lccccc} \n');
      fprintf(fidTeX,'\\hline\\hline \\\\ \n');
      fprintf(fidTeX,'  & Prior distribution & Prior mean  & Prior s.d. & Post. mean & HPD inf & HPD sup  \\\\ \n');
      fprintf(fidTeX,'\\hline \\\\ \n');
      for i=1:nvx
	k = estim_params_.var_exo(i,1);
	fprintf(fidTeX,' $%s$ & %s & %7.3f & %6.4f & %8.4f & %7.4f & %7.4f \\\\ \n',...
		deblank(lgx_TeX_(k,:)),...
		deblank(pnames(bayestopt_.pshape(ip)+1,:)), ...
		bayestopt_.pmean(ip), ...
		bayestopt_.pstdev(ip), ...
		post_mean(ip), ...
		min_interval(ip,1), ...
		min_interval(ip,1));
	ip = ip+1;
      end
      fprintf(fidTeX,'\\hline\\hline \n');
      fprintf(fidTeX,'\\end{tabular}\n ');    
      fprintf(fidTeX,'\\caption{Results from Metropolis Hastings (standard deviation of structural shocks)}\n ');
      fprintf(fidTeX,'\\label{Table:MhPosterior:2}\n');
      fprintf(fidTeX,'\\end{table}\n');
      fprintf(fidTeX,'} \n');
      fprintf(fidTeX,'%% End of TeX file.\n');
      fclose(fidTeX);
    end
    if nvn
      ip = nvx+1;
      fidTeX = fopen([fname_ '_MH_Posterior_3.TeX'],'w');
      fprintf(fidTeX,'%% TeX-table generated by metropolis.m (Dynare).\n');
      fprintf(fidTeX,'%% RESULTS FROM METROPOLIS HASTINGS (standard deviation of measurement errors)\n');
      fprintf(fidTeX,['%% ' datestr(now,0)]);
      fprintf(fidTeX,' \n');
      fprintf(fidTeX,' \n');
      fprintf(fidTeX,'{\\tiny \n');
      fprintf(fidTeX,'\\begin{table}[H]\n');
      fprintf(fidTeX,'\\centering\n');
      fprintf(fidTeX,'\\begin{tabular}{l|lccccc} \n');
      fprintf(fidTeX,'\\hline\\hline \\\\ \n');
      fprintf(fidTeX,'  & Prior distribution & Prior mean  & Prior s.d. & Post. mean & HPD inf & HPD sup  \\\\ \n');
      fprintf(fidTeX,'\\hline \\\\ \n');
      for i=1:nvn
	fprintf(fidTeX,' $%s$ & %s & %7.3f & %6.4f & %8.4f & %7.4f & %7.4f \\\\ \n',...
		deblank(options_.varobs_TeX(estim_params_.var_endo(i,1),:)), ...
		deblank(pnames(bayestopt_.pshape(ip)+1,:)), ...
		bayestopt_.pmean(ip), ...
		bayestopt_.pstdev(ip), ...
		post_mean(ip), ...
		min_interval(ip,1), ...
		min_interval(ip,2));
	p = ip+1;
      end
      fprintf(fidTeX,'\\hline\\hline \n');
      fprintf(fidTeX,'\\end{tabular}\n ');    
      fprintf(fidTeX,'\\caption{Results from Metropolis Hastings (standard deviation of measurement errors)}\n ');
      fprintf(fidTeX,'\\label{Table:MhPosterior:3}\n');
      fprintf(fidTeX,'\\end{table}\n');
      fprintf(fidTeX,'} \n');
      fprintf(fidTeX,'%% End of TeX file.\n');
      fclose(fidTeX);
    end
    if ncx
      ip = nvx+nvn+1;
      fidTeX = fopen([fname_ '_MH_Posterior_4.TeX'],'w');
      fprintf(fidTeX,'%% TeX-table generated by metropolis.m (Dynare).\n');
      fprintf(fidTeX,'%% RESULTS FROM METROPOLIS HASTINGS (correlation of structural shocks)\n');
      fprintf(fidTeX,['%% ' datestr(now,0)]);
      fprintf(fidTeX,' \n');
      fprintf(fidTeX,' \n');
      fprintf(fidTeX,'{\\tiny \n');
      fprintf(fidTeX,'\\begin{table}[H]\n');
      fprintf(fidTeX,'\\centering\n');
      fprintf(fidTeX,'\\begin{tabular}{l|lccccc} \n');
      fprintf(fidTeX,'\\hline\\hline \\\\ \n');
      fprintf(fidTeX,'  & Prior distribution & Prior mean  & Prior s.d. & Post. mean & HPD inf & HPD sup  \\\\ \n');
      fprintf(fidTeX,'\\hline \\\\ \n');
      for i=1:ncx
	k1 = estim_params_.corrx(i,1);
	k2 = estim_params_.corrx(i,2);
	name = [deblank(lgx_TeX_(k1,:)) ',' deblank(lgx_TeX_(k2,:))];
	fprintf(fidTeX,' $%s$ & %s & %7.3f & %6.4f & %8.4f & %7.4f & %7.4f \\\\ \n',...
		name, ...
		deblank(pnames(bayestopt_.pshape(ip)+1,:)), ...
		bayestopt_.pmean(ip), ...
		bayestopt_.pstdev(ip), ...
		post_mean(ip), ...
		min_interval(ip,1), ...
		min_interval(ip,2));
	ip = ip+1;
      end
      fprintf(fidTeX,'\\hline\\hline \n');
      fprintf(fidTeX,'\\end{tabular}\n ');    
      fprintf(fidTeX,'\\caption{Results from Metropolis Hastings (correlation of structural shocks)}\n ');
      fprintf(fidTeX,'\\label{Table:MhPosterior:4}\n');
      fprintf(fidTeX,'\\end{table}\n');
      fprintf(fidTeX,'} \n');
      fprintf(fidTeX,'%% End of TeX file.\n');
      fclose(fidTeX);
    end
    if ncn
      ip = nvx+nvn+ncx+1;
      fidTeX = fopen([fname_ '_MH_Posterior_5.TeX'],'w');
      fprintf(fidTeX,'%% TeX-table generated by metropolis.m (Dynare).\n');
      fprintf(fidTeX,'%% RESULTS FROM METROPOLIS HASTINGS (correlation of measurement errors)\n');
      fprintf(fidTeX,['%% ' datestr(now,0)]);
      fprintf(fidTeX,' \n');
      fprintf(fidTeX,' \n');
      fprintf(fidTeX,'{\\tiny \n');
      fprintf(fidTeX,'\\begin{table}[H]\n');
      fprintf(fidTeX,'\\centering\n');
      fprintf(fidTeX,'\\begin{tabular}{l|lccccc} \n');
      fprintf(fidTeX,'\\hline\\hline \\\\ \n');
      fprintf(fidTeX,'  & Prior distribution & Prior mean  & Prior s.d. & Post. mean & HPD inf & HPD sup  \\\\ \n');
      fprintf(fidTeX,'\\hline \\\\ \n');
      for i=1:ncn
	k1 = estim_params_.corrn(i,1);
	k2 = estim_params_.corrn(i,2);
	name = [deblank(lgy_TeX_(k1,:)) ',' deblank(lgy_TeX_(k2,:))];
	fprintf(fidTeX,' $%s$ & %s & %7.3f & %6.4f & %8.4f & %7.4f & %7.4f \\\\ \n',...
		name, ...
		deblank(pnames(bayestopt_.pshape(ip)+1,:)), ...
		bayestopt_.pmean(ip), ...
		bayestopt_.pstdev(ip), ...
		post_mean(ip), ...
		min_interval(ip,1), ...
		min_interval(ip,2));
	ip = ip+1;
      end
      fprintf(fidTeX,'\\hline\\hline \n');
      fprintf(fidTeX,'\\end{tabular}\n ');    
      fprintf(fidTeX,'\\caption{Results from Metropolis Hastings (correlation of structural shocks)}\n ');
      fprintf(fidTeX,'\\label{Table:MhPosterior:5}\n');
      fprintf(fidTeX,'\\end{table}\n');
      fprintf(fidTeX,'} \n');
      fprintf(fidTeX,'%% End of TeX file.\n');
      fclose(fidTeX);
    end
  end % if TeX
  %                                                               
  %%                                                              
  %%%                                                             
  %%%%                                                            
  %%%%% Plot posterior distributions
  %%%%                                                            
  %%%                                                             
  %%                                                              
  %                                                               
  figurename = 'Priors and posteriors';
  if TeX    
    fidTeX = fopen([fname_ '_PriorsAndPosteriors.TeX'],'w');
    fprintf(fidTeX,'%% TeX eps-loader file generated by metropolis.m (Dynare).\n');
    fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
    fprintf(fidTeX,' \n');
  end
  [nbplt,nr,nc,lr,lc,nstar] = pltorg(npar);
  if nbplt == 1
    h1 = figure('Name',figurename);
    if TeX
      TeXNAMES = [];
    end    
    NAMES    = []; 
    for i=1:npar
      [borneinf,bornesup,x1,x2,f1,f2,top,nam,texnam] = ...
      posterior_distribution(i,nfile,ffil,ifil,...
                 nblck,n,number_of_simulations_per_file,TeX);
      eval(['oo_.posterior_density.' deblank(nam) ' = [x1,f1];']);
      eval(['oo_.prior_density.' deblank(nam) ' = [x2,f2];']); 
      if TeX
	TeXNAMES = strvcat(TeXNAMES,texnam);
      end    
      NAMES = strvcat(NAMES,nam);
      subplot(nr,nc,i);
      hh = plot(x2,f2,'-k','linewidth',2);
      set(hh,'color',[0.7 0.7 0.7]);
      hold on;
      plot(x1,f1,'-k','linewidth',2);
      plot( [xparam1(i) xparam1(i)], [0,1.1*top], '--g', 'linewidth', 2);
      box on;
      axis([borneinf bornesup 0 1.1*top]);
      title(nam,'Interpreter','none');
      hold off;
      drawnow
    end
    eval(['print -depsc2 ' fname_ '_PriorsAndPosteriors' int2str(1)]);
    eval(['print -dpdf ' fname_ '_PriorsAndPosteriors' int2str(1)]);
    saveas(h1,[fname_ '_PriorsAndPosteriors' int2str(1) '.fig']);
    if TeX
      fprintf(fidTeX,'\\begin{figure}[H]\n');
      for jj = 1:npar
	fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
      end    
      fprintf(fidTeX,'\\centering\n');
      fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_PriorsAndPosteriors%s}\n',fname_,int2str(1));
      fprintf(fidTeX,'\\caption{Priors and posteriors.}');
      fprintf(fidTeX,'\\label{Fig:PriorsAndPosteriors:%s}\n',int2str(1));
      fprintf(fidTeX,'\\end{figure}\n');
      fprintf(fidTeX,' \n');
      fprintf(fidTeX,'%% End of TeX file.\n');
      fclose(fidTeX);
    end
    if options_.nograph, close(h1), end
  else
    for plt = 1:nbplt-1
      hplt = figure('Name',figurename);
      if TeX
	TeXNAMES = [];
      end    
      NAMES    = []; 
      for index=1:nstar
	names = [];
	i = (plt-1)*nstar + index;
	[borneinf,bornesup,x1,x2,f1,f2,top,nam,texnam] = ...
	    posterior_distribution(i,nfile,ffil,ifil,...
				   nblck,n,number_of_simulations_per_file,TeX);
	eval(['oo_.posterior_density.' deblank(nam) ' = [x1,f1];']);
	eval(['oo_.prior_density.' deblank(nam) ' = [x2,f2];']);                     
	if TeX
	  TeXNAMES = strvcat(TeXNAMES,texnam);
	end    
	NAMES = strvcat(NAMES,nam);
	subplot(nr,nc,index);
	hh = plot(x2,f2,'-k','linewidth',2);
	set(hh,'color',[0.7 0.7 0.7]);
	hold on;
	plot(x1,f1,'-k','linewidth',2);
	plot( [xparam1(i) xparam1(i)], [0,1.1*top], '--g', 'linewidth', 2);
	box on;
	axis([borneinf bornesup 0 1.1*top]);
	title(nam,'Interpreter','none');
	hold off;
	drawnow;
      end  % index=1:nstar
      eval(['print -depsc2 ' fname_ '_PriorsAndPosteriors' int2str(plt)]);
      eval(['print -dpdf ' fname_ '_PriorsAndPosteriors' int2str(plt)]);
      saveas(hplt,[fname_ '_PriorsAndPosteriors' int2str(plt) '.fig']);
      if TeX
	fprintf(fidTeX,'\\begin{figure}[H]\n');
	for jj = 1:nstar
	  fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
	end    
	fprintf(fidTeX,'\\centering\n');
	fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_PriorsAndPosteriors%s}\n',fname_,int2str(plt));
	fprintf(fidTeX,'\\caption{Priors and posteriors.}');
	fprintf(fidTeX,'\\label{Fig:PriorsAndPosteriors:%s}\n',int2str(plt));
	fprintf(fidTeX,'\\end{figure}\n');
	fprintf(fidTeX,' \n');
      end    
      if options_.nograph, close(hplt), end
    end % plt = 1:nbplt-1
    hplt = figure('Name',figurename);
    if TeX
      TeXNAMES = [];
    end    
    NAMES    = []; 
    for index=1:npar-(nbplt-1)*nstar
      i = (nbplt-1)*nstar +  index;
      [borneinf,bornesup,x1,x2,f1,f2,top,nam,texnam] = ...
	  posterior_distribution(i,nfile,ffil,ifil,...
				 nblck,n,number_of_simulations_per_file,TeX);
      eval(['oo_.posterior_density.' deblank(nam) ' = [x1,f1];']);
      eval(['oo_.prior_density.' deblank(nam) ' = [x2,f2];']);             
      if TeX
	TeXNAMES = strvcat(TeXNAMES,texnam);
      end
      NAMES = strvcat(NAMES,nam);
      if lr
	subplot(lc,lr,index);
      else
	subplot(nr,nc,index);
      end    
      hh = plot(x2,f2,'-k','linewidth',2);
      set(hh,'color',[0.7 0.7 0.7]);
      hold on;
      plot(x1,f1,'-k','linewidth',2);
      plot( [xparam1(i) xparam1(i)], [0,1.1*top], '--g', 'linewidth', 2);
      box on;
      axis([borneinf bornesup 0 1.1*top]);
      title(nam,'Interpreter','none');
      hold off;
      drawnow;
    end  % index=1:npar-(nbplt-1)*nstar
    eval(['print -depsc2 ' fname_ '_PriorsAndPosteriors' int2str(nbplt)]);
    eval(['print -dpdf ' fname_ '_PriorsAndPosteriors' int2str(nbplt)]);
    saveas(hplt,[fname_ '_PriorsAndPosteriors' int2str(nbplt) '.fig']);
    if TeX
      fprintf(fidTeX,'\\begin{figure}[H]\n');
      for jj = 1:npar-(nbplt-1)*nstar
	fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
      end    
      fprintf(fidTeX,'\\centering\n');
      fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_PriorsAndPosteriors%s}\n',fname_,int2str(nbplt));
      fprintf(fidTeX,'\\caption{Priors and posteriors.}');
      fprintf(fidTeX,'\\label{Fig:PriorsAndPosteriors:%s}\n',int2str(nbplt));
      fprintf(fidTeX,'\\end{figure}\n');
      fprintf(fidTeX,' \n');
      fprintf(fidTeX,'%% End of TeX file.\n');
      fclose(fidTeX);
    end
    if options_.nograph, close(hplt), end
  end
  %                                                               
  %%                                                              
  %%%                                                             
  %%%%                                                            
  %%%%% Je (re)fais mes comptes... I should be able to skip this part (already done) 
  %%%%                                                            
  %%%                                                             
  %%                                                              
  %                                                               
  FLN = zeros(nfile-ffil+1,3);% Describes the number of lines in each file
  if nblck == 1
    instr1 = [fname_ '_mh'];
    instr2 = '';
  else  % I only consider draws from the first chain. This is correct
    % if and only if the metropolis-hastings did converge.
    instr1 = [fname_ '_mh'];
    instr2 = '_blck1';
  end   
  eval(['load ' instr1 int2str(ffil) instr2]);
  clear post2 x2;
  FLN(1,1) = ffil;                            % File number   
  FLN(1,2) = size(logpo2(ifil:end,:),1);      % Number of simulations in this file (density) 
  FLN(1,3) = FLN(1,2);                        % Cumulative Distribution Function
  if nfile-ffil+1>1
    linee = 1;
    for n = ffil+1:nfile
      linee = linee+1;
      instr = [instr1 int2str(n) instr2];
      eval(['load ' instr]);
      clear post2 x2;
      FLN(linee,1) = n;
      FLN(linee,2) = size(logpo2,1);
      FLN(linee,3) = FLN(linee-1,3) + FLN(linee,2);  
    end
    clear logpo2
    nruns = FLN(linee,3);
  else
    nruns = FLN(1,3);    
  end
  FLN(:,3) = FLN(:,3)/nruns;% I'm scaling the CDF
  nvar     = endo_nbr;
  B        = min(500,nruns);
  deciles = [round(0.1*B) ...
	     round(0.2*B)...
	     round(0.3*B)...
	     round(0.4*B)...
	     round(0.5*B)...
	     round(0.6*B)...
	     round(0.7*B)...
	     round(0.8*B)...
	     round(0.9*B)];
  %                                                               
  %%                                                              
  %%%                                                             
  %%%%                                                            
  %%%%% SDGE-based forecasts, smooth and filtered variables, IRFs and theoretical moments 
  %%%%                                                            
  %%%                                                             
  %%                                                              
  %                                                               
  if options_.forecast | options_.smoother | options_.filtered_vars
    deep = MU;
    subindx = subset();
    % [1] I delete some old files...    
    disp(' ')
    disp(' ')
    if options_.forecast
      files = eval(['dir(''' fname_ '_forecast*.mat'');']);
      if size(files,1)
	delete([fname_ '_forecast*.mat']);
	disp(['MH: Old ' fname_ '_forecast files deleted! '])
      end
    end
    if options_.smoother        
      files = eval(['dir(''' fname_ '_smooth*.mat'');']);
      if size(files,1)
	delete([fname_ '_smooth*.mat']);
	disp(['MH: Old ' fname_ '_smooth files deleted! '])
      end
      files = eval(['dir(''' fname_ '_innovation*.mat'');']);
      if size(files,1)
	delete([fname_ '_innovation*.mat']);
	disp(['MH: Old ' fname_ '_innovation files deleted! '])
      end
      files = eval(['dir(''' fname_ '_error*.mat'');']);
      if size(files,1)
	delete([fname_ '_error*.mat']);
	disp(['MH: Old ' fname_ '_error files deleted! '])
      end
    end
    if options_.filtered_vars
      files = eval(['dir(''' fname_ '_filter*.mat'');']);     
      if size(files,1)                                        
	delete([fname_ '_filter*.mat']);                    
	disp(['MH: Old ' fname_ '_filter files deleted! ']) 
      end                                                         
    end
    disp(' ')
    disp(' ')
    % [2] Initialization...    
    ex_      = zeros(horizon+xkmin_+xkmax_,exo_nbr);
    yyyy     = zeros(nvar,ykmin_);
    IdObs    = bayestopt_.mfys;
    if options_.forecast 
      if B <= MAX_nforc
	stock_forcst = zeros(options_.forecast+ykmin_,nvar,B);
	stock_forcst1 = zeros(options_.forecast+ykmin_,nvar,B);
      else
	stock_forcst = zeros(options_.forecast+ykmin_,nvar,MAX_nforc);
	stock_forcst1 = zeros(options_.forecast+ykmin_,nvar,MAX_nforc);
      end
    end 
    if options_.smoother
      if B <= MAX_nsmoo
	stock_smooth = zeros(endo_nbr,gend,B);
      else
	stock_smooth = zeros(endo_nbr,gend,MAX_nsmoo);
      end
      if B <= MAX_ninno 
	stock_innov  = zeros(exo_nbr,gend,B);
      else
	stock_innov  = zeros(exo_nbr,gend,MAX_ninno);
      end
      if nvn & B <= MAX_nerro
	%stock_error = zeros(gend,nvobs,B);
	stock_error = zeros(nvobs,gend,B);
      else nvn & B > MAX_nerro
	%stock_error = zeros(gend,nvobs,MAX_nerro);
	stock_error = zeros(nvobs,gend,MAX_nerro);
      end
    end
    if options_.filtered_vars
      if B <= MAX_nfilt
	stock_filter = zeros(endo_nbr,gend+1,B);
      else
	stock_filter = zeros(endo_nbr,gend+1,MAX_nfilt);
      end
    end
    h = waitbar(0,'SDGE model based forecasts...');
    % [3]   CoRe    
    % [3.1] First we consider the case with measurement error
    if nvn
      % [3.1.1] More than one _mh file 
      if nfile-ffil+1>1			
	if options_.forecast
	  sfil_forc = 0;
	  irun_forc = 0;  			
	end
	if options_.smoother
	  sfil_smoo = 0;
	  sfil_inno = 0;
	  sfil_erro = 0;
	  irun_smoo = 0;
	  irun_inno = 0;
	  irun_erro = 0;
	end
	if options_.filtered_vars
	  sfil_filt = 0;
	  irun_filt = 0;  			
	end
	% [3.1.1.1] Loop in the metropolis
	for b = 1:B;
	  if options_.forecast
	    irun_forc = irun_forc+1;
	  end
	  if options_.smoother
	    irun_smoo = irun_smoo+1;
	    irun_inno = irun_inno+1;
	    irun_erro = irun_erro+1;
	  end
	  if options_.filtered_vars
	    irun_filt = irun_filt+1;  			
	  end    			
	  % FIRST, I choose an _mh file (where the posterior distribution is stored)
	  choose_an_mh_file = rand;
	  mh_file_number = FLN(find(choose_an_mh_file>=FLN(:,3)),1);
	  if isempty(mh_file_number)
	    mh_file_number = ffil;
	  else    
	    mh_file_number = mh_file_number(1);
	  end    
	  eval(['load ' instr1 int2str(mh_file_number) instr2]);
	  clear post2 logpo2;
	  % SECOND, I choose a vector of structural parameters (a line in the _mh file) 
	  DEEP  = x2(floor(rand*FLN(find(mh_file_number == FLN(:,1)),2))+1,:);
	  deep(subindx) = DEEP(subindx);
	  % THIRD, I estimate the smooth and filtered variables. I need the smoothed variables
	  % to estimate the state of the model at the end of the sample. 
	  [atT,innov,obs_err,filtered_state_vector,ys,trend_coeff] = DsgeSmoother(transpose(deep),gend,data);
	  % FOURTH, smoothed and filtered variables are saved if needed
	  if options_.smoother
	    if irun_erro < MAX_nerro
	      stock_error(:,:,irun_erro) = obs_err;
	    else
	      stock_error(:,:,irun_erro) = obs_err;
	      sfil_erro = sfil_erro + 1;
	      instr = [fname_ '_error' int2str(sfil_erro) ' stock_error;'];
	      eval(['save ' instr]);
	      irun_erro = 0;
	      stock_error  = zeros(gend,nvobs,MAX_nerro);
	    end
	    if irun_smoo < MAX_nsmoo
	      stock_smooth(:,:,irun_smoo) = atT(1:endo_nbr,1:gend);
	    else
	      stock_smooth(:,:,irun_smoo) = atT(1:endo_nbr,1:gend);
	      sfil_smoo = sfil_smoo + 1;
	      instr = [fname_ '_smooth' int2str(sfil_smoo) ' stock_smooth;'];
	      eval(['save ' instr]);
	      irun_smoo = 0;
	      stock_smooth = ...
		  zeros(endo_nbr,gend,MAX_nsmoo);
	    end	
	    if irun_inno < MAX_ninno
	      stock_innov(:,:,irun_inno) = innov;
	    else
	      stock_innov(:,:,irun_inno) = innov;
	      sfil_inno = sfil_inno + 1;
	      instr = [fname_ '_innovation' int2str(sfil_inno) ' stock_innov;'];
	      eval(['save ' instr]);
	      irun_inno = 0;
	      stock_innov  = zeros(exo_nbr,gend,MAX_ninno);
	    end	
	  end
	  if options_.filtered_vars
	    if irun_filt < MAX_nfilt
	      stock_filter(:,:,irun_filt) = filtered_state_vector;
	    else
	      stock_filter(:,:,irun_filt) = filtered_state_vector;
	      sfil_filt = sfil_filt + 1;
	      instr = [fname_ '_filter' int2str(sfil_filt) ' stock_filter;'];
	      eval(['save ' instr]);
	      irun_filt = 0;
	      stock_filter = ...
		  zeros(endo_nbr,gend+1,MAX_nfilt);
	    end	
	  end    			
	  if options_.forecast
	    % FIFTH, I update variable dr_ 
	    dr_ = resol(ys_,0);
	    % SIXTH, I do and save the forecasts (for all the endogenous variables)
	    % The state of the economy at the end of the sample 
	    % depends on the structural parameters.	    
	    yyyy(:,1:ykmin_) = atT(1:endo_nbr,size(atT,2)-ykmin_+1:size(atT,2));
	    yf = forcst2a(yyyy,dr_,ex_);
	    if options_.prefilter == 1
	      yf(:,IdObs) = yf(:,IdObs)+repmat(bayestopt_.mean_varobs', ...
					       horizon+ykmin_,1);
	    end
	    yf(:,IdObs) = yf(:,IdObs)+(gend+[1-ykmin_:horizon]')*trend_coeff';
	    if options_.loglinear == 1
	      yf = yf+repmat(log(ys'),horizon+ykmin_,1);
	      yf = exp(yf);
	    else
	      yf = yf+repmat(ys',horizon+ykmin_,1);
	    end
	    stock_forcst(:,:,irun_forc) = yf;
	    yf1 = forcst2(yyyy,horizon,dr_,1);
	    if options_.prefilter == 1
	      yf1(:,IdObs,:) = yf1(:,IdObs,:)+ ...
		  repmat(bayestopt_.mean_varobs',[horizon+ykmin_,1,1]);
	    end
	    yf1(:,IdObs,:) = yf1(:,IdObs,:)+repmat((gend+[1-ykmin_:horizon]')* ...
		trend_coeff',[1,1,1]);
	    if options_.loglinear == 1
	      yf1 = yf1 + repmat(log(ys'),[horizon+ykmin_,1,1]);
	      yf1 = exp(yf1);
	    else
	      yf1 = yf1 + repmat(ys',[horizon+ykmin_,1,1]);
	    end
	    stock_forcst1(:,:,irun_forc) = yf1;
	    if irun_forc == MAX_nforc
	      sfil_forc = sfil_forc + 1;
	      save([fname_ '_forecast' int2str(sfil_forc)],'stock_forcst','stock_forcst1');
	      irun_forc = 0;
	      stock_forcst = zeros(horizon+ykmin_,nvar,MAX_nforc);
	      stock_forcst1 = zeros(horizon+ykmin_,nvar,MAX_nforc);
	    end
	  end		
	  waitbar(b/B,h);    
	end % of loop [3.1.1.1]
	if options_.smoother
	  if irun_smoo
	    stock_smooth = stock_smooth(:,:,1:irun_smoo);
	    sfil_smoo = sfil_smoo + 1;
	    instr = [fname_ '_smooth' int2str(sfil_smoo) ' stock_smooth;'];
	    eval(['save ' instr]);
	  end
	  clear stock_smooth;
	  if irun_inno
	    stock_innov = stock_innov(:,:,1:irun_inno);
	    sfil_inno = sfil_inno + 1;
	    instr = [fname_ '_innovation' int2str(sfil_inno) ' stock_innov;'];
	    eval(['save ' instr]);
	  end
	  clear stock_innov;
	  if irun_erro
	    stock_error = stock_error(:,:,1:irun_erro);
	    sfil_erro = sfil_erro + 1;
	    instr = [fname_ '_error' int2str(sfil_erro) ' stock_error;'];
	    eval(['save ' instr]);
	  end
	  clear stock_error;
	end
	if options_.forecast	
	  if irun_forc
	    stock_forcst = stock_forcst(:,:,1:irun_forc);  
	    stock_forcst1 = stock_forcst1(:,:,1:irun_forc);  
	    sfil_forc = sfil_forc + 1;
	    save([fname_ '_forecast' int2str(sfil_forc)],'stock_forcst','stock_forcst1');
	  end
	  clear stock_forcst stock_forcst1
	end
	if options_.filtered_vars	
	  if irun_filt
	    stock_filter = stock_filter(:,:,1:irun_filt);
	    sfil_filt = sfil_filt + 1;
	    instr = [fname_ '_filter' int2str(sfil_filt) ' stock_filter;'];
	    eval(['save ' instr]);
	  end
	  clear stock_filter;
	end
      else % [3.1.2] Just one _mh file
	if options_.forecast
	  sfil_forc = 0;
	  irun_forc = 0;  			
	end
	if options_.smoother
	  sfil_smoo = 0;
	  sfil_inno = 0;
	  sfil_erro = 0;
	  irun_smoo = 0;
	  irun_inno = 0;
	  irun_erro = 0;
	end
	if options_.filtered_vars
	  sfil_filt = 0;
	  irun_filt = 0;  			
	end
	eval(['load ' instr1 int2str(ffil) instr2]);
	NumberOfSimulations = length(logpo2);
	clear post2 logpo2;
	for b = 1:B;
	  if options_.forecast
	    irun_forc = irun_forc+1;
	  end
	  if options_.smoother
	    irun_smoo = irun_smoo+1;
	    irun_inno = irun_inno+1;
	    irun_erro = irun_erro+1;
	  end
	  if options_.filtered_vars
	    irun_filt = irun_filt+1;            
	  end
	  DEEP  = x2(floor(rand*NumberOfSimulations)+1,:); 
	  deep(subindx) = DEEP(subindx);
	  [atT,innov,obs_err,filtered_state_vector,ys,trend_coeff] = DsgeSmoother(transpose(deep),gend,data);
	  if options_.smoother
	    if irun_erro < MAX_nerro
	      stock_error(:,:,irun_erro) = obs_err;
	    else
	      stock_error(:,:,irun_erro) = obs_err;
	      sfil_erro = sfil_erro + 1;
	      instr = [fname_ '_error' int2str(sfil_erro) ' stock_error;'];
	      eval(['save ' instr]);
	      irun_erro = 0;
	      stock_error  = zeros(gend,nvobs,MAX_nerro);
	    end
	    if irun_smoo < MAX_nsmoo
	      stock_smooth(:,:,irun_smoo) = atT(1:endo_nbr,1:gend);
	    else
	      stock_smooth(:,:,irun_smoo) = atT(1:endo_nbr,1:gend);
	      sfil_smoo = sfil_smoo + 1;
	      instr = [fname_ '_smooth' int2str(sfil_smoo) ' stock_smooth;'];
	      eval(['save ' instr]);
	      irun_smoo = 0;
	      stock_smooth = ...
		  zeros(endo_nbr,gend,MAX_nsmoo);
	    end 
	    if irun_inno < MAX_ninno
	      stock_innov(:,:,irun_inno) = innov;
	    else
	      stock_innov(:,:,irun_inno) = innov;
	      sfil_inno = sfil_inno + 1;
	      instr = [fname_ '_innovation' int2str(sfil_inno) ' stock_innov;'];
	      eval(['save ' instr]);
	      irun_inno = 0;
	      stock_innov  = zeros(exo_nbr,gend,MAX_ninno);
	    end
	  end
	  if options_.filtered_vars
	    if irun_filt < MAX_nfilt                                             
	      stock_filter(:,:,irun_filt) = filtered_state_vector;             
	    else                                                                 
	      stock_filter(:,:,irun_filt) = filtered_state_vector;             
	      sfil_filt = sfil_filt + 1;                                       
	      instr = [fname_ '_filter' int2str(sfil_filt) ' stock_filter;'];  
	      eval(['save ' instr]);                                           
	      irun_filt = 0;                                                   
	      stock_filter = ...                                               
		  zeros(endo_nbr,gend+1,MAX_nfilt);     
	    end                                                                  
	  end
	  if options_.forecast
	    dr_ = resol(ys_,0);
	    yyyy(:,1:ykmin_) = atT(1:endo_nbr,size(atT,2)-ykmin_+1:size(atT,2));
	    yf = forcst2a(yyyy,dr_,ex_);
	    if options_.prefilter == 1
	      yf(:,IdObs) = yf(:,IdObs)+repmat(bayestopt_.mean_varobs', ...
					       horizon+ykmin_,1);
	    end
	    yf(:,IdObs) = yf(:,IdObs)+(gend+[1-ykmin_:horizon]')*trend_coeff';
	    if options_.loglinear == 1
	      yf = yf+repmat(log(ys'),horizon+ykmin_,1);
	      yf = exp(yf);
	    else
	      yf = yf+repmat(ys',horizon+ykmin_,1);
	    end
	    stock_forcst(:,:,irun_forc) = yf;
	    yf1 = forcst2(yyyy,horizon,dr_,1);
	    if options_.prefilter == 1
	      yf1(:,IdObs,:) = yf1(:,IdObs,:)+ ...
		  repmat(bayestopt_.mean_varobs',[horizon+ykmin_,1,1]);
	    end
	    yf1(:,IdObs,:) = yf1(:,IdObs,:)+repmat((gend+[1-ykmin_:horizon]')* ...
						   trend_coeff',[1,1,1]);
	    if options_.loglinear == 1
	      yf1 = yf1 + repmat(log(ys'),[horizon+ykmin_,1,1]);
	      yf1 = exp(yf1);
	    else
	      yf1 = yf1 + repmat(ys',[horizon+ykmin_,1,1]);
	    end
	    stock_forcst1(:,:,irun_forc) = yf1;
	    if irun_forc == MAX_nforc
	      sfil_forc = sfil_forc + 1;
	      save([fname_ '_forecast' int2str(sfil_forc)],'stock_forcst','stock_forcst1');
	      irun_forc = 0;
	      stock_forcst = zeros(horizon+ykmin_,nvar,MAX_nforc);
	      stock_forcst1 = zeros(horizon+ykmin_,nvar,MAX_nforc);
	    end
	  end   
	  waitbar(b/B,h);    
	end % of the loop over the metropolis simulations
	if options_.smoother
	  if irun_smoo
	    stock_smooth = stock_smooth(:,:,1:irun_smoo);
	    sfil_smoo = sfil_smoo + 1;
	    instr = [fname_ '_smooth' int2str(sfil_smoo) ' stock_smooth;'];
	    eval(['save ' instr]);
	  end
	  clear stock_smooth;
	  if irun_inno
	    stock_innov = stock_innov(:,:,1:irun_inno);
	    sfil_inno = sfil_inno + 1;
	    instr = [fname_ '_innovation' int2str(sfil_inno) ' stock_innov;'];
	    eval(['save ' instr]);
	  end
	  clear stock_innov;
	  if irun_erro
	    stock_error = stock_error(:,:,1:irun_erro);
	    sfil_erro = sfil_erro + 1;
	    instr = [fname_ '_error' int2str(sfil_erro) ' stock_error;'];
	    eval(['save ' instr]);
	  end
	  clear stock_error;
	end
	if options_.forecast    
	  if irun_forc
	    stock_forcst = stock_forcst(:,:,1:irun_forc);  
	    stock_forcst1 = stock_forcst1(:,:,1:irun_forc);  
	    sfil_forc = sfil_forc + 1;
	    save([fname_ '_forecast' int2str(sfil_forc)],'stock_forcst','stock_forcst1');
	  end
	  clear stock_forcst stock_forcst1
	end
	if options_.filtered_vars
	  if irun_filt
	    stock_filter = stock_filter(:,:,1:irun_filt);
	    sfil_filt = sfil_filt + 1;
	    instr = [fname_ '_filter' int2str(sfil_filt) ' stock_filter;'];
	    eval(['save ' instr]);
	  end
	  clear stock_filter;           
	end
      end
    else % [3.2]    Second we consider the case without measurement error
      if nfile-ffil+1>1
	if options_.forecast
	  sfil_forc = 0;
	  irun_forc = 0;  			
	end
	if options_.smoother
	  sfil_smoo = 0;
	  sfil_inno = 0;
	  sfil_erro = 0;
	  irun_smoo = 0;
	  irun_inno = 0;
	  irun_erro = 0;
	end
	if options_.filtered_vars
	  sfil_filt = 0;
	  irun_filt = 0;  			
	end
	for b = 1:B;
	  if options_.forecast
	    irun_forc = irun_forc+1;
	  end
	  if options_.smoother
	    irun_smoo = irun_smoo+1;
	    irun_inno = irun_inno+1;
	    irun_erro = irun_erro+1;
	  end
	  if options_.filtered_vars
	    irun_filt = irun_filt+1;  			
	  end	    
	  choose_an_mh_file = rand;
	  mh_file_number = FLN(find(choose_an_mh_file>=FLN(:,3)),1);
	  if isempty(mh_file_number)
	    mh_file_number = ffil;
	  else    
	    mh_file_number = mh_file_number(1);
	  end    
	  eval(['load ' instr1 int2str(mh_file_number) instr2]);
	  clear post2 logpo2;
	  DEEP  = x2(floor(rand*FLN(find(mh_file_number == FLN(:,1)),2))+1,:);
	  deep(subindx) = DEEP(subindx);
	  [atT,innov,obs_err,filtered_state_vector,ys,trend_coeff] = DsgeSmoother(transpose(deep),gend,data);
	  if options_.smoother
	    %if irun_erro < MAX_nerro
	    %	stock_error(:,:,irun_erro) = obs_err;
	    %else
	    %	stock_error(:,:,irun_erro) = obs_err;
	    %	instr = [fname_ '_error' int2str(sfil_erro) ' stock_error;'];
	    %	eval(['save ' instr]);
	    %	sfil_erro = sfil_erro + 1;
	    %	irun_erro = 0;
	    %	stock_error  = zeros(gend,nvobs,MAX_nerro);
	    %end
	    if irun_smoo < MAX_nsmoo
	      stock_smooth(:,:,irun_smoo) = atT(1:endo_nbr,1:gend);
	      if options_.prefilter == 1
		stock_smooth(bayestopt_.mf,:,irun_smoo) = atT(bayestopt_.mf,:)+repmat(bayestopt_.mean_varobs',1,gend);
	      elseif options_.loglinear == 1
		stock_smooth(bayestopt_.mf,:,irun_smoo) = atT(bayestopt_.mf,:)+...
		    repmat(log(ys(bayestopt_.mfys)),1,gend)+...
		     trend_coeff*[1:gend];
	      else
		stock_smooth(bayestopt_.mf,:,irun_smoo) = atT(bayestopt_.mf,:)+...
		    repmat(ys(bayestopt_.mfys),1,gend)+...
		     trend_coeff*[1:gend];
	      end
	    else
	      stock_smooth(:,:,irun_smoo) = atT(1:endo_nbr,1:gend);
	      if options_.prefilter == 1
		stock_smooth(bayestopt_.mf,:,irun_smoo) = atT(bayestopt_.mf,:)+repmat(bayestopt_.mean_varobs',1,gend);
	      elseif options_.loglinear == 1
		stock_smooth(bayestopt_.mf,:,irun_smoo) = atT(bayestopt_.mf,:)+...
		    repmat(log(ys(bayestopt_.mfys)),1,gend)+...
		     trend_coeff*[1:gend];
	      else
		stock_smooth(bayestopt_.mf,:,irun_smoo) = atT(bayestopt_.mf,:)+...
		    repmat(ys(bayestopt_.mfys),1,gend)+...
		     trend_coeff*[1:gend];
	      end
	      sfil_smoo = sfil_smoo + 1;
	      instr = [fname_ '_smooth' int2str(sfil_smoo) ' stock_smooth;'];
	      eval(['save ' instr]);
	      irun_smoo = 0;
	      stock_smooth = ...
		zeros(endo_nbr,gend,MAX_nsmoo);
	    end	
	    if irun_inno < MAX_ninno
	      stock_innov(:,:,irun_inno) = innov;
	    else
	      stock_innov(:,:,irun_inno) = innov;
	      sfil_inno = sfil_inno + 1;
	      instr = [fname_ '_innovation' int2str(sfil_inno) ' stock_innov;'];
	      eval(['save ' instr]);
	      irun_inno = 0;
	      stock_innov  = zeros(exo_nbr,gend,MAX_ninno);
	    end	
	  end
	  if options_.filtered_vars
	    if irun_filt < MAX_nfilt
	      stock_filter(:,:,irun_filt) = filtered_state_vector;
	    else
	      stock_filter(:,:,irun_filt) = filtered_state_vector;
	      sfil_filt = sfil_filt + 1;
	      instr = [fname_ '_filter' int2str(sfil_filt) ' stock_filter;'];
	      eval(['save ' instr]);
	      irun_filt = 0;
	      stock_filter = ...
		 zeros(endo_nbr,gend+1,MAX_nfilt);
	    end
	  end
	  if options_.forecast	    
	    dr_ = resol(ys_,0);
	    for j = 1:nvar 
	      yyyy(j,1:ykmin_) = atT(j,size(atT,2)-ykmin_+1:size(atT,2));
	    end
	    yf = forcst2a(yyyy,dr_,ex_);
	    if options_.prefilter == 1
	      yf(:,IdObs) = yf(:,IdObs)+repmat(bayestopt_.mean_varobs', ...
					       horizon+ykmin_,1);
	    end
	    yf(:,IdObs) = yf(:,IdObs)+(gend+[1-ykmin_:horizon]')*trend_coeff';
	    if options_.loglinear == 1
	      yf = yf+repmat(log(ys'),horizon+ykmin_,1);
	      yf = exp(yf);
	    else
	      yf = yf+repmat(ys',horizon+ykmin_,1);
	    end
	    stock_forcst(:,:,irun_forc) = yf;
	    yf1 = forcst2(yyyy,horizon,dr_,1);
	    if options_.prefilter == 1
	      yf1(:,IdObs,:) = yf1(:,IdObs,:)+ ...
		  repmat(bayestopt_.mean_varobs',[horizon+ykmin_,1,1]);
	    end
	    yf1(:,IdObs,:) = yf1(:,IdObs,:)+repmat((gend+[1-ykmin_:horizon]')* ...
		trend_coeff',[1,1,1]);
	    if options_.loglinear == 1
	      yf1 = yf1 + repmat(log(ys'),[horizon+ykmin_,1,1]);
	      yf1 = exp(yf1);
	    else
	      yf1 = yf1 + repmat(ys',[horizon+ykmin_,1,1]);
	    end
	    stock_forcst1(:,:,irun_forc) = yf1;
	    if irun_forc == MAX_nforc
	      sfil_forc = sfil_forc + 1;
	      save([fname_ '_forecast' int2str(sfil_forc)],'stock_forcst','stock_forcst1');
	      irun_forc = 0;
	      stock_forcst = zeros(horizon+ykmin_,nvar,MAX_nforc);
	      stock_forcst1 = zeros(horizon+ykmin_,nvar,MAX_nforc);
	    end
	  end
	  waitbar(b/B,h);    
	end
	if options_.smoother
	  if irun_smoo
	    stock_smooth = stock_smooth(:,:,1:irun_smoo);
	    sfil_smoo = sfil_smoo + 1;
	    instr = [fname_ '_smooth' int2str(sfil_smoo) ' stock_smooth;'];
	    eval(['save ' instr]);
	  end
	  clear stock_smooth;
	  if irun_inno
	    stock_innov = stock_innov(:,:,1:irun_inno);
	    sfil_inno = sfil_inno + 1;
	    instr = [fname_ '_innovation' int2str(sfil_inno) ' stock_innov;'];
	    eval(['save ' instr]);
	  end
	  clear stock_innov;
	end
	if options_.forecast	
	  if irun_forc
	    stock_forcst = stock_forcst(:,:,1:irun_forc);  
	    stock_forcst1 = stock_forcst1(:,:,1:irun_forc);  
	    sfil_forc = sfil_forc + 1;
	    save([fname_ '_forecast' int2str(sfil_forc)],'stock_forcst','stock_forcst1');
	  end
	  clear stock_forcst stock_forcst1
	end
	if options_.filtered_vars	
	  if irun_filt
	    stock_filter = stock_filter(:,:,1:irun_filt);
	    sfil_filt = sfil_filt + 1;
	    instr = [fname_ '_filter' int2str(sfil_filt) ' stock_filter;'];
	    eval(['save ' instr]);
	  end
	  clear stock_filter;
	end
      else % just one _mh file
	if options_.forecast
	  sfil_forc = 0;
	  irun_forc = 0;  			
	end
	if options_.smoother
	  sfil_smoo = 0;
	  sfil_inno = 0;
	  %sfil_erro = 1;
	  irun_smoo = 0;
	  irun_inno = 0;
	  %irun_erro = 0;
	end
	if options_.filtered_vars
	  sfil_filt = 0;
	  irun_filt = 0;  			
	end	  		
	eval(['load ' instr1 int2str(ffil) instr2]);
	NumberOfSimulations = length(logpo2);
	clear post2 logpo2;
	for b = 1:B;
	  if options_.forecast
	    irun_forc = irun_forc+1;
	  end
	  if options_.smoother
	    irun_smoo = irun_smoo+1;
	    irun_inno = irun_inno+1;
	    %irun_erro = irun_erro+1;
	  end
	  if options_.filtered_vars
	    irun_filt = irun_filt+1;  			
	  end	    
	  DEEP  = x2(floor(rand*NumberOfSimulations)+1,:); 
	  deep(subindx) = DEEP(subindx);
	  [atT,innov,obs_err,filtered_state_vector,ys,trend_coeff] = DsgeSmoother(deep',gend,data);
           % removing lagged variables when ykmin_ > 1
           filtered_state_vector = filtered_state_vector(1:endo_nbr,:);
	  if options_.smoother
	    if irun_smoo < MAX_nsmoo
	      stock_smooth(:,:,irun_smoo) = atT(1:endo_nbr,1:gend);
	      if options_.prefilter == 1
		stock_smooth(bayestopt_.mf,:,irun_smoo) = atT(bayestopt_.mf,:)+repmat(bayestopt_.mean_varobs',1,gend);
	      elseif options_.loglinear == 1
		stock_smooth(bayestopt_.mf,:,irun_smoo) = atT(bayestopt_.mf,:)+...
		    repmat(log(ys(bayestopt_.mfys)),1,gend)+...
		     trend_coeff*[1:gend];
	      else
		stock_smooth(bayestopt_.mf,:,irun_smoo) = atT(bayestopt_.mf,:)+...
		    repmat(ys(bayestopt_.mfys),1,gend)+...
		     trend_coeff*[1:gend];
	      end
	    else
	      stock_smooth(:,:,irun_smoo) = atT(1:endo_nbr,1:gend);
	      if options_.prefilter == 1
		stock_smooth(bayestopt_.mf,:,irun_smoo) = atT(bayestopt_.mf,:)+repmat(bayestopt_.mean_varobs',1,gend);
	      elseif options_.loglinear == 1
		stock_smooth(bayestopt_.mf,:,irun_smoo) = atT(bayestopt_.mf,:)+...
		    repmat(log(ys(bayestopt_.mfys)),1,gend)+...
		     trend_coeff*[1:gend];
	      else
		stock_smooth(bayestopt_.mf,:,irun_smoo) = atT(bayestopt_.mf,:)+...
		    repmat(ys(bayestopt_.mfys),1,gend)+...
		     trend_coeff*[1:gend];
	      end
	      sfil_smoo = sfil_smoo + 1;
	      instr = [fname_ '_smooth' int2str(sfil_smoo) ' stock_smooth;'];
	      eval(['save ' instr]);
	      irun_smoo = 0;
	      stock_smooth = ...
		  zeros(endo_nbr,gend,MAX_nsmoo);
	    end	
	    if irun_inno < MAX_ninno
	      stock_innov(:,:,irun_inno) = innov;
	    else
	      stock_innov(:,:,irun_inno) = innov;
	      sfil_inno = sfil_inno + 1;
	      instr = [fname_ '_innovation' int2str(sfil_inno) ' stock_innov;'];
	      eval(['save ' instr]);
	      irun_inno = 0;
	      stock_innov  = zeros(exo_nbr,gend,MAX_ninno);
	    end	
	  end
	  if options_.filtered_vars
	    if irun_filt < MAX_nfilt
	      stock_filter(:,:,irun_filt) = filtered_state_vector;
	    else
	      stock_filter(:,:,irun_filt) = filtered_state_vector;
	      sfil_filt = sfil_filt + 1;
	      instr = [fname_ '_filter' int2str(sfil_filt) ' stock_filter;'];
	      eval(['save ' instr]);
	      irun_filt = 0;
	      stock_filter = ...
		  zeros(endo_nbr,gend+1,MAX_nfilt);
	    end	
	  end    			
	  if options_.forecast	    
	    dr_ = resol(ys_,0);
	    yyyy(:,1:ykmin_) = atT(1:endo_nbr,size(atT,2)-ykmin_+1:size(atT,2));
	    yf = forcst2a(yyyy,dr_,ex_);
	    if options_.prefilter == 1
	      yf(:,IdObs) = yf(:,IdObs)+repmat(bayestopt_.mean_varobs', ...
					       horizon+ykmin_,1);
	    end
	    yf(:,IdObs) = yf(:,IdObs)+(gend+[1-ykmin_:horizon]')*trend_coeff';
	    if options_.loglinear == 1
	      yf = yf+repmat(log(ys'),horizon+ykmin_,1);
	      yf = exp(yf);
	    else
	      yf = yf+repmat(ys',horizon+ykmin_,1);
	    end
	    stock_forcst(:,:,irun_forc) = yf;
	    yf1 = forcst2(yyyy,horizon,dr_,1);
	    if options_.prefilter == 1
	      yf1(:,IdObs,:) = yf1(:,IdObs,:)+ ...
		  repmat(bayestopt_.mean_varobs',[horizon+ykmin_,1,1]);
	    end
	    yf1(:,IdObs,:) = yf1(:,IdObs,:)+repmat((gend+[1-ykmin_:horizon]')* ...
						   trend_coeff',[1,1,1]);
	    if options_.loglinear == 1
	      yf1 = yf1 + repmat(log(ys'),[horizon+ykmin_,1,1]);
	      yf1 = exp(yf1);
	    else
	      yf1 = yf1 + repmat(ys',[horizon+ykmin_,1,1]);
	    end
	    stock_forcst1(:,:,irun_forc) = yf1;
	    if irun_forc == MAX_nforc
	      sfil_forc = sfil_forc + 1;
	      save([fname_ '_forecast' int2str(sfil_forc)],'stock_forcst','stock_forcst1');
	      irun_forc = 0;
	      stock_forcst = zeros(horizon+ykmin_,nvar,MAX_nforc);
	      stock_forcst1 = zeros(horizon+ykmin_,nvar,MAX_nforc);
	    end
	  end   
	  waitbar(b/B,h);    
	end
	if options_.smoother
	  if irun_smoo
	    stock_smooth = stock_smooth(:,:,1:irun_smoo);
	    sfil_smoo = sfil_smoo + 1;
	    instr = [fname_ '_smooth' int2str(sfil_smoo) ' stock_smooth;'];
	    eval(['save ' instr]);
	  end
	  clear stock_smooth;
	  if irun_inno
	    stock_innov = stock_innov(:,:,1:irun_inno);
	    sfil_inno = sfil_inno + 1;
	    instr = [fname_ '_innovation' int2str(sfil_inno) ' stock_innov;'];
	    eval(['save ' instr]);
	  end
	  clear stock_innov;
	end
	if options_.forecast    
	  if irun_forc
	    stock_forcst = stock_forcst(:,:,1:irun_forc);  
	    stock_forcst1 = stock_forcst1(:,:,1:irun_forc);  
	    sfil_forc = sfil_forc + 1;
	    save([fname_ '_forecast' int2str(sfil_forc)],'stock_forcst','stock_forcst1');
	  end
	  clear stock_forcst stock_forcst1
	end
	if options_.filtered_vars   
	  if irun_filt
	    stock_filter = stock_filter(:,:,1:irun_filt);
	    sfil_filt = sfil_filt + 1;
	    instr = [fname_ '_filter' int2str(sfil_filt) ' stock_filter;'];
	    eval(['save ' instr]);
	  end
	  clear stock_filter;
	end   
      end
    end
    close(h);
  end
  %%
  %% Only a subset of variables may be treated    
  %%
  varlist = options_.varlist;
  if isempty(varlist)
    varlist = lgy_;
    nvar    = size(lgy_,1);
    SelecVariables = transpose(1:nvar);
  else
    nvar = size(varlist,1);
    SelecVariables = [];
    for i=1:nvar
      if ~isempty(strmatch(deblank(varlist(i,:)),lgy_,'exact'))
        SelecVariables = [ SelecVariables ; strmatch(deblank(varlist(i,:)),lgy_,'exact') ];
      end
    end
    IdObs    = zeros(nvobs,1);
    for j=1:nvobs
      for i=1:nvar
        iobs = strmatch(options_.varobs(j,:),varlist,'exact');
      end
      if ~isempty(iobs)
        IdObs(j,1) = iobs;
      end    
    end 
  end
  if TeX
    varlist_TeX = [];
    for i=1:nvar
      varlist_TeX = strvcat(varlist_TeX,lgy_TeX_(SelecVariables(i),:));
    end
  end
  %%                                    %%
  %% Now I treat the forecasts (plots)  %%   
  %%                                    %%
  if options_.forecast
    tmp = zeros(B,1);
    tmp_big = zeros(B,1);
    fprintf('MH: Out of sample forecasts...\n');
    MeanForecast = zeros(options_.forecast,nvar);
    MedianForecast = zeros(options_.forecast,nvar);
    StdForecast = zeros(options_.forecast,nvar);
    HPD   = zeros(options_.forecast,nvar,2);
    StdForecast_total = zeros(options_.forecast,nvar);
    HPD_total   = zeros(options_.forecast,nvar,2);
    for step = 1:options_.forecast % ... Suffering is one very long moment.
      truestep = step+ykmin_;
      for i = 1:nvar;
	StartLine = 0;
	StartLine1 = 0;
	for file = 1:sfil_forc;
	  load([fname_ '_forecast' int2str(file)]);
	  MeanForecast(step,i) = MeanForecast(step,i)+sum(stock_forcst(truestep,SelecVariables(i),:),3);
	  DeProfundis = size(stock_forcst,3); 
	  DeProfundis1 = size(stock_forcst1,3); 
	  tmp(StartLine+1:StartLine+DeProfundis) = squeeze(stock_forcst(truestep,SelecVariables(i),:)); 
	  tmp_big(StartLine1+1:StartLine1+DeProfundis1) = squeeze(stock_forcst1(truestep,SelecVariables(i),:)); 
	  StartLine = StartLine+DeProfundis;
	  StartLine1 = StartLine1+DeProfundis1;
	end
	tmp = sort(tmp);
	tmp_big = sort(tmp_big);
	MedianForecast(step,i) = tmp(round(B*0.5));
	StdForecast(step,i) = std(tmp);
	StdForecast_total(step,i) = std(tmp_big);
	t = floor(options_.mh_conf_sig*B);
	a = 1; 
	b = t;
	tmp2 = [1;t;tmp(t)-tmp(1)];
	while b <= B
	  tmp1 = [a;b;tmp(b)-tmp(a)];
	  a = a + 1;
	  b = b + 1;
	  if tmp1(3) < tmp2(3)
	    tmp2 = tmp1;     
	  end    
	end
	HPD(step,i,1) = tmp(tmp2(1));
	HPD(step,i,2) = tmp(tmp2(2));
	t = floor(options_.mh_conf_sig*B);
	a = 1; 
	b = t;
	tmp2_big = [1;t;tmp_big(t)-tmp_big(1)];
	while b <= B
	  tmp1_big = [a;b;tmp_big(b)-tmp_big(a)];
	  a = a + 1;
	  b = b + 1;
	  if tmp1_big(3) < tmp2_big(3)
	    tmp2_big = tmp1_big;     
	  end    
	end
	HPD_total(step,i,1) = tmp_big(tmp2_big(1));
	HPD_total(step,i,2) = tmp_big(tmp2_big(2));
      end
      disp(['    Period = ' int2str(step)]);
    end
    MeanForecast = MeanForecast/B;
    [nbplt,nr,nc,lr,lc,nstar] = pltorg(nvar);
    if TeX
      fidTeX = fopen([fname_ '_BayesianForecasts.TeX'],'w');
      fprintf(fidTeX,'%% TeX eps-loader file generated by metropolis.m (Dynare).\n');
      fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
      fprintf(fidTeX,' \n');
      NAMES = [];
      TEXNAMES = [];
    end    
    for plt = 1:nbplt
      if TeX
        NAMES = [];
        TEXNAMES = [];
      end
      hfig = figure('Name','Out of sample forecasts');
      for i = 1:nstar
        k = (plt-1)*nstar+i;
        if k > nvar
          break
        end
        subplot(nr,nc,i)
        hold on
        if any(k==IdObs)
          idx = find(k==IdObs);
          if options_.loglinear == 1
            plot(1:10+options_.forecast,...
                 [exp(data(idx,size(data,2)-10+1:end))';...
                  MeanForecast(:,k)],'-b','linewidth',2)
          else
            plot(1:10+options_.forecast,[data(idx,size(data,2)-10+1:end)';MeanForecast(:,k)],'-b','linewidth',2)
          end
          offsetx = 10;
        else
          plot(1:options_.forecast,MeanForecast(:,k),'-b', ...
               'linewidth',2)
          offsetx = 0;
        end   
        plot(offsetx+[1:options_.forecast],HPD(:,k,1),'--g', ...
             'linewidth',1.5)
        plot(offsetx+[1:options_.forecast],HPD(:,k,2),'--g', ...
             'linewidth',1.5)
        plot(offsetx+[1:options_.forecast],HPD_total(:,k,1),'--r', ...
             'linewidth',1.5)
        plot(offsetx+[1:options_.forecast],HPD_total(:,k,2),'--r','linewidth',1.5)
        set(gca,'XTick',offsetx+[1 10 20 30 40 50 60 70 80 90]);
        set(gca,'XTickLabel',{'1';'10';'20';'30';'40';'50';'60';'70';'80';'90'});
        %   xlim([1 options_.forecast+10]);
        if any(k==IdObs)
          plot([11 11],ylim,'-c')
        end
        box on
        title(deblank(varlist(k,:)),'Interpreter','none')
        hold off
        eval(['oo_.Forecast.Mean.' deblank(varlist(k,:)) ' = MeanForecast(:,k)'';']);
        eval(['oo_.Forecast.Median.' deblank(varlist(k,:)) ' = MedianForecast(:,k)'';']);
        eval(['oo_.Forecast.Std.' deblank(varlist(k,:)) ' = StdForecast(:,k)'';']);
        eval(['oo_.Forecast.HPDinf.' deblank(varlist(k,:)) ' = squeeze(HPD(:,k,1))'';']);
        eval(['oo_.Forecast.HPDsup.' deblank(varlist(k,:)) ' = squeeze(HPD(:,k,2))'';']);
        eval(['oo_.Forecast.HPDTotalinf.' deblank(varlist(k,:)) ' = squeeze(HPD_total(:,k,1))'';']);
        eval(['oo_.Forecast.HPDTotalsup.' deblank(varlist(k,:)) ' = squeeze(HPD_total(:,k,2))'';']);
        if TeX
          NAMES = strvcat(NAMES,deblank(varlist(k,:)));
          TEXNAMES = strvcat(TEXNAMES,['$ ' deblank(varlist_TeX(k,:)) ' $']);
        end
      end
      eval(['print -depsc2 ' fname_ '_Forecasts' int2str(plt)]);
      eval(['print -dpdf ' fname_ '_Forecasts' int2str(plt)]);
      saveas(hfig,[fname_ '_Forecasts' int2str(plt) '.fig']);
      if options_.nograph, close(hfig), end
      if TeX
        fprintf(fidTeX,'\\begin{figure}[H]\n');
        for jj = 1:nstar
          fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TEXNAMES(jj,:)));
        end    
        fprintf(fidTeX,'\\centering \n');
        fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_Forecasts%s}\n',fname_,int2str(plt));
        fprintf(fidTeX,'\\caption{DSGE posterior mean forecats with HPD intervals.}');
        fprintf(fidTeX,'\\label{Fig:Forecasts:%s}\n',int2str(plt));
        fprintf(fidTeX,'\\end{figure}\n');
        fprintf(fidTeX,' \n');
      end
    end
    fprintf('MH: Out of sample forecasts, done!\n')
    disp(' ')
  end
  %%
  %% Smooth variables and Filtered variables (all endogenous variables are considered here)        
  %%
  if options_.smoother
    [MeanSmooth,MedianSmooth,StdSmooth,DistribSmooth,HPDSmooth] = GetPosteriorStatistics(gend,B,'SmoothedVariables');    
    [MeanInnov,MedianInnov,StdInnov,DistribInnov,HPDInnov] = GetPosteriorStatistics(gend,B,'SmoothedShocks');
    if nvn
      [MeanError,MedianError,StdError,DistribError,HPDError] = GetPosteriorStatistics(gend,B,'SmoothedObservationErrors');
    end
    %%
    %% Now I plot the smooth -structural- shocks
    %%
    MakeSmoothVariablesPlots('SmoothedShocks',DistribInnov,MeanInnov,gend)
    %%
    %% Smoothed variables (observed and unobserved)
    %%
    MakeSmoothVariablesPlots('SmoothedVariables',DistribSmooth,MeanSmooth,gend)
    %%
    %% Smoothed observation error
    %%
    if nvn
      MakeSmoothVariablesPlots('SmoothedObservationErrors',DistribError,MeanError,gend)
      %%
      %% Historical and smoothed variabes
      %%
      MakeSmoothVariablesPlots('Historical&SmoothedObservableVariables',...
                   MeanSmooth,rawdata(options_.first_obs+(0:gend-1),:),gend)
    end  
    %%
    %%
    %%
  end % options_.smoother
  if options_.filtered_vars % Filtered variables.
    [MeanFilter,MedianFilter,StdFilter,DistribFilter,HPDFilter] = GetPosteriorStatistics(gend,B,'FilteredVariables');
    MakeSmoothVariablesPlots('FilteredVariables',DistribFilter,MeanFilter,gend)
  end          
  %%
  %%    Posterior IRFs. Instead of displaying the IRFs associated to the posterior mean
  %%    of the structural parameters (by calling stoch_simul after estimation), 
  %%    metropolis.m will display the posterior mean of the IRFs and the deciles of 
  %%    the IRFs' posterior distribution. All the results are saved in the global 
  %%    structure oo_ (posterior medians, posterior standard deviations and posterior HPD   
  %%    intervals are also computed and saved).
  %%
  if options_.bayesian_irf
    deep = MU;
    subindx = subset();
    nirfs = options_.irf;
    if ~isempty(dsge_prior_weight)
      files = eval(['dir(''' fname_ '_irf_dsgevar*.mat'');']);     
      if size(files,1)                                        
	delete([fname_ '_irf_dsgevar*.mat']);                    
	disp(['MH: Old ' fname_ '_irf_dsgevar files deleted! ']) 
      end
      files = eval(['dir(''' fname_ '_irf_dsge*.mat'');']);     
      if size(files,1)                                        
	delete([fname_ '_irf_dsge*.mat']);                    
	disp(['MH: Old ' fname_ '_irf_dsge files deleted! ']) 
      end
    else
      files = eval(['dir(''' fname_ '_irf_dsge*.mat'');']);     
      if size(files,1)                                        
	delete([fname_ '_irf_dsge*.mat']);                    
	disp(['MH: Old ' fname_ '_irf_dsge files deleted! ']) 
      end      
    end    
    if B <= MAX_nirfs_dsge
      stock_irf_dsge = zeros(nirfs,size(lgy_,1),exo_nbr,B);
    elseif nvn & B > MAX_nirfs_dsge
      stock_irf_dsge = zeros(nirfs,size(lgy_,1),exo_nbr,MAX_nirfs_dsge);
    end
    if ~isempty(dsge_prior_weight)
      if B <= MAX_nirfs_dsgevar
	stock_irf_dsgevar = zeros(nirfs,nvobs,exo_nbr,B);
      else
	stock_irf_dsgevar = zeros(nirfs,nvobs,exo_nbr,MAX_nirfs_dsgevar);
      end
      [mYY,mXY,mYX,mXX,Ydata,Xdata] = ...
	  VarSampleMoments(options_.first_obs,options_.first_obs+options_.nobs-1,options_.varlag,-1);
      NumberOfLags = options_.varlag;
      NumberOfLagsTimesNvobs = NumberOfLags*nvobs;
      COMP_draw = diag(ones(nvobs*(NumberOfLags-1),1),-nvobs);
    end
    h = waitbar(0,'Bayesian IRFs...');
    if nfile-ffil+1>1
      sfil_irf_dsge = 0;
      irun_irf_dsge = 0;
      sfil_irf_dsgevar = 0;
      irun_irf_dsgevar = 0;
      for b = 1:B;
	irun_irf_dsge = irun_irf_dsge+1;
	tmp_dsge = zeros(nirfs,size(lgy_,1),exo_nbr);
	if ~isempty(dsge_prior_weight)
	  irun_irf_dsgevar = irun_irf_dsgevar+1;
	  tmp_dsgevar = zeros(nirfs,nvobs*exo_nbr);
	end
	choose_an_mh_file = rand;
	mh_file_number = FLN(find(choose_an_mh_file>=FLN(:,3)),1);
	if isempty(mh_file_number)
	  mh_file_number = ffil;
	else    
	  mh_file_number = mh_file_number(1);
	end
	eval(['load ' instr1 int2str(mh_file_number) instr2]);
	clear post2 logpo2;
	DEEP  = x2(floor(rand*FLN(find(mh_file_number == FLN(:,1)),2))+1,:);
	deep(subindx) = DEEP(subindx);
	% dsge
	set_parameters(deep);
	dr_ = resol(ys_,0);
	SS(lgx_orig_ord_,lgx_orig_ord_)=Sigma_e_+1e-14* ...
	    eye(exo_nbr);
	cs = transpose(chol(SS));
	tit(lgx_orig_ord_,:) = lgx_;
	for i = 1:exo_nbr
	  if SS(i,i) > 1e-13
	    y=irf(dr_,cs(lgx_orig_ord_,i),nirfs,options_.drop,options_.replic, options_.order);
	    if options_.relative_irf
	      y = 100*y/cs(i,i); 
	    end
	    for j = 1:size(lgy_,1)
	      if max(y(j,:)) - min(y(j,:)) > 1e-10 
		tmp_dsge(:,j,i) = transpose(y(j,:));
	      end
	    end
	  end
	end
	if irun_irf_dsge < MAX_nirfs_dsge
	  stock_irf_dsge(:,:,:,irun_irf_dsge) = tmp_dsge;
	else
	  stock_irf_dsge(:,:,:,irun_irf_dsge) = tmp_dsge;
	  sfil_irf_dsge = sfil_irf_dsge + 1;
	  instr = [fname_ '_irf_dsge' int2str(sfil_irf_dsge) ' stock_irf_dsge;'];
	  eval(['save ' instr]);
	  irun_irf_dsge = 0;
	  stock_irf_dsge = zeros(nirfs,size(lgy_,1),exo_nbr,MAX_nirfs_dsge);
	end
	% bvar-dsge 
	if ~isempty(dsge_prior_weight)
	  [fval,cost_flag,ys,trend_coeff,info,PHI,SIGMAu,iXX] = DsgeVarLikelihood(deep',gend);
	  DSGE_PRIOR_WEIGHT = floor(gend*(1+dsge_prior_weight));
	  tmp1 = SIGMAu*gend*(dsge_prior_weight+1);
          val  = 1;
          tmp1 = chol(inv(tmp1))'; 
	  while val;
	    % draw from the marginal posterior of sig
	    tmp2 = tmp1*randn(nvobs,DSGE_PRIOR_WEIGHT-NumberOfLagsTimesNvobs);
	    SIGMAu_draw = inv(tmp2*tmp2');
            % draw from the conditional posterior of PHI
	    VARvecPHI = kron(SIGMAu_draw,iXX);
	    PHI_draw  = PHI(:) + chol(VARvecPHI)'*randn(nvobs*NumberOfLagsTimesNvobs,1);
	    COMP_draw(1:nvobs,:) = reshape(PHI_draw,NumberOfLagsTimesNvobs,nvobs)';
	    % Check for stationarity
	    tests = find(abs(eig(COMP_draw))>0.9999999999);
	    if isempty(tests)
	      val=0;
	    end
	  end
	  % Get rotation
	  if dsge_prior_weight > 0
	    Atheta(dr_.order_var,lgx_orig_ord_) = dr_.ghu*sqrt(Sigma_e_);
	    A0 = Atheta(bayestopt_.mfys,:);
	    [OMEGAstar,SIGMAtr] = qr2(A0');
	  end
	  SIGMAu_chol = chol(SIGMAu_draw)';
	  SIGMAtrOMEGA = SIGMAu_chol*OMEGAstar';
	  PHIpower = eye(NumberOfLagsTimesNvobs);
	  irfs = zeros (nirfs,nvobs*exo_nbr);
	  tmp3 = PHIpower(1:nvobs,1:nvobs)*SIGMAtrOMEGA;
	  irfs(1,:) = tmp3(:)';
	  for t = 2:nirfs
	    PHIpower = COMP_draw*PHIpower;
	    tmp3 = PHIpower(1:nvobs,1:nvobs)*SIGMAtrOMEGA;
	    irfs(t,:)  = tmp3(:)';
	  end            
	  for j = 1:(nvobs*exo_nbr)
	    if max(irfs(:,j)) - min(irfs(:,j)) > 1e-10 
	      tmp_dsgevar(:,j) = (irfs(:,j));
	    end
	  end
	  if irun_irf_dsgevar < MAX_nirfs_dsgevar
	    stock_irf_dsgevar(:,:,:,irun_irf_dsgevar) = reshape(tmp_dsgevar,nirfs,nvobs,exo_nbr);
	  else
	    stock_irf_dsgevar(:,:,:,irun_irf_dsgevar) = reshape(tmp_dsgevar,nirfs,nvobs,exo_nbr);
	    sfil_irf_dsgevar = sfil_irf_dsgevar + 1;
	    instr = [fname_ '_irf_dsgevar' int2str(sfil_irf_dsgevar) ' stock_irf_dsgevar;'];,
	    eval(['save ' instr]);
	    irun_irf_dsgevar  = 0;
	    stock_irf_dsgevar = zeros(nirfs,nvobs,exo_nbr,MAX_nirfs_dsgevar);
	  end
	end
	waitbar(b/B,h);    
      end
      clear tmp;
      if irun_irf_dsge
	stock_irf_dsge = stock_irf_dsge(:,:,:,1:irun_irf_dsge);
	sfil_irf_dsge = sfil_irf_dsge + 1;
	instr = [fname_ '_irf_dsge' int2str(sfil_irf_dsge) ' stock_irf_dsge;'];
	eval(['save ' instr]);
      end
      clear stock_irf_dsge;
      if ~isempty(dsge_prior_weight)
	if irun_irf_dsgevar
	  stock_irf_dsgevar = stock_irf_dsgevar(:,:,:,1:irun_irf_dsgevar);
	  sfil_irf_dsgevar = sfil_irf_dsgevar + 1;
	  instr = [fname_ '_irf_dsgevar' int2str(sfil_irf_dsgevar) ' stock_irf_dsgevar;'];
	  eval(['save ' instr]);
	end
	clear stock_irf_dsgevar;
      end
    else
      sfil_irf_dsge = 0;
      irun_irf_dsge = 0;
      if ~isempty(dsge_prior_weight)
	sfil_irf_dsgevar = 0;
	irun_irf_dsgevar = 0;
      end
      eval(['load ' instr1 int2str(ffil) instr2]);
      NumberOfSimulations = length(logpo2);
      clear post2 logpo2;
      for b = 1:B;
	irun_irf_dsge = irun_irf_dsge+1;
	tmp_dsge = zeros(nirfs,size(lgy_,1),exo_nbr);
	if ~isempty(dsge_prior_weight)
	  irun_irf_dsgevar = irun_irf_dsgevar+1;
	  tmp_dsgevar = zeros(nirfs,nvobs*exo_nbr);	  
	end
	DEEP = x2(floor(rand*NumberOfSimulations)+1,:);
	deep(subindx) = DEEP(subindx);
	% dsge
	set_parameters(deep);
	dr_ = resol(ys_,0);
	SS(lgx_orig_ord_,lgx_orig_ord_) = Sigma_e_+1e-14*eye(exo_nbr);
	SS = transpose(chol(SS));
	tit(lgx_orig_ord_,:) = lgx_;
	for i = 1:exo_nbr
	  if SS(i,i) > 1e-13
	    y=irf(dr_,SS(lgx_orig_ord_,i), options_.irf, options_.drop, ...
		  options_.replic, options_.order);
	    if options_.relative_irf
	      y = 100*y/cs(i,i); 
	    end
	    for j = 1:size(lgy_,1)
	      if max(y(j,:)) - min(y(j,:)) > 1e-10 
		tmp_dsge(:,j,i) = transpose(y(j,:));
	      end
	    end
	  end
	end
	if irun_irf_dsge < MAX_nirfs_dsge
	  stock_irf_dsge(:,:,:,irun_irf_dsge) = tmp_dsge;
	else
	  stock_irf_dsge(:,:,:,irun_irf_dsge) = tmp_dsge;
	  sfil_irf_dsge = sfil_irf_dsge + 1;
	  instr = [fname_ '_irf_dsge' int2str(sfil_irf_dsge) ' stock_irf_dsge;'];
	  eval(['save ' instr]);
	  irun_irf_dsge = 0;
	  stock_irf_dsge = zeros(nirfs,size(lgy_,1),exo_nbr,MAX_nirfs_dsge);
	end
	% bvar-dsge
	if ~isempty(dsge_prior_weight)
	  [fval,cost_flag,ys,trend_coeff,info,PHI,SIGMAu,iXX] = DsgeVarLikelihood(deep',gend);
	  DSGE_PRIOR_WEIGHT = floor(gend*(1+dsge_prior_weight));
          tmp1 = SIGMAu*gend*(1+dsge_prior_weight);
	  tmp1 = chol(inv(tmp1))';
          val = 1;
	  while val;
	    % draw from the marginal posterior of sig
	    tmp2 = tmp1*randn(nvobs,DSGE_PRIOR_WEIGHT-NumberOfLagsTimesNvobs);
	    SIGMAu_draw = inv(tmp2*tmp2');
	    % draw from the conditional posterior of PHI
	    VARvecPHI = kron(SIGMAu_draw,iXX);
	    PHI_draw  = PHI(:) + chol(VARvecPHI)'*randn(nvobs*NumberOfLagsTimesNvobs,1);
	    COMP_draw(1:nvobs,:) = reshape(PHI_draw,NumberOfLagsTimesNvobs,nvobs)';
	    % Check for stationarity
	    tests = find(abs(eig(COMP_draw))>0.9999999999);
	    if isempty(tests)
	      val=0;
	    end
	  end
	  % Get rotation
	  if dsge_prior_weight > 0
	    Atheta(dr_.order_var,lgx_orig_ord_) = dr_.ghu*sqrt(Sigma_e_);%%%
	    A0 = Atheta(bayestopt_.mfys,:);
	    [OMEGAstar,SIGMAtr] = qr2(A0');
	  end
	  SIGMAu_chol = chol(SIGMAu_draw)';
	  SIGMAtrOMEGA = SIGMAu_chol*OMEGAstar';
	  PHIpower = eye(NumberOfLagsTimesNvobs);
	  irfs = zeros (nirfs,nvobs*exo_nbr);
	  tmp3 = PHIpower(1:nvobs,1:nvobs)*SIGMAtrOMEGA;
	  irfs(1,:) = tmp3(:)';
	  for t = 2:nirfs
	    PHIpower = COMP_draw*PHIpower;
	    tmp3 = PHIpower(1:nvobs,1:nvobs)*SIGMAtrOMEGA;
	    irfs(t,:)  = tmp3(:)';
	  end
	  for j = 1:(nvobs*exo_nbr)
	    if max(irfs(:,j)) - min(irfs(:,j)) > 1e-10 
	      tmp_dsgevar(:,j) = (irfs(:,j));
	    end	
	  end	
	  if irun_irf_dsgevar < MAX_nirfs_dsgevar
	    stock_irf_dsgevar(:,:,:,irun_irf_dsgevar) = reshape(tmp_dsgevar,nirfs,nvobs,exo_nbr);
	  else
	    stock_irf_dsgevar(:,:,:,irun_irf_dsgevar) = reshape(tmp_dsgevar,nirfs,nvobs,exo_nbr);
	    sfil_irf_dsgevar = sfil_irf_dsgevar + 1;
	    instr = [fname_ '_irf_dsgevar' int2str(sfil_irf_dsgevar) ' stock_irf_dsgevar;'];,
	    eval(['save ' instr]);
	    irun_irf_dsgevar = 0;
	    stock_irf_dsgevar = zeros(nirfs,nvobs,exo_nbr,MAX_nirfs_dsgevar);
	  end
	end
	waitbar(b/B,h);
      end
      if irun_irf_dsge
	stock_irf_dsge = stock_irf_dsge(:,:,:,1:irun_irf_dsge);
	sfil_irf_dsge = sfil_irf_dsge + 1;
	instr = [fname_ '_irf_dsge' int2str(sfil_irf_dsge) ' stock_irf_dsge;'];
	eval(['save ' instr]);
      end
      clear stock_irf_dsge;
      if ~isempty(dsge_prior_weight)
	if irun_irf_dsgevar
	  stock_irf_dsgevar = stock_irf_dsgevar(:,:,:,1:irun_irf_dsgevar);
	  sfil_irf_dsgevar = sfil_irf_dsgevar + 1;
	  instr = [fname_ '_irf_dsgevar' int2str(sfil_irf_dsgevar) ' stock_irf_dsgevar;'];
	  eval(['save ' instr]);
	end
	clear stock_irf_dsgevar;    
      end
    end
    close(h)
    %%
    %%  Now i compute some statistics (mean, median, std, deciles, HPD intervals)
    %%
    %%  DSGE:
    tmp = zeros(B,1);
    MeanIRF_dsge = zeros(nirfs,nvar,exo_nbr);
    MedianIRF_dsge = zeros(nirfs,nvar,exo_nbr);
    StdIRF_dsge = zeros(nirfs,nvar,exo_nbr);
    DistribIRF_dsge = zeros(nirfs,nvar,exo_nbr,9);
    HPDIRF_dsge = zeros(nirfs,nvar,exo_nbr,2);
    if ~isempty(dsge_prior_weight)
      SelecVariables = bayestopt_.mfys;
      nvar = length(SelecVariables);
    end
    fprintf('MH: Posterior IRFs (dsge)...\n')
    for i = 1:exo_nbr
      for j = 1:nvar
	for k = 1:nirfs
	  StartLine = 0;
	  for file = 1:sfil_irf_dsge;
	    instr = [fname_ '_irf_dsge' int2str(file)];
	    eval(['load ' instr]);
	    MeanIRF_dsge(k,j,i) = MeanIRF_dsge(k,j,i)+sum(stock_irf_dsge(k,SelecVariables(j),i,:),4);
	    DeProfundis = size(stock_irf_dsge,4); 
	    tmp(StartLine+1:StartLine+DeProfundis) = squeeze(stock_irf_dsge(k,SelecVariables(j),i,:)); 
	    StartLine = StartLine+DeProfundis;
	  end
	  tmp = sort(tmp);
	  MedianIRF_dsge(k,j,i) = tmp(round(B*0.5));
	  StdIRF_dsge(k,j,i) = std(tmp);
	  DistribIRF_dsge(k,j,i,:) = reshape(tmp(deciles),1,1,1,9);
	  tt = floor(options_.mh_conf_sig*B);
	  a = 1; 
	  b = tt;
	  tmp2 = [1;tt;tmp(tt)-tmp(1)];
	  while b <= B
	    tmp1 = [a;b;tmp(b)-tmp(a)];
	    a = a + 1;
	    b = b + 1;
	    if tmp1(3,1) < tmp2(3,1)
	      tmp2 = tmp1;     
	    end
	  end
	  HPDIRF_dsge(k,j,i,1) = tmp(tmp2(1,1));
	  HPDIRF_dsge(k,j,i,2) = tmp(tmp2(2,1));
	end
	disp(['    Variable: ' deblank(lgy_(SelecVariables(j),:)) ', orthogonalized shock to ' deblank(tit(i,:))])  
      end
    end
    clear stock_irf_dsge;
    MeanIRF_dsge = MeanIRF_dsge/B;
    for i = 1:exo_nbr
      for j = 1:nvar
	eval(['oo_.PosteriorIRF.Dsge.Mean.' ...
	      deblank(lgy_(SelecVariables(j),:)) '_' deblank(tit(i,:)) ' = MeanIRF_dsge(:,j,i);']);
	eval(['oo_.PosteriorIRF.Dsge.Median.' ...
	      deblank(lgy_(SelecVariables(j),:)) '_' deblank(tit(i,:)) ' = MedianIRF_dsge(:,j,i);']);
	eval(['oo_.PosteriorIRF.Dsge.Std.' ...
	      deblank(lgy_(SelecVariables(j),:)) '_' deblank(tit(i,:)) ' = StdIRF_dsge(:,j,i);']);
	eval(['oo_.PosteriorIRF.Dsge.Distribution.' ...
	      deblank(lgy_(SelecVariables(j),:)) '_' deblank(tit(i,:)) ' = squeeze(DistribIRF_dsge(:,j,i,:));']);
	eval(['oo_.PosteriorIRF.Dsge.HPDinf.' ...
	      deblank(lgy_(SelecVariables(j),:)) '_' deblank(tit(i,:)) ' = squeeze(HPDIRF_dsge(:,j,i,1));']);
	eval(['oo_.PosteriorIRF.Dsge.HPDsup.' ...
	      deblank(lgy_(SelecVariables(j),:)) '_' deblank(tit(i,:)) ' = squeeze(HPDIRF_dsge(:,j,i,2));']);
      end
    end
    if ~isempty(dsge_prior_weight)% BVAR-DSGE:
      tmp = zeros(B,1);
      MeanIRF_dsgevar = zeros(nirfs,nvar,exo_nbr);
      MedianIRF_dsgevar = zeros(nirfs,nvar,exo_nbr);
      StdIRF_dsgevar = zeros(nirfs,nvar,exo_nbr);
      DistribIRF_dsgevar = zeros(nirfs,nvar,exo_nbr,9);
      HPDIRF_dsgevar = zeros(nirfs,nvar,exo_nbr,2);
      disp('')
      fprintf('MH: Posterior IRFs (bvar-dsge)...\n')
      for i = 1:exo_nbr
	for j = 1:nvobs
	  for k = 1:nirfs
	    StartLine = 0;
	    for file = 1:sfil_irf_dsgevar
	      instr = [fname_ '_irf_dsgevar' int2str(file)];
	      eval(['load ' instr]);
	      MeanIRF_dsgevar(k,j,i) = MeanIRF_dsgevar(k,j,i)+sum(stock_irf_dsgevar(k,j,i,:),4);
	      DeProfundis = size(stock_irf_dsgevar,4); 
	      tmp(StartLine+1:StartLine+DeProfundis) = squeeze(stock_irf_dsgevar(k,j,i,:)); 
	      StartLine = StartLine+DeProfundis;
	    end
	    tmp = sort(tmp);
	    MedianIRF_dsgevar(k,j,i) = tmp(round(B*0.5));
	    StdIRF_dsgevar(k,j,i) = std(tmp);
	    DistribIRF_dsgevar(k,j,i,:) = reshape(tmp(deciles),1,1,1,9);
	    tt = floor(options_.mh_conf_sig*B);
	    a = 1; 
	    b = tt;
	    tmp2 = [1;tt;tmp(tt)-tmp(1)];
	    while b <= B
	      tmp1 = [a;b;tmp(b)-tmp(a)];
	      a = a + 1;
	      b = b + 1;
	      if tmp1(3,1) < tmp2(3,1)
		tmp2 = tmp1;     
	      end    
	    end
	    HPDIRF_dsgevar(k,j,i,1) = tmp(tmp2(1,1));
	    HPDIRF_dsgevar(k,j,i,2) = tmp(tmp2(2,1));
	  end
	  disp(['    Variable: ' deblank(options_.varobs(j,:)) ', orthogonalized shock to ' deblank(tit(i,:))])  
	end   
      end
      clear stock_irf_dsgevar;
      MeanIRF_dsgevar = MeanIRF_dsgevar/B;
      for i = 1:exo_nbr
	for j = 1:nvobs
	  eval(['oo_.PosteriorIRF.BvarDsge.Mean.' ...
		deblank(options_.varobs(j,:)) '_' deblank(tit(i,:)) ' = MeanIRF_dsgevar(:,j,i);']);
	  eval(['oo_.PosteriorIRF.BvarDsge.Median.' ...
		deblank(options_.varobs(j,:)) '_' deblank(tit(i,:)) ' = MedianIRF_dsgevar(:,j,i);']);
	  eval(['oo_.PosteriorIRF.BvarDsge.Std.' ...
		deblank(options_.varobs(j,:)) '_' deblank(tit(i,:)) ' = StdIRF_dsgevar(:,j,i);']);
	  eval(['oo_.PosteriorIRF.BvarDsge.Distribution.' ...
		deblank(options_.varobs(j,:)) '_' deblank(tit(i,:)) ' = squeeze(DistribIRF_dsgevar(:,j,i,:));']);
	  eval(['oo_.PosteriorIRF.BvarDsge.HPDinf.' ...
		deblank(options_.varobs(j,:)) '_' deblank(tit(i,:)) ' = squeeze(HPDIRF_dsgevar(:,j,i,1));']);
	  eval(['oo_.PosteriorIRF.BvarDsge.HPDsup.' ...
		deblank(options_.varobs(j,:)) '_' deblank(tit(i,:)) ' = squeeze(HPDIRF_dsgevar(:,j,i,2));']);
	end
      end
    end
    %%
    %%  Finally I build the plots.
    %%
    if TeX
      fidTeX = fopen([fname_ '_BayesianIRF.TeX'],'w');
      fprintf(fidTeX,'%% TeX eps-loader file generated by metropolis.m (Dynare).\n');
      fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
      fprintf(fidTeX,' \n');
    end
    tit(lgx_orig_ord_,:) = lgx_;
    if TeX; titTeX(lgx_orig_ord_,:) = lgx_TeX_; end;
    for i=1:exo_nbr
      number_of_plots_to_draw = 0;
      index = [];
      if ~~isempty(dsge_prior_weight) 
	for j=1:nvar
	  if MeanIRF_dsge(1,j,i)
	    number_of_plots_to_draw = number_of_plots_to_draw + 1;
	    index = cat(1,index,j);
	  end
	end
      else% BVAR-DSGE
	number_of_plots_to_draw = nvar;
	index = (1:nvar)';
      end
      [nbplt,nr,nc,lr,lc,nstar] = pltorg(number_of_plots_to_draw);  
      if nbplt == 1
	if options_.relative_irf
	  hh = figure('Name',['Relative response to orthogonalized shock to ' tit(i,:)]);
	else
	  hh = figure('Name',['Orthogonalized shock to ' tit(i,:)]);
	end
	NAMES = [];
	if TeX; TEXNAMES = []; end;
	for j=1:number_of_plots_to_draw
	  set(0,'CurrentFigure',hh)
	  subplot(nr,nc,j);
	  plot([1 nirfs],[0 0],'-r','linewidth',0.5);% zero line.
	  hold on
	  for k = 1:9
	    plot(1:nirfs,DistribIRF_dsge(:,index(j),i,k),'-k','linewidth',0.5,'Color',[0 0 0])
	    if ~isempty(dsge_prior_weight)
	      plot(1:nirfs,DistribIRF_dsgevar(:,j,i,k),'-k','linewidth',0.5,'Color',[0.80 0.80 0.80])
	    end
	  end
	  plot(1:nirfs,MeanIRF_dsge(:,index(j),i),'-k','linewidth',3,'Color',[0 0 0])
	  if ~isempty(dsge_prior_weight)
	    plot(1:nirfs,MeanIRF_dsgevar(:,j,i),'-k','linewidth',3,'Color',[0.80 0.80 0.80])
	  end
	  xlim([1 nirfs]);
	  hold off
	  name = deblank(lgy_(SelecVariables(index(j)),:));
	  NAMES = strvcat(NAMES,name);
	  if TeX
	    texname = deblank(lgy_TeX_(SelecVariables(index(j)),:));
	    TEXNAMES = strvcat(TEXNAMES,['$' texname '$']);
	  end
	  title(name,'Interpreter','none')
	end
	if isempty(dsge_prior_weight)
	  eval(['print -depsc2 ' fname_ '_Bayesian_IRFdsge_' deblank(tit(i,:))]);
	  eval(['print -dpdf ' fname_  '_Bayesian_IRFdsge_' deblank(tit(i,:))]);
	  saveas(hh,[fname_  '_Bayesian_IRFdsge_' deblank(tit(i,:)) '.fig']);
	else
	  eval(['print -depsc2 ' fname_ '_Bayesian_IRFbvardsge_' deblank(tit(i,:))]);
	  eval(['print -dpdf ' fname_  '_Bayesian_IRFbvardsge_' deblank(tit(i,:))]);
	  saveas(hh,[fname_  '_Bayesian_IRFbvardsge_' deblank(tit(i,:)) '.fig']);	  
	end
	if options_.nograph, close(hh), end
	if TeX
	  fprintf(fidTeX,'\\begin{figure}[H]\n');
	  for jj = 1:number_of_plots_to_draw
	    fprintf(fidTeX,['\\psfrag{%s}[1][][0.5][0]{%s}\n'],deblank(NAMES(jj,:)),deblank(TEXNAMES(jj,:)));
	  end    
	  fprintf(fidTeX,'\\centering \n');
	  if ~~isempty(dsge_prior_weight)
	    fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_Bayesian_IRFdsge_%s}\n',fname_,deblank(tit(i,:)));
	  else
	    fprintf(fidTeX,['\\includegraphics[scale=0.5]{%s_Bayesian_IRFbvardsge_%s}\n',fname_,deblank(tit(i,:))]);
	  end  
	  if options_.relative_irf
	    fprintf(fidTeX,['\\caption{Bayesian relative IRF.}']);
	  else
	    if ~~isempty(dsge_prior_weight)
	      fprintf(fidTeX,'\\caption{Bayesian IRF (DSGE model).}');
	    else
	      fprintf(fidTeX,'\\caption{Bayesian IRF (BVAR-DSGE model).}');
	    end
	  end
	  if ~~isempty(dsge_prior_weight)
	    fprintf(fidTeX,'\\label{Fig:BayesianIRFdsge:%s}\n',deblank(tit(i,:)));
	  else
	    fprintf(fidTeX,'\\label{Fig:BayesianIRFbvardsge:%s}\n',deblank(tit(i,:)));
	  end  
	  fprintf(fidTeX,'\\end{figure}\n');
	  fprintf(fidTeX,' \n');
	end    
      elseif nbplt > 1
	for fig = 1:nbplt-1
	  if options_.relative_irf
	    hh = figure('Name',['Relative response to orthogonalized' ...
				' shock to ' tit(i,:) ' figure ' int2str(fig) '.']);
	  else
	    hh = figure('Name',['Orthogonalized shock to ' tit(i,:) ...
				' figure ' int2str(fig) '.']);
	  end
	  NAMES = [];
	  if TeX; TEXNAMES = []; end;
	  for j=1:nstar
	    jj = (fig-1)*nstar + j;
	    subplot(nr,nc,j);
	    plot([1 nirfs],[0 0],'-r','linewidth',0.5)
	    hold on
	    for k = 1:9
	      plot(1:options_.irf,DistribIRF_dsge(:,index(jj),i,k),'-k','linewidth',0.5,'Color',[0 0 0])
	      if ~isempty(dsge_prior_weight)
		plot(1:nirfs,DistribIRF_dsgevar(:,jj,i,k),'-k','linewidth',0.5,'Color',[0.80 0.80 0.80])
	      end
	    end
	    plot(1:nirfs,MeanIRF_dsge(:,index(jj),i),'-k','linewidth',3,'Color',[0 0 0])
	    if ~isempty(dsge_prior_weight)
	      plot(1:nirfs,MeanIRF_dsgevar(:,jj,i),'-k','linewidth',3,'Color',[0.80 0.80 0.80])
	    end
	    xlim([1 nirfs]);
	    hold off
	    name = deblank(lgy_(SelecVariables(index(jj)),:));
	    NAMES = strvcat(NAMES,name);
	    if TeX
	      texname = deblank(lgy_TeX_(SelecVariables(index(jj)),:));
	      TEXNAMES   = strvcat(TEXNAMES,['$' texname '$']);
	    end
	    title(name,'Interpreter','none')
	  end
	  if ~~isempty(dsge_prior_weight)
	    eval(['print -depsc2 ' fname_ '_Bayesian_IRFdsge_' deblank(tit(i,:)) int2str(fig)]);
	    eval(['print -dpdf ' fname_  '_Bayesian_IRFdsge_' deblank(tit(i,:)) int2str(fig)]);
	    saveas(hh,[fname_  '_Bayesian_IRFdsge_' deblank(tit(i,:)) int2str(fig) '.fig']);
	  else
	    eval(['print -depsc2 ' fname_ '_Bayesian_IRFbvardsge_' deblank(tit(i,:)) int2str(fig)]);
	    eval(['print -dpdf ' fname_  '_Bayesian_IRFbvardsge_' deblank(tit(i,:)) int2str(fig)]);
	    saveas(hh,[fname_  '_Bayesian_IRFbvardsge_' deblank(tit(i,:)) int2str(fig) '.fig']);
	  end
	  if options_.nograph, close(hh), end
	  if TeX
	    fprintf(fidTeX,'\\begin{figure}[H]\n');
	    for jj = 1:nstar
	      fprintf(fidTeX,['\\psfrag{%s}[1][][0.5][0]{%s}\n'],deblank(NAMES(jj,:)),deblank(TEXNAMES(jj,:)));
	    end    
	    fprintf(fidTeX,'\\centering \n');
	    if ~~isempty(dsge_prior_weight)
	      fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_Bayesian_IRFdsge_%s%s}\n',fname_,deblank(tit(i,:)),int2str(fig));
	    else
	      fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_Bayesian_IRFbvardsge_%s%s}\n',fname_,deblank(tit(i,:)),int2str(fig));
	    end
	    if options_.relative_irf == 1
	      fprintf(fidTeX,['\\caption{Bayesian relative IRF.}']);
	    else
	      if ~~isempty(dsge_prior_weight)
		fprintf(fidTeX,'\\caption{Bayesian IRF (DSGE model).}');
	      else
		fprintf(fidTeX,['\\caption{Bayesian IRF (BVAR-DSGE model).}']);
	      end
	    end
	    if ~~isempty(dsge_prior_weight)
	      fprintf(fidTeX,'\\label{Fig:BayesianIRFdsge:%s:%s}\n',deblank(tit(i,:)), int2str(fig));
	    else
	      fprintf(fidTeX,'\\label{Fig:BayesianIRFbvardsge:%s:%s}\n',deblank(tit(i,:)), int2str(fig));
	    end
	    fprintf(fidTeX,'\\end{figure}\n');
	    fprintf(fidTeX,' \n');
	  end    
	end
	hh = figure('Name',['Orthogonalized shock to ' tit(i,:) ' figure ' int2str(nbplt) '.']);
	NAMES = [];
	if TeX; TEXNAMES = []; end;
	for j=1:number_of_plots_to_draw -(nbplt-1)*nstar
	  jj = (nbplt-1)*nstar + j;
	  subplot(nr,nc,j);
	  plot([1 nirfs],[0 0],'-r','linewidth',0.5);
	  hold on
	  for k = 1:9
	    plot(1:options_.irf,DistribIRF_dsge(:,index(jj),i,k),'-k','linewidth',0.5,'Color',[0 0 0])
	    if ~isempty(dsge_prior_weight)
	      plot(1:nirfs,DistribIRF_dsgevar(:,jj,i,k),'-k','linewidth',0.5,'Color',[0.80 0.80 0.80])
	    end
	  end
	  plot(1:nirfs,MeanIRF_dsge(:,index(jj),i),'-k','linewidth',3,'Color',[0 0 0])
	  if ~isempty(dsge_prior_weight)
	    plot(1:nirfs,MeanIRF_dsgevar(:,jj,i),'-k','linewidth',3,'Color',[0.80 0.80 0.80])
	  end  
	  xlim([1 nirfs]);
	  hold off
	  name = deblank(lgy_(SelecVariables(index(jj)),:));
	  NAMES = strvcat(NAMES,name);
	  if TeX
	    texname = deblank(lgy_TeX_(SelecVariables(index(jj)),:));
	    TEXNAMES   = strvcat(TEXNAMES,['$' texname '$']);
	  end
	  title(name,'Interpreter','none')
	end
	if isempty(dsge_prior_weight)
	  eval(['print -depsc2 ' fname_ '_Bayesian_IRFdsge_' deblank(tit(i,:)) int2str(nbplt)]);
	  eval(['print -dpdf ' fname_  '_Bayesian_IRFdsge_' deblank(tit(i,:)) int2str(nbplt)]);
	  saveas(hh,[fname_  '_Bayesian_IRFdsge_' deblank(tit(i,:)) int2str(nbplt) '.fig']);
	else
	  eval(['print -depsc2 ' fname_ '_Bayesian_IRFbvardsge_' deblank(tit(i,:)) int2str(nbplt)]);
	  eval(['print -dpdf ' fname_  '_Bayesian_IRFbvardsge_' deblank(tit(i,:)) int2str(nbplt)]);
	  saveas(hh,[fname_  '_Bayesian_IRFbvardsge_' deblank(tit(i,:)) int2str(nbplt) '.fig']);	  
	end  
	  if options_.nograph, close(hh), end
	if TeX
	  fprintf(fidTeX,'\\begin{figure}[H]\n');
	  for jj = 1:nstar
	    fprintf(fidTeX,['\\psfrag{%s}[1][][0.5][0]{%s}\n'],deblank(NAMES(jj,:)),deblank(TEXNAMES(jj,:)));
	  end    
	  fprintf(fidTeX,'\\centering \n');
	  if ~~isempty(dsge_prior_weight)
	    fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_Bayesian_IRFdsge_%s%s}\n',fname_,deblank(tit(i,:)),int2str(nbplt));
	    fprintf(fidTeX,'\\caption{Bayesian IRF (DSGE model).}');
	    fprintf(fidTeX,'\\label{Fig:BayesianIRFdsge:%s:%s}\n',deblank(tit(i,:)), int2str(nbplt));
	  else
	    fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_Bayesian_IRFbvardsge_%s%s}\n',fname_,deblank(tit(i,:)),int2str(nbplt));
	    fprintf(fidTeX,'\\caption{Bayesian IRF (BVAR-DSGE model).}');
	    fprintf(fidTeX,'\\label{Fig:BayesianIRF:%s:%s}\n',deblank(tit(i,:)), int2str(nbplt));
	  end
	  fprintf(fidTeX,'\\end{figure}\n');
	  fprintf(fidTeX,' \n');
	end    
      else % nbplt = 0
	disp('There''s nothing to plot here!')
      end
    end
    if TeX
      fprintf(fidTeX,'%% End of TeX file.\n');
      fclose(fidTeX);
    end
    fprintf('MH: Posterior IRFs, done!\n');
  end
  %%
  %%    Posterior theoretical moments. Instead of displaying the posterior moments 
  %%  associated to the posterior mean of the structural parameters (by calling 
  %%  stoch_simul after estimation), metropolis.m will display the posterior mean 
  %%  of the theoretical moments and the posterior HPD intervals of theoretical
  %%  moments. All the results are saved in the global structure oo_ (posterior 
  %%  medians, posterior standard deviations and posterior deciles are also
  %%    computed and saved).
  %%
  if ~isempty(options_.unit_root_vars)
    vartan = []; 
    for i=1:nvar
      if isempty(strmatch(deblank(varlist(i,:)),options_.unit_root_vars,'exact'))       
	vartan = strvcat(vartan,varlist(i,:));
      end
    end
  else
    vartan = varlist;
  end
  if options_.moments_varendo & ~isempty(vartan)
    deep = MU;
    subindx = subset();
    nvar    = size(vartan,1);
    ivar = zeros(nvar,1);
    for i = 1:nvar
      ivar(i) = strmatch(vartan(i,:),lgy_,'exact');
    end
    nar = options_.ar;
    if B <= MAX_nthm1
      stock_thm1 = zeros(nvar,B);
    elseif B > MAX_nthm1
      stock_thm1 = zeros(nvar,MAX_nthm1);
    end
    if B <= MAX_nthm2
      stock_thm2 = zeros(nvar,nvar,B);
    elseif B > MAX_nthm2
      stock_thm2 = zeros(nvar,nvar,MAX_nthm2);
    end
    if B <= MAX_nthm3
      stock_thm3 = zeros(nvar,exo_nbr,B);
    elseif B > MAX_nthm3
      stock_thm3 = zeros(nvar,exo_nbr,MAX_nthm3);
    end
    if B <= MAX_nthm4
      stock_thm4 = zeros(nvar,nar,B);
    elseif B > MAX_nthm4
      stock_thm4 = zeros(nvar,nar,MAX_nthm4);
    end
    h = waitbar(0,'Posterior theoretical moments...');
    if nfile-ffil+1>1
      sfil_thm1 = 0;
      irun_thm1 = 0;
      sfil_thm2 = 0;
      irun_thm2 = 0;
      sfil_thm3 = 0;
      irun_thm3 = 0;
      sfil_thm4 = 0;
      irun_thm4 = 0;
      for b = 1:B;
	irun_thm1 = irun_thm1+1;
	irun_thm2 = irun_thm2+1;
	irun_thm3 = irun_thm3+1;
	irun_thm4 = irun_thm4+1;
	choose_an_mh_file = rand;
	mh_file_number = ...
	    FLN(find(choose_an_mh_file>=FLN(:,3)),1);
	if isempty(mh_file_number)
	  mh_file_number = ffil;
	else    
	  mh_file_number = mh_file_number(1);
	end    
	eval(['load ' instr1 int2str(mh_file_number) instr2]);
	clear post2 logpo2;
	DEEP  = x2(floor(rand*FLN(find(mh_file_number == FLN(:,1)),2))+1,:);
	deep(subindx) = DEEP(subindx);
	set_parameters(deep);
	dr_ = resol(ys_,0);
	Gamma_y = th_autocovariances(dr_,ivar);
	if options_.order == 2
	  m_mean = dr_.ys(ivar) + Gamma_y{options_.ar+3};
	else
	  m_mean = dr_.ys(ivar);
	end
	variance =  Gamma_y{1};
	if irun_thm1 < MAX_nthm1
	  stock_thm1(:,irun_thm1) = m_mean;
	else
	  stock_thm1(:,irun_thm1) = m_mean;
	  sfil_thm1 = sfil_thm1 + 1;
	  instr = [fname_ '_thm1' int2str(sfil_thm1) ' stock_thm1;'];
	  eval(['save ' instr]);
	  irun_thm1 = 0;
	  stock_thm1 = zeros(nvar,MAX_nthm1);
	end
	if irun_thm2 < MAX_nthm2
	  stock_thm2(:,:,irun_thm2) = variance;
	else
	  stock_thm2(:,:,irun_thm2) = variance;
	  sfil_thm2 = sfil_thm2 + 1;
	  instr = [fname_ '_thm2' int2str(sfil_thm2) ' stock_thm2;'];
	  eval(['save ' instr]);
	  irun_thm2 = 0;
	  stock_thm2 = zeros(nvar,nvar,MAX_nthm2);
	end
	if irun_thm3 < MAX_nthm3
	  stock_thm3(:,:,irun_thm3) = Gamma_y{nar+2};
	else
	  stock_thm3(:,:,irun_thm3) = Gamma_y{nar+2};
	  sfil_thm3 = sfil_thm3 + 1;
	  instr = [fname_ '_thm3' int2str(sfil_thm3) ' stock_thm3;'];
	  eval(['save ' instr]);
	  irun_thm3 = 0;
	  stock_thm3 = zeros(nvar,exo_nbr,MAX_nthm3);
	end
	if irun_thm4 < MAX_nthm4
	  for lag = 1:nar
	    stock_thm4(:,lag,irun_thm4) = diag(Gamma_y{1+lag});
	  end	
	else
	  for lag = 1:nar
	    stock_thm4(:,lag,irun_thm4) = diag(Gamma_y{1+lag});
	  end	
	  sfil_thm4 = sfil_thm4 + 1;
	  instr = [fname_ '_thm4' int2str(sfil_thm4) ' stock_thm4;'];
	  eval(['save ' instr]);
	  irun_thm4 = 0;
	  stock_thm4 = zeros(nvar,nar,MAX_nthm4);
	end
	waitbar(b/B,h);    
      end
      clear m_mean variance Gamma_y;
      if irun_thm1
	stock_thm1 = stock_thm1(:,1:irun_thm1);
	sfil_thm1 = sfil_thm1 + 1;
	instr = [fname_ '_thm1' int2str(sfil_thm1) ' stock_thm1;'];
	eval(['save ' instr]);
      end
      clear stock_thm1;
      if irun_thm2
	stock_thm2 = stock_thm2(:,:,1:irun_thm2);
	sfil_thm2 = sfil_thm2 + 1;
	instr = [fname_ '_thm2' int2str(sfil_thm2) ' stock_thm2;'];
	eval(['save ' instr]);
      end
      clear stock_thm2;
      if irun_thm3
	stock_thm3 = stock_thm3(:,:,1:irun_thm3);
	sfil_thm3 = sfil_thm3 + 1;
	instr = [fname_ '_thm3' int2str(sfil_thm3) ' stock_thm3;'];
	eval(['save ' instr]);
      end
      clear stock_thm3;
      if irun_thm4
	stock_thm4 = stock_thm4(:,:,1:irun_thm4);
	sfil_thm4 = sfil_thm4 + 1;
	instr = [fname_ '_thm4' int2str(sfil_thm4) ' stock_thm4;'];
	eval(['save ' instr]);
      end
      clear stock_thm4;
    else        
      sfil_thm1 = 0;
      irun_thm1 = 0;
      sfil_thm2 = 0;
      irun_thm2 = 0;
      sfil_thm3 = 0;
      irun_thm3 = 0;
      sfil_thm4 = 0;
      irun_thm4 = 0;
      eval(['load ' instr1 int2str(ffil) instr2]);
      NumberOfSimulations = length(logpo2);
      clear post2 logpo2;
      ivar1 = find(ismember(ivar,bayestopt_.i_var_stable));
      ivar1 = ivar(ivar1);
      for b = 1:B;
	irun_thm1 = irun_thm1+1;
	irun_thm2 = irun_thm2+1;
	irun_thm3 = irun_thm3+1;
	irun_thm4 = irun_thm4+1;
	DEEP  = x2(floor(rand*NumberOfSimulations)+1,:);
	deep(subindx) = DEEP(subindx);
	set_parameters(deep);
	dr_ = resol(ys_,0);
	Gamma_y = th_autocovariances(dr_,ivar1);
	if options_.order == 2
	  m_mean = dr_.ys(ivar) + Gamma_y{options_.ar+3};
	else
	  m_mean = dr_.ys(ivar);
	end
	variance = Gamma_y{1};
	if irun_thm1 < MAX_nthm1
	  stock_thm1(:,irun_thm1) = m_mean;
	else
	  stock_thm1(:,irun_thm1) = m_mean;
	  sfil_thm1 = sfil_thm1 + 1;
	  instr = [fname_ '_thm1' int2str(sfil_thm1) ' stock_thm1;'];
	  eval(['save ' instr]);
	  irun_thm1 = 0;
	  stock_thm1 = zeros(nvar,MAX_nthm1);
	end
	if irun_thm2 < MAX_nthm2
	  stock_thm2(:,:,irun_thm2) = variance;
	else
	  stock_thm2(:,:,irun_thm2) = variance;
	  sfil_thm2 = sfil_thm2 + 1;
	  instr = [fname_ '_thm2' int2str(sfil_thm2) ' stock_thm2;'];
	  eval(['save ' instr]);
	  irun_thm2 = 0;
	  stock_thm2 = zeros(nvar,nvar,MAX_nthm2);
	end
	if irun_thm3 < MAX_nthm3
	  stock_thm3(:,:,irun_thm3) = Gamma_y{nar+2};
	else
	  stock_thm3(:,:,irun_thm3) = Gamma_y{nar+2};
	  sfil_thm3 = sfil_thm3 + 1;
	  instr = [fname_ '_thm3' int2str(sfil_thm3) ' stock_thm3;'];
	  eval(['save ' instr]);
	  irun_thm3 = 0;
	  stock_thm3 = zeros(nvar,exo_nbr,MAX_nthm3);
	end
	if irun_thm4 < MAX_nthm4
	  for lag = 1:nar
	    stock_thm4(:,lag,irun_thm4) = diag(Gamma_y{1+lag});
	  end	
	else
	  for lag = 1:nar
	    stock_thm4(:,lag,irun_thm4) = diag(Gamma_y{1+lag});
	  end	
	  sfil_thm4 = sfil_thm4 + 1;
	  instr = [fname_ '_thm4' int2str(sfil_thm4) ' stock_thm4;'];
	  eval(['save ' instr]);
	  irun_thm4 = 0;
	  stock_thm4 = zeros(nvar,nar,MAX_nthm4);
	end
	waitbar(b/B,h);    
      end
      clear m_mean variance Gamma_y;
      if irun_thm1
	stock_thm1 = stock_thm1(:,1:irun_thm1);
	sfil_thm1 = sfil_thm1 + 1;
	instr = [fname_ '_thm1' int2str(sfil_thm1) ' stock_thm1;'];
	eval(['save ' instr]);
      end
      clear stock_thm1;
      if irun_thm2
	stock_thm2 = stock_thm2(:,:,1:irun_thm2);
	sfil_thm2 = sfil_thm2 + 1;
	instr = [fname_ '_thm2' int2str(sfil_thm2) ' stock_thm2;'];
	eval(['save ' instr]);
      end
      clear stock_thm2;
      if irun_thm3
	stock_thm3 = stock_thm3(:,:,1:irun_thm3);
	sfil_thm3 = sfil_thm3 + 1;
	instr = [fname_ '_thm3' int2str(sfil_thm3) ' stock_thm3;'];
	eval(['save ' instr]);
      end
      clear stock_thm3;
      if irun_thm4
	stock_thm4 = stock_thm4(:,:,1:irun_thm4);
	sfil_thm4 = sfil_thm4 + 1;
	instr = [fname_ '_thm4' int2str(sfil_thm4) ' stock_thm4;'];
	eval(['save ' instr]);
      end
      clear stock_thm4;
    end     
    close(h)
    %%
    %%  Now i compute some statistics (mean, median, std, deciles, HPD intervals)
    %%
    MeanMean = zeros(nvar,1);
    MedianMean = zeros(nvar,1);
    StdMean = zeros(nvar,1);
    DistribMean = zeros(nvar,9);
    HPDMean = zeros(nvar,2);
    tmp = zeros(B,1);
    for i = 1:nvar
      StartLine = 0;
      for file = 1:sfil_thm1 
	instr = [fname_ '_thm1' int2str(file)];
	eval(['load ' instr]);
	DeProfundis = size(stock_thm1,2);
	tmp(StartLine+1:StartLine+DeProfundis) = squeeze(stock_thm1(i,:));
	StartLine = StartLine+DeProfundis;
      end
      tmp = sort(tmp);
      MeanMean(i) = mean(tmp);
      MedianMean(i) = tmp(round(B*0.5));
      StdMean(i) = std(tmp);
      DistribMean(i,:) = reshape(tmp(deciles),1,9);
      tt = floor(options_.mh_conf_sig*B);
      a = 1; 
      b = tt;
      tmp2 = [1;tt;tmp(tt)-tmp(1)];
      while b <= B
	tmp1 = [a;b;tmp(b)-tmp(a)];
	a = a + 1;
	b = b + 1;
	if tmp1(3,1) < tmp2(3,1)
	  tmp2 = tmp1;     
	end    
      end
      HPDMean(i,1) = tmp(tmp2(1,1));
      HPDMean(i,2) = tmp(tmp2(2,1));
    end
    disp(' ')
    disp(' ')
    disp('POSTERIOR THEORETICAL EXPECTATION')
    disp(' ')
    titre = sprintf('%15s \t %6s \t %6s \t %6s \t %6s \t %6s\n',...
            'Variables',...
            'mean  ',...
            'median',...
            'std   ',...
            'HPDinf',...
            'HPDsup');
    disp(titre)
    for i=1:nvar
      disp(sprintf('%15s \t %6.3f \t %6.3f \t %6.3f \t %6.3f \t %6.3f',...
           deblank(lgy_(ivar(i),:)), ...
           MeanMean(i),...
           MedianMean(i),...
           StdMean(i),...
           HPDMean(i,1),...
           HPDMean(i,2)));
      eval(['oo_.PosteriorTheoreticalMoment.Expectation.Mean.' deblank(lgy_(ivar(i),:)) ' = MeanMean(i);']);
      eval(['oo_.PosteriorTheoreticalMoment.Expectation.Median.' deblank(lgy_(ivar(i),:)) ' = MedianMean(i);']);
      eval(['oo_.PosteriorTheoreticalMoment.Expectation.Std.' deblank(lgy_(ivar(i),:)) ' = StdMean(i);']);
      eval(['oo_.PosteriorTheoreticalMoment.Expectation.HPDinf.' deblank(lgy_(ivar(i),:)) ' = HPDMean(i,1);']);
      eval(['oo_.PosteriorTheoreticalMoment.Expectation.HPDsup.' deblank(lgy_(ivar(i),:)) ' = HPDMean(i,2);']);
    end
    if TeX
      fidTeX = fopen([fname_ '_PosteriorTheoreticalExpectation.TeX'],'w');
      fprintf(fidTeX,'%% TeX-table generated by metropolis.m (Dynare).\n');
      fprintf(fidTeX,['%% ' datestr(now,0)]);
      fprintf(fidTeX,' \n');
      fprintf(fidTeX,' \n');
      fprintf(fidTeX,'{\\tiny \n');
      fprintf(fidTeX,'\\begin{table}[H]\n');
      fprintf(fidTeX,'\\centering\n');
      fprintf(fidTeX,'\\begin{tabular}{l|ccccc} \n');
      fprintf(fidTeX,'\\hline\\hline \\\\ \n');
      fprintf(fidTeX,' Variables & mean & median  & std & HPD inf & HPD sup  \\\\ \n');
      fprintf(fidTeX,'\\hline \\\\ \n');
      for i=1:nvar
	fprintf(fidTeX,' $%s$ & %6.3f & %6.3f & %6.3f & %6.3f & %6.3f \\\\ \n',...
		deblank(lgy_TeX_(ivar(i),:)), ...
		MeanMean(i),...
		MedianMean(i),...
		StdMean(i),...
		HPDMean(i,1),...
		HPDMean(i,2));
      end   
      fprintf(fidTeX,'\\hline\\hline \n');
      fprintf(fidTeX,'\\end{tabular}\n ');    
      fprintf(fidTeX,'\\caption{Posterior theoretical expectation.}\n ');
      fprintf(fidTeX,'\\label{Table:PosteriorTheoreticalExpectation}\n');
      fprintf(fidTeX,'\\end{table}\n');
      fprintf(fidTeX,'} \n');
      fprintf(fidTeX,'%% End of TeX file.\n');
      fclose(fidTeX);
    end 
    MeanVariance = zeros(nvar,nvar,1);
    MedianVariance = zeros(nvar,nvar,1);
    StdVariance = zeros(nvar,nvar,1);
    DistribVariance = zeros(nvar,nvar,9);
    HPDVariance = zeros(nvar,nvar,2);
    for i = 1:nvar
      for j=1:nvar
	StartLine = 0;
	tmp = zeros(B,1);
	for file = 1:sfil_thm2 
	  instr = [fname_ '_thm2' int2str(file)];
	  eval(['load ' instr]);
	  DeProfundis = size(stock_thm2,3);
	  tmp(StartLine+1:StartLine+DeProfundis) = squeeze(stock_thm2(i,j,:));
	  StartLine = StartLine+DeProfundis;
	end 
	tmp = sort(tmp);
	MeanVariance(i,j) = mean(tmp);
	MedianVariance(i,j) = tmp(round(B*0.5));
	StdVariance(i,j) = std(tmp);
	DistribVariance(i,j,:) = reshape(tmp(deciles),1,1,9);
	tt = floor(options_.mh_conf_sig*B);
	a = 1; 
	b = tt;
	tmp2 = [1;tt;tmp(tt)-tmp(1)];
	while b <= B
	  tmp1 = [a;b;tmp(b)-tmp(a)];
	  a = a + 1;
	  b = b + 1;
	  if tmp1(3,1) < tmp2(3,1)
	    tmp2 = tmp1;     
	  end    
	end
	HPDVariance(i,j,1) = tmp(tmp2(1,1));
	HPDVariance(i,j,2) = tmp(tmp2(2,1));
      end   
    end
    disp(' ')
    disp(' ')
    disp('POSTERIOR THEORETICAL VARIANCES AND COVARIANCES')
    disp(' ')
    titre = sprintf('%15s \t %15s \t %9s \t %9s \t %9s \t %9s \t %9s\n',...
            'Variables',...
            'Variables',...
            'mean  ',...
            'median',...
            'std   ',...
            'HPDinf',...
            'HDPsup');
    disp(titre)
    for i=1:nvar
      for j=i:nvar
	disp(sprintf('%15s \t %15s \t %9.3g \t %9.3g \t %9.3g \t %9.3g \t %9.3g',...
		     deblank(lgy_(ivar(i),:)), ...
		     deblank(lgy_(ivar(j),:)), ...
		     MeanVariance(i,j),...
		     MedianVariance(i,j),...
		     StdVariance(i,j),...
		     HPDVariance(i,j,1),...
		     HPDVariance(i,j,2)));
	eval(['oo_.PosteriorTheoreticalMoment.Variance.Mean.' deblank(lgy_(ivar(i),:)) '_' deblank(lgy_(ivar(j),:)) ' = MeanVariance(i,j);']);
	eval(['oo_.PosteriorTheoreticalMoment.Variance.Median.' deblank(lgy_(ivar(i),:)) '_' deblank(lgy_(ivar(j),:)) ' = MedianVariance(i,j);']);
	eval(['oo_.PosteriorTheoreticalMoment.Variance.Std.' deblank(lgy_(ivar(i),:)) '_' deblank(lgy_(ivar(j),:)) ' = StdVariance(i,j);']);
	eval(['oo_.PosteriorTheoreticalMoment.Variance.HPDinf.' deblank(lgy_(ivar(i),:)) '_' deblank(lgy_(ivar(j),:)) ' = HPDVariance(i,j,1);']);
	eval(['oo_.PosteriorTheoreticalMoment.Variance.HPDsup.' deblank(lgy_(ivar(i),:)) '_' deblank(lgy_(ivar(j),:)) ' = HPDVariance(i,j,2);']);
      end       
    end
    if TeX
      fidTeX = fopen([fname_ '_PosteriorTheoreticalVariance.TeX'],'w');
      fprintf(fidTeX,'%% TeX-table generated by metropolis.m (Dynare).\n');
      fprintf(fidTeX,['%% ' datestr(now,0)]);
      fprintf(fidTeX,' \n');
      fprintf(fidTeX,' \n');
      fprintf(fidTeX,'{\\tiny \n');
      fprintf(fidTeX,'\\begin{table}[H]\n');
      fprintf(fidTeX,'\\centering\n');
      fprintf(fidTeX,'\\begin{tabular}{ll|ccccc} \n');
      fprintf(fidTeX,'\\hline\\hline \\\\ \n');
      fprintf(fidTeX,' Variables & Variables & mean & median  & std & HPD inf & HPD sup  \\\\ \n');
      fprintf(fidTeX,'\\hline \\\\ \n');
      for i=1:nvar
	for j=i:nvar
	  fprintf(fidTeX,' $%s$ & $%s$ & %6.3f & %6.3f & %6.3f & %6.3f & %6.3f \\\\ \n',...
		  deblank(lgy_TeX_(ivar(i),:)), ...
		  deblank(lgy_TeX_(ivar(j),:)), ...
		  MeanVariance(i,j),...
		  MedianVariance(i,j),...
		  StdVariance(i,j),...
		  HPDVariance(i,j,1),...
		  HPDVariance(i,j,2));
	end     
      end   
      fprintf(fidTeX,'\\hline\\hline \n');
      fprintf(fidTeX,'\\end{tabular}\n ');    
      fprintf(fidTeX,'\\caption{Posterior theoretical variances and covariances.}\n ');
      fprintf(fidTeX,'\\label{Table:PosteriorTheoreticalVariances}\n');
      fprintf(fidTeX,'\\end{table}\n');
      fprintf(fidTeX,'} \n');
      fprintf(fidTeX,'%% End of TeX file.\n');
      fclose(fidTeX);
    end 
    MeanCorrelation = zeros(nvar,nvar,1);
    MedianCorrelation = zeros(nvar,nvar,1);
    StdCorrelation = zeros(nvar,nvar,1);
    DistribCorrelation = zeros(nvar,nvar,9);
    HPDCorrelation = zeros(nvar,nvar,2);
    tmpp    = zeros(B,1);
    tmppp   = zeros(B,1);
    for i = 1:nvar
      for j=1:nvar
	StartLine = 0;
	tmp = zeros(B,1);
	for file = 1:sfil_thm2 
	  instr = [fname_ '_thm2' int2str(file)];
	  eval(['load ' instr]);
	  DeProfundis = size(stock_thm2,3);
	  tmp(StartLine+1:StartLine+DeProfundis) = squeeze(stock_thm2(i,j,:));
	  tmpp(StartLine+1:StartLine+DeProfundis) = squeeze(stock_thm2(i,i,:));
	  tmppp(StartLine+1:StartLine+DeProfundis) = squeeze(stock_thm2(j,j,:));
	  StartLine = StartLine+DeProfundis;
	end
	tmp = sort(tmp./sqrt(tmpp.*tmppp));
	MeanCorrelation(i,j) = mean(tmp);
	MedianCorrelation(i,j) = tmp(round(B*0.5));
	StdCorrelation(i,j) = std(tmp);
	DistribCorrelation(i,j,:) = reshape(tmp(deciles),1,1,9);
	tt = floor(options_.mh_conf_sig*B);
	a = 1; 
	b = tt;
	tmp2 = [1;tt;tmp(tt)-tmp(1)];
	while b <= B
	  tmp1 = [a;b;tmp(b)-tmp(a)];
	  a = a + 1;
	  b = b + 1;
	  if tmp1(3,1) < tmp2(3,1)
	    tmp2 = tmp1;     
	  end    
	end
	HPDCorrelation(i,j,1) = tmp(tmp2(1,1));
	HPDCorrelation(i,j,2) = tmp(tmp2(2,1));
      end   
    end
    clear tmpp tmppp;
    disp(' ')
    disp(' ')
    disp('POSTERIOR THEORETICAL CORRELATIONS')
    disp(' ')
    titre = sprintf('%15s \t %15s \t %6s \t %6s \t %6s \t %6s \t %6s\n',...
            'Variables',...
            'Variables',...
            'mean  ',...
            'median',...
            'std   ',...
            'HPDinf',...
            'HPDsup');
    disp(titre)
    for i=1:nvar-1
      for j=i+1:nvar
	disp(sprintf('%15s \t %15s \t %6.3f \t %6.3f \t %6.3f \t %6.3f \t %6.3f',...
		     deblank(lgy_(ivar(i),:)), ...
		     deblank(lgy_(ivar(j),:)), ...
		     MeanCorrelation(i,j),...
		     MedianCorrelation(i,j),...
		     StdCorrelation(i,j),...
		     HPDCorrelation(i,j,1),...
		     HPDCorrelation(i,j,2)));
	eval(['oo_.PosteriorTheoreticalMoment.Correlation.Mean.' deblank(lgy_(ivar(i),:)) '_' deblank(lgy_(ivar(j),:)) ' = MeanCorrelation(i,j);']);
	eval(['oo_.PosteriorTheoreticalMoment.Correlation.Median.' deblank(lgy_(ivar(i),:)) '_' deblank(lgy_(ivar(j),:)) ' = MedianCorrelation(i,j);']);
	eval(['oo_.PosteriorTheoreticalMoment.Correlation.Std.' deblank(lgy_(ivar(i),:)) '_' deblank(lgy_(ivar(j),:)) ' = StdCorrelation(i,j);']);
	eval(['oo_.PosteriorTheoreticalMoment.Correlation.HPDinf.' deblank(lgy_(ivar(i),:)) '_' deblank(lgy_(ivar(j),:)) ' = HPDCorrelation(i,j,1);']);
	eval(['oo_.PosteriorTheoreticalMoment.Correlation.HPDsup.' deblank(lgy_(ivar(i),:)) '_' deblank(lgy_(ivar(j),:)) ' = HPDCorrelation(i,j,2);']);
      end       
    end
    if TeX
      fidTeX = fopen([fname_ '_PosteriorTheoreticalCorrelation.TeX'],'w');
      fprintf(fidTeX,'%% TeX-table generated by metropolis.m (Dynare).\n');
      fprintf(fidTeX,['%% ' datestr(now,0)]);
      fprintf(fidTeX,' \n');
      fprintf(fidTeX,' \n');
      fprintf(fidTeX,'{\\tiny \n');
      fprintf(fidTeX,'\\begin{table}[H]\n');
      fprintf(fidTeX,'\\centering\n');
      fprintf(fidTeX,'\\begin{tabular}{ll|ccccc} \n');
      fprintf(fidTeX,'\\hline\\hline \\\\ \n');
      fprintf(fidTeX,' Variables & Variables & mean & median  & std & HPD inf & HPD sup  \\\\ \n');
      fprintf(fidTeX,'\\hline \\\\ \n');
      for i=1:nvar-1
	for j=i+1:nvar
	  fprintf(fidTeX,' $%s$ & $%s$ & %6.3f & %6.3f & %6.3f & %6.3f & %6.3f \\\\ \n',...
		  deblank(lgy_TeX_(ivar(i),:)), ...
		  deblank(lgy_TeX_(ivar(j),:)), ...
		  MeanCorrelation(i,j),...
		  MedianCorrelation(i,j),...
		  StdCorrelation(i,j),...
		  HPDCorrelation(i,j,1),...
		  HPDCorrelation(i,j,2));
	end
      end   
      fprintf(fidTeX,'\\hline\\hline \n');
      fprintf(fidTeX,'\\end{tabular}\n ');    
      fprintf(fidTeX,'\\caption{Posterior theoretical correlations.}\n ');
      fprintf(fidTeX,'\\label{Table:PosteriorTheoreticalCorrelations}\n');
      fprintf(fidTeX,'\\end{table}\n');
      fprintf(fidTeX,'} \n');
      fprintf(fidTeX,'%% End of TeX file.\n');
      fclose(fidTeX);
    end 
    MeanDecomp = zeros(nvar,exo_nbr,1);
    MedianDecomp = zeros(nvar,exo_nbr,1);
    StdDecomp = zeros(nvar,exo_nbr,1);
    DistribDecomp = zeros(nvar,exo_nbr,9);
    HPDDecomp = zeros(nvar,exo_nbr,2);
    for i = 1:nvar
      for j=1:exo_nbr
	StartLine = 0;
	for file = 1:sfil_thm3 
	  instr = [fname_ '_thm3' int2str(file)];
	  eval(['load ' instr]);
	  DeProfundis = size(stock_thm3,3);
	  tmp(StartLine+1:StartLine+DeProfundis) = squeeze(stock_thm3(i,j,:));
	  StartLine = StartLine+DeProfundis;
	end
	tmp = sort(tmp);
	MeanDecomp(i,j) = mean(tmp);
	MedianDecomp(i,j) = tmp(round(B*0.5));
	StdDecomp(i,j) = std(tmp);
	DistribDecomp(i,j,:) = reshape(tmp(deciles),1,1,9);
	tt = floor(options_.mh_conf_sig*B);
	a = 1; 
	b = tt;
	tmp2 = [1;tt;tmp(tt)-tmp(1)];
	while b <= B
	  tmp1 = [a;b;tmp(b)-tmp(a)];
	  a = a + 1;
	  b = b + 1;
	  if tmp1(3,1) < tmp2(3,1)
	    tmp2 = tmp1;     
	  end    
	end
	HPDDecomp(i,j,1) = tmp(tmp2(1,1));
	HPDDecomp(i,j,2) = tmp(tmp2(2,1));
      end   
    end
    disp(' ')
    disp(' ')
    disp('POSTERIOR THEORETICAL VARIANCE DECOMPOSITION')
    disp(' ')
    titre = sprintf('%15s \t %15s \t %6s \t %6s \t %6s \t %6s \t %6s\n',...
            'Variables',...
            'Sources',...
            'mean  ',...
            'median',...
            'std   ',...
            'HPDinf',...
            'HDPsup');
    disp(titre)
    lgx1(lgx_orig_ord_,:) = lgx_;
    for i=1:nvar
      for j=1:exo_nbr
	disp(sprintf('%15s \t %15s \t %6.3f \t %6.3f \t %6.3f \t %6.3f \t %6.3f',...
		     deblank(lgy_(ivar(i),:)), ...
		     deblank(lgx1(j,:)), ...
		     MeanDecomp(i,j),...
		     MedianDecomp(i,j),...
		     StdDecomp(i,j),...
		     HPDDecomp(i,j,1),...
		     HPDDecomp(i,j,2)));
	eval(['oo_.PosteriorTheoreticalMoment.Decomp.Mean.' deblank(lgy_(ivar(i),:)) '_' deblank(lgx1(j,:)) ' = MeanDecomp(i,j);']);
	eval(['oo_.PosteriorTheoreticalMoment.Decomp.Median.' deblank(lgy_(ivar(i),:)) '_' deblank(lgx1(j,:)) ' = MedianDecomp(i,j);']);
	eval(['oo_.PosteriorTheoreticalMoment.Decomp.Std.' deblank(lgy_(ivar(i),:)) '_' deblank(lgx1(j,:)) ' = StdDecomp(i,j);']);
	eval(['oo_.PosteriorTheoreticalMoment.Decomp.HPDinf.' deblank(lgy_(ivar(i),:)) '_' deblank(lgx1(j,:)) ' = HPDDecomp(i,j,1);']);
	eval(['oo_.PosteriorTheoreticalMoment.Decomp.HPDsup.' deblank(lgy_(ivar(i),:)) '_' deblank(lgx1(j,:)) ' = HPDDecomp(i,j,2);']);
      end       
    end
    if TeX
      fidTeX = fopen([fname_ '_PosteriorTheoreticalVarianceDecomposition.TeX'],'w');
      fprintf(fidTeX,'%% TeX-table generated by metropolis.m (Dynare).\n');
      fprintf(fidTeX,['%% ' datestr(now,0)]);
      fprintf(fidTeX,' \n');
      fprintf(fidTeX,' \n');
      fprintf(fidTeX,'{\\tiny \n');
      fprintf(fidTeX,'\\begin{table}[H]\n');
      fprintf(fidTeX,'\\centering\n');
      fprintf(fidTeX,'\\begin{tabular}{ll|ccccc} \n');
      fprintf(fidTeX,'\\hline\\hline \\\\ \n');
      fprintf(fidTeX,' Variables & Sources & mean & median  & std & HPD inf & HPD sup  \\\\ \n');
      fprintf(fidTeX,'\\hline \\\\ \n');
      lgx_TeX1(lgx_orig_ord_,:) = lgx_TeX_;
      for i=1:nvar
	for j=1:exo_nbr
	  fprintf(fidTeX,' $%s$ & $%s$ & %6.3f & %6.3f & %6.3f & %6.3f & %6.3f \\\\ \n',...
		  deblank(lgy_TeX_(ivar(i),:)), ...
		  deblank(lgx_TeX1(j,:)), ...
		  MeanDecomp(i,j),...
		  MedianDecomp(i,j),...
		  StdDecomp(i,j),...
		  HPDDecomp(i,j,1),...
		  HPDDecomp(i,j,2));
	end     
      end   
      fprintf(fidTeX,'\\hline\\hline \n');
      fprintf(fidTeX,'\\end{tabular}\n ');    
      fprintf(fidTeX,'\\caption{Posterior theoretical variance decomposition.}\n ');
      fprintf(fidTeX,'\\label{Table:PosteriorTheoreticalVarianceDecomposition}\n');
      fprintf(fidTeX,'\\end{table}\n');
      fprintf(fidTeX,'} \n');
      fprintf(fidTeX,'%% End of TeX file.\n');
      fclose(fidTeX);
    end 
    MeanAutoCorr = zeros(nvar,nar,1);
    MedianAutoCorr = zeros(nvar,nar,1);
    StdAutoCorr = zeros(nvar,nar,1);
    DistribAutoCorr = zeros(nvar,nar,9);
    HPDAutoCorr = zeros(nvar,nar,2);
    for i = 1:nvar
      for j=1:nar
	StartLine = 0;
	for file = 1:sfil_thm4 
	  instr = [fname_ '_thm4' int2str(file)];
	  eval(['load ' instr]);
	  DeProfundis = size(stock_thm4,3);
	  tmp(StartLine+1:StartLine+DeProfundis) = squeeze(stock_thm4(i,j,:));
	  StartLine = StartLine+DeProfundis;
	end
	tmp = sort(tmp);
	MeanAutoCorr(i,j) = mean(tmp);
	MedianAutoCorr(i,j) = tmp(round(B*0.5));
	StdAutoCorr(i,j) = std(tmp);
	DistribAutoCorr(i,j,:) = reshape(tmp(deciles),1,1,9);
	tt = floor(options_.mh_conf_sig*B);
	a = 1; 
	b = tt;
	tmp2 = [1;tt;tmp(tt)-tmp(1)];
	while b <= B
	  tmp1 = [a;b;tmp(b)-tmp(a)];
	  a = a + 1;
	  b = b + 1;
	  if tmp1(3,1) < tmp2(3,1)
	    tmp2 = tmp1;     
	  end    
	end
	HPDAutoCorr(i,j,1) = tmp(tmp2(1,1));
	HPDAutoCorr(i,j,2) = tmp(tmp2(2,1));
      end   
    end
    disp(' ')
    disp(' ')
    disp('POSTERIOR THEORETICAL AUTOCORRELATION')
    disp(' ')
    titre = sprintf('%15s \t %3s \t %6s \t %6s \t %6s \t %6s \t %6s\n',...
            'Variables',...
            'Lag',...
            'mean  ',...
            'median',...
            'std   ',...
            'HPDinf',...
            'HDPsup');
    disp(titre)
    for i=1:nvar
      for j=1:nar
	disp(sprintf('%15s \t %3s \t %6.3f \t %6.3f \t %6.3f \t %6.3f \t %6.3f',...
		     deblank(lgy_(ivar(i),:)), ...
		     [int2str(j) ' '], ...
		     MeanAutoCorr(i,j),...
		     MedianAutoCorr(i,j),...
		     StdAutoCorr(i,j),...
		     HPDAutoCorr(i,j,1),...
		     HPDAutoCorr(i,j,2)));
	eval(['oo_.PosteriorTheoreticalMoment.AutoCorrelation.Mean.' deblank(lgy_(ivar(i),:)) '_lag' int2str(j) ' = MeanAutoCorr(i,j);']);
	eval(['oo_.PosteriorTheoreticalMoment.AutoCorrelation.Median.' deblank(lgy_(ivar(i),:)) '_lag' int2str(j) ' = MedianAutoCorr(i,j);']);
	eval(['oo_.PosteriorTheoreticalMoment.AutoCorrelation.Std.' deblank(lgy_(ivar(i),:)) '_lag' int2str(j) ' = StdAutoCorr(i,j);']);
	eval(['oo_.PosteriorTheoreticalMoment.AutoCorrelation.HPDinf.' deblank(lgy_(ivar(i),:)) '_lag' int2str(j) ' = HPDAutoCorr(i,j,1);']);
	eval(['oo_.PosteriorTheoreticalMoment.AutoCorrelation.HPDsup.' deblank(lgy_(ivar(i),:)) '_lag' int2str(j) ' = HPDAutoCorr(i,j,2);']);
      end       
    end
    if TeX
      fidTeX = fopen([fname_ '_PosteriorTheoreticalAutocorrelation.TeX'],'w');
      fprintf(fidTeX,'%% TeX-table generated by metropolis.m (Dynare).\n');
      fprintf(fidTeX,['%% ' datestr(now,0)]);
      fprintf(fidTeX,' \n');
      fprintf(fidTeX,' \n');
      fprintf(fidTeX,'{\\tiny \n');
      fprintf(fidTeX,'\\begin{table}[H]\n');
      fprintf(fidTeX,'\\centering\n');
      fprintf(fidTeX,'\\begin{tabular}{ll|ccccc} \n');
      fprintf(fidTeX,'\\hline\\hline \\\\ \n');
      fprintf(fidTeX,' Variables & Lag & mean & median  & std & HPD inf & HPD sup  \\\\ \n');
      fprintf(fidTeX,'\\hline \\\\ \n');
      for i=1:nvar
	for j=1:nar
	  fprintf(fidTeX,' $%s$ & $%s$ & %6.3f & %6.3f & %6.3f & %6.3f & %6.3f \\\\ \n',...
		  deblank(lgy_TeX_(ivar(i),:)), ...
		  int2str(j), ...
		  MeanAutoCorr(i,j),...
		  MedianAutoCorr(i,j),...
		  StdAutoCorr(i,j),...
		  HPDAutoCorr(i,j,1),...
		  HPDAutoCorr(i,j,2));
	end     
      end   
      fprintf(fidTeX,'\\hline\\hline \n');
      fprintf(fidTeX,'\\end{tabular}\n ');    
      fprintf(fidTeX,'\\caption{Posterior theoretical auto-correlation.}\n ');
      fprintf(fidTeX,'\\label{Table:PosteriorTheoreticalAutoCorrelation}\n');
      fprintf(fidTeX,'\\end{table}\n');
      fprintf(fidTeX,'} \n');
      fprintf(fidTeX,'%% End of TeX file.\n');
      fclose(fidTeX);
    end 
  end % options_.moments_varendo
  
  %% Un dernier petit coup de DsgeLikelihood juste pour remettre les parametres
  %% structurels et la matrice de variance-covariance aux valeurs qui
  %% correspondent a la moyenne posterieure (en vue d'une utilisation éventuelle
  %% de stoch_simul après le Metropolis-Hastings).
  [lnprior,cost_flag,ys,trend_coeff] = DsgeLikelihood(post_mean,gend,data);
  %% Now I save the seeds (If the user wants to start another MH, he can start from the
  %% previous state of the random number generator by using the command "LoadPreviousSeed"
  %% before the estimation command)
  Seed.NormalDeviates  = randn('state');
  Seed.UniformDeviates = rand('state');
  save LastSeed Seed;
  %% That's All!


  % SA 08-18-2004       * Corrected a bug in forecasts (HPD intervals).
  %                 * metropolis.m now displays "true bayesian" smooth shocks. The mean
  %                 - across the metropolis draws - of the smooth shocks instead of the 
  %                 smooth shocks obtained from the posterior mean are displayed.
  %                 * Added "true bayesian" smooth measurement error.
  %                 * Added "true bayesian" smooth variables (all the variables in the 
  %                 state vector).
  %                 * Added deciles for the posterior distribution of the smooth shocks,
  %                 variables and measurement errors (green curves).
  % SA 08-19-2004       * Added posterior IRFs.
  % SA 08-21-2004       * Added posterior theoretical moments.
  % SA 08-23-2004     * Added correction to the modified harmonic mean estimator of the
  %                   log-marginal density. The variance of the weighting distribution
  %                   automatically increases if needed (to be revised).)
  % SA 12-02-2004       * Changed architecture for the global structure oo_ 
