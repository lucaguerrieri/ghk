clear


% will hold ISP shock with no adj costs for investment 
!copy param_a2.m param_reset.m
dynare ghk3_prices    % 
save simulation1 dr_ ys_ lgy_ lgx_

% will hold MFP shock with only incomplete specialization and no adj. costs
% for investment
!copy param_h2.m param_reset.m
dynare ghk2_prices    % 
save simulation2 dr_ ys_ lgy_ lgx_

% will hold MFP shock with all departures from equivalence and no adj.
% costs for investment
!copy param_k2.m param_reset.m
dynare ghk2_prices    % 
save simulation3 dr_ ys_ lgy_ lgx_
      

% simulation horizon
nperiods = 5000;

nperiods_long = 200;
% shorter simulation horizon for plotting (needs to be less than or
% equal to nperiods)
nperiods_short = 20;


load simulation1 
shock =  [0            %eps_ah
         0            %eps_al
         0            %eps_gh
         0            %eps_gl
         0.01*9.800905865074714         %eps_zeh
         0.01*9.800905865074714         %eps_zel
         0            %eps_zsh
         0]* 0.261810709000957/0.9162321537865;%*1.108239152179862;          %eps_zsl     
dset = 'f1';        
makeirf

load simulation2

shock = [0            %eps_ah   labor-augmenting mfp
         0            %eps_al
         0.02/.5615*2.751579484874512   %eps_gh   Hicks-neutral mfp
         0            %eps_gl
         0            %eps_zeh   
         0            %eps_zel
         0            %eps_zsh
         0]*0.282969138755187/ 0.9301710793305;    
     
dset = 'f2';        
makeirf


load simulation3

shock = [0            %eps_ah   labor-augmenting mfp
         0            %eps_al
         0.02/.5615*2.751579484874512   %eps_gh   Hicks-neutral mfp
         0            %eps_gl
         0            %eps_zeh   
         0            %eps_zel
         0            %eps_zsh
         0]*0.282969138755187*0.760677891047739/0.9301710793305/0.8986536154680/ 1.0750710511445;%*1.382651915542789;          %eps_zsl
  %eps_zsl

     
dset = 'f3';        
makeirf

%%
legendlist = char('MEI shock in a one-sector model (aggregate equivalence)','MFP shock in the machinery sector (incomplete specialization)','MFP shock in the machinery sector (all departures)');
figtitle = [''];
plotirf3linesprices