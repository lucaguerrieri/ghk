clear

%% Run model

%%description of param10
%%Model with identical CD production
%%assembly of consumption and of structures all l
%%assembly of equiqpment all h
%%perfect capital mobility
%%no investment adjustment costs
%%Hicks neutral Shock in h sector only
%%No captial moves because labor shifts t fully offset
%%effect of productivity increase
!copy param_a.m param_reset.m
dynare ghk3_prices    % NB to get equivalence, I had to modify the way in which
               % the adjustment cost on investment are specified
               % use ghk2 to go back to original adjustment costs
save simulation1 dr_ ys_ lgy_ lgx_

!copy param_a.m param_reset.m
dynare ghk2_prices    % 
save simulation2 dr_ ys_ lgy_ lgx_

!copy param_k.m param_reset.m
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
         0]* 0.261810709000957* 1.091426446974682;%*1.108239152179862;          %eps_zsl     
dset = 'f1';        
makeirf

load simulation2 
shock =  [0            %eps_ah
         0            %eps_al
         0            %eps_gh
         0            %eps_gl
         0.01*9.800905865074714         %eps_zeh
         0.01*9.800905865074714         %eps_zel
         0            %eps_zsh
         0]* 0.261810709000957* 1.091426446974682;%*1.108239152179862;          %eps_zsl     
dset = 'f2';        
makeirf


load simulation3

shock = [          0
                   0
   0.023494373177152
                   0
                   0
                   0
                   0
                   0];          %eps_zsl

     
dset = 'f3';        
makeirf



%%
legendlist = char('MEI shock in a one-sector model (aggregate equivalence)',...
                  'MEI shock in a one-sector model (inv. adj. costs in mixed units)',...
                  'MFP shock in the machinery sector (all departures)');
figtitle = [''];
plotirf3linesprices

