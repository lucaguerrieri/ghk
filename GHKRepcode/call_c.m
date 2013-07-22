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
dynare ghk3    % 
save simulation1 dr_ ys_ lgy_ lgx_

!copy param_c.m param_reset.m
dynare ghk3    % 
save simulation2 dr_ ys_ lgy_ lgx_



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
         0.01         %eps_zeh
         0.01         %eps_zel
         0            %eps_zsh
         0]*2.565982108451157/0.9162321481851%*0.684153001841835;          %eps_zsl     
dset = 'f1';        
makeirf

load simulation2

shock = [0            %eps_ah   labor-augmenting mfp
         0            %eps_al
         0.02*0.877844730479991   %eps_gh   Hicks-neutral mfp
         0            %eps_gl
         0            %eps_zeh   
         0            %eps_zel
         0            %eps_zsh
         0]/ 0.8737704207027;          %eps_zsl

     
dset = 'f2';        
makeirf

%%
legendlist = char('MEI Shock','MFP Shock (Different production functions)');
figtitle = [''];
plotirf2linesoriginal