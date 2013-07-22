%%Designed for ISP shock. Perfect capital mobility so like one-good economy
%%identical CD production
%%assembly of consumption and of structures all l
%%assembly of equiqpment all h
%%perfect capital mobility
%%investment adjustment costs

% parameter values
beta = 0.99;
eta = 0;
gamma = 1.00001;    % governs intertemporal elasticity of substitution
phi_ch = 0.00001;   % weight of high productivity sector in consumption
phi_cl = 1-phi_ch;  % weight of low productivity sector in consumption

phi_eh = .9999;     % weight of high productivity sector in equipment investment
phi_el = 1-phi_eh;  % weight of low productivity sector in equipement investment
phi_sh = 0.00001;   % weight of high productivity sector in structures investment
phi_sl = 1-phi_sh;  % weight of low productivity sector in structures investment

% free parameter7765 for steady state
g_c = 1;
alpha_eh = 0.17;    % CD shares in high productivity sector
alpha_nh = 0.70; 

alpha_el = 0.17;    % CD shares in low productivity sector
alpha_nl = 0.70; 

% free parameter for steady state
g_e = 1;
g_s = 1;

% free parameter for steady state
nbar = 1;           % time spent working

% production elasticity (hooks to sigma_h and sigma_l)
prodmetasigma = 0.99;

% assembly elasticity (hooks to sigma_c, Asigma_e, Asigma_s)
assemmetasigma = 0.5;


% depreciation  rates
delta_eh = 0.124/4;
delta_el = 0.124/4;
delta_sh = 0.056/4;
delta_sl = 0.056/4;

% mobility of capital ( if close to zero, capital perfectly mobile
%                       if high value, capital predetermined)
% hooks to omega_eh, omega_sh, omega_el, omega_sl
metaomega = .0000001;

% hooks to nu_eh, nu_el, nu_sh, nu_sl
% investment adjustment costs
%metanu = 2;
metanu = 2;
%metanu = 0;
% hooks to all persistence coefficients for shocks' AR processes
metarho = 1.0;









