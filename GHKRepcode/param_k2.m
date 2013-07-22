phi_ch = 0.04;   % weight of high productivity sector in consumption
phi_cl = 1-phi_ch;  % weight of low productivity sector in consumption

phi_eh = .85;     % weight of high productivity sector in equipment investment
phi_el = 1-phi_eh;  % weight of low productivity sector in equipement investment
phi_sh = 0.00001;   % weight of high productivity sector in structures investment
phi_sl = 1-phi_sh;  % weight of low productivity sector in structures investment

% change CD production function in H and L sectors
alpha_eh = 0.43; 
alpha_nh = 0.46; 

alpha_el = 0.15; 
alpha_nl = 0.72;

metanu = 0;
metaomega=100;  %makes capital predetermined (immobile in the first period).