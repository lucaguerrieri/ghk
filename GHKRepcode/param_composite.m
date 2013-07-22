sigma_h = prodmetasigma; %1.001; % production elasticity
sigma_l = prodmetasigma; %1.001; % production elasticity


sigma_c =  assemmetasigma;%1.001; % assembly elasticity
Asigma_e = assemmetasigma;%1.001; % assembly elasticity
Asigma_s = assemmetasigma;%1.001; % assembly elasticity


omega_eh = metaomega;
omega_sh = metaomega;
omega_el = metaomega;
omega_sl = metaomega;  

nu_eh = metanu;
nu_el = metanu;
nu_sh = metanu;
nu_sl = metanu;

rho_ah = metarho; 
rho_al = metarho; 
rho_gh = metarho;
rho_gl = metarho;
rho_zeh = metarho; 
rho_zel = metarho; 
rho_zsh = metarho; 
rho_zsl = metarho; 

% from the first-order conditions with respect to labor
n_h2y_h = alpha_nh;
n_l2y_l = alpha_nl;


% from the technology of production
alpha_sh=(1-alpha_nh-alpha_eh*(beta*g_e/(1-beta*(1-delta_eh)))^(sigma_h-1))/(beta*g_s/(1-beta*(1-delta_sh)))^(sigma_h-1);
alpha_sl=(1-alpha_nl-alpha_el*(beta*g_e/(1-beta*(1-delta_el)))^(sigma_l-1))/(beta*g_s/(1-beta*(1-delta_sl)))^(sigma_l-1);

% from the first-order conditions with respect to capital
k_eh2y_h = alpha_eh*(beta*g_e/(1-beta*(1-delta_eh)))^(sigma_h);
k_el2y_l = alpha_el*(beta*g_e/(1-beta*(1-delta_el)))^(sigma_l);
k_sh2y_h = alpha_sh*(beta*g_s/(1-beta*(1-delta_sh)))^(sigma_h);
k_sl2y_l = alpha_sl*(beta*g_s/(1-beta*(1-delta_sl)))^(sigma_l);

% sectoral investment ratios from the capital accumulation equations
s_eh = (delta_eh)*k_eh2y_h;
s_sh = (delta_sh)*k_sh2y_h;
s_el = (delta_el)*k_el2y_l;
s_sl = (delta_sl)*k_sl2y_l;


% define
a11= s_eh *(phi_ch- phi_eh)+ s_el *(phi_cl-phi_el); 
a12= s_eh *(phi_ch- phi_sh)+ s_el *(phi_cl-phi_sl);
a21= s_sh *(phi_ch- phi_eh)+ s_sl *(phi_cl-phi_el); 
a22= s_sh *(phi_ch- phi_sh)+ s_sl *(phi_cl-phi_sl);

tmp = inv([1+g_e^(-1)*a11, g_s^(-1)*a12;... 
          g_e^(-1)*a21 , 1+g_s^(-1)*a22 ])*[phi_ch*s_eh+phi_cl*s_el; phi_ch*s_sh + phi_cl*s_sl];
se=tmp(1);
ss=tmp(2);
% calculate the aggregate savings rate
s = g_e^(-1)*se+g_s^(-1)*ss; 

y_h2y = phi_ch-s*phi_ch+phi_eh*g_e^(-1)*se+phi_sh*g_s^(-1)*ss;
y_l2y = phi_cl*(1-s)+phi_el*g_e^(-1)*se+phi_sl*g_s^(-1)*ss;

% calculate labor allocation shares
n2y = n_h2y_h*y_h2y+n_l2y_l*y_l2y;
y2n = 1/n2y;
y = y2n*nbar;

check_k(y_h2y,y_l2y,y, k_eh2y_h, k_el2y_l, k_sh2y_h, k_sl2y_l)
