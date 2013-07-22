//Adjustment costs are "standard" which is taken to mean that they
//depend on physical not effective units. For in which adjustment costs
// depend on effective units see ghk_standardaltinvadj.mod
//Same as 4 except we put in Hicks-Neutral technical progress using g_h and g_l


// model's variables  
// we matched the names to Dale's notes. variable names come first, then superscripts (if any), then subsricpts.

var lambda_c, lambda_n, lambda_yh, lambda_yl, lambda_ke, lambda_ks,
 lambda_deh, lambda_del, lambda_dsh, lambda_dsl, lambda_je, lambda_js,
 c, c_h, c_l, n_h, n_l, y_h, y_l,
 k_eh, d_eh, j_eh, k_el, d_el, j_el, k_sh, d_sh, j_sh, k_sl, d_sl, j_sl,
 i_eh, i_el, i_sh, i_sl,  j_e, j_s,
 a_h, a_l, z_eh, z_el, z_sh, z_sl, g_h, g_l,
 jB1_eh, jB1_el,  jB1_sh, jB1_sl,
 lambdaF1_deh, lambdaF1_del, lambdaF1_dsh, lambdaF1_dsl,
 j_eh_cp, j_el_cp, j_sh_cp, j_sl_cp,
 j_cp, c_cp, y_cp, j_cp_share, c_cp_share, p_ce, r, y_qa_h, y_qa_l,
 y_cp_share_h, y_cp_share_l rpjepc;

// innovations to shock processes
varexo eps_ah, eps_al, eps_zeh, eps_zel, eps_zsh, eps_zsl, eps_gh, eps_gl;

// model's parameters
 parameters beta, eta, gamma, phi_ch, g_c, sigma_c, nbar;
 parameters alpha_nh, alpha_eh, alpha_sh, alpha_nl, alpha_el, alpha_sl;
 parameters g_e, g_s;
 parameters sigma_h, sigma_l, Asigma_e, Asigma_s;
 parameters delta_eh, delta_el, delta_sh, delta_sl;
 parameters phi_eh, phi_sh;
 parameters omega_eh, omega_sh, omega_el, omega_sl;  
 parameters nu_eh, nu_el, nu_sh, nu_sl;
 parameters rho_ah, rho_al, rho_zeh, rho_zel, rho_zsh, rho_zsl, rho_gh, rho_gl;


// the file param_reset is overwritten with the
// relevant parameter_reset file by the call program
param_simple
param_reset
param_composite


model;
/////////////////////////////////////////////////////////////////


// 1
lambda_c = 1/(1-eta) * ((c-eta*c(-1))/(1-eta))^(-gamma) 
         - eta/(1-eta)*((c(1)-eta*c)/(1-eta))^(-gamma);
// 2
c^((sigma_c-1)/sigma_c) = g_c * ( phi_ch *(c_h/phi_ch)^((sigma_c-1)/sigma_c)
           +(1-phi_ch)*(c_l/(1-phi_ch))^((sigma_c-1)/sigma_c));

// 3 - 4
c_h = (g_c)^(sigma_c-1) * phi_ch * c*(lambda_c/lambda_yh)^sigma_c ;
 c_l = (g_c)^(sigma_c-1) *(1-phi_ch) * c*(lambda_c/lambda_yl)^sigma_c ;

// 5
nbar = n_h + n_l;

// 6 - 7
n_h = (g_h*a_h)^(sigma_h-1)* alpha_nh*y_h*(lambda_yh/lambda_n)^(sigma_h);
n_l = (g_l*a_l)^(sigma_l-1)* alpha_nl*y_l*(lambda_yl/lambda_n)^(sigma_l);

// 8 - 9
y_h = c_h + i_eh + i_sh;
y_l = c_l + i_el + i_sl;

// 10 - 11



y_h^((sigma_h-1)/sigma_h) = g_h^((sigma_h-1)/sigma_h) * ( alpha_nh * (a_h*n_h/alpha_nh)^((sigma_h-1)/sigma_h)
             + alpha_eh * (k_eh/alpha_eh)^((sigma_h-1)/sigma_h)
             + alpha_sh * (k_sh/alpha_sh)^((sigma_h-1)/sigma_h) );


y_l^((sigma_l-1)/sigma_l) = g_l^((sigma_l-1)/sigma_l) * ( alpha_nl * (a_l*n_l/alpha_nl)^((sigma_l-1)/sigma_l)
             + alpha_el * (k_el/alpha_el)^((sigma_l-1)/sigma_l)
             + alpha_sl * (k_sl/alpha_sl)^((sigma_l-1)/sigma_l) );

// 12 - 13
k_eh + k_el =  d_eh*(1 - omega_eh/2*(k_eh/d_eh-1)^2) + 
              d_el*(1 - omega_el/2*(k_el/d_el-1)^2);

k_sh + k_sl =  d_sh *(1 - omega_sh/2*(k_sh/d_sh-1)^2) + 
              d_sl *(1 - omega_sl/2*(k_sl/d_sl-1)^2);

// 14 - 17
d_eh = (1-delta_eh)*k_eh(-1) + z_eh(-1)*j_eh(-1) * (1 -nu_eh/2*(j_eh(-1)*z_eh(-1)/jB1_eh(-1)/z_eh(-2) -1 )^2 ); 
d_el = (1-delta_el)*k_el(-1) + z_el(-1)*j_el(-1) * (1 -nu_el/2*(j_el(-1)*z_el(-1)/jB1_el(-1)/z_el(-2) -1 )^2 );
d_sh = (1-delta_sh)*k_sh(-1) + z_sh(-1)*j_sh(-1) * (1 -nu_sh/2*(j_sh(-1)*z_sh(-1)/jB1_sh(-1)/z_sh(-2) -1 )^2 );
d_sl = (1-delta_sl)*k_sl(-1) + z_sl(-1)*j_sl(-1) * (1 -nu_sl/2*(j_sl(-1)*z_sl(-1)/jB1_sl(-1)/z_sl(-2) -1 )^2 ); 

// 18 - 21
jB1_eh = j_eh(-1);
jB1_el = j_el(-1);
jB1_sh = j_sh(-1);
jB1_sl = j_sl(-1);

// 22 - 25
lambda_ke*(1+omega_eh*(k_eh/d_eh - 1) ) = lambda_deh(1)*beta* (1- delta_eh)
          + lambda_yh *g_h^((sigma_h-1)/sigma_h)*(alpha_eh*y_h/k_eh)^(1/sigma_h);

lambda_ks*(1+omega_sh*(k_sh/d_sh - 1) ) = lambda_dsh(1)*beta* (1- delta_sh)
          + lambda_yh *g_h^((sigma_h-1)/sigma_h)*(alpha_sh*y_h/k_sh)^(1/sigma_h);

lambda_ke*(1+omega_el*(k_el/d_el - 1) ) = lambda_del(1)*beta* (1- delta_el)
          + lambda_yl *g_l^((sigma_l-1)/sigma_l)*(alpha_el*y_l/k_el)^(1/sigma_l);

lambda_ks*(1+omega_sl*(k_sl/d_sl - 1) ) = lambda_dsl(1)*beta* (1- delta_sl)
          + lambda_yl *g_l^((sigma_l-1)/sigma_l)*(alpha_sl*y_l/k_sl)^(1/sigma_l);
    
// 26 - 29
lambda_ke * (1 + omega_eh/2 * (( k_eh/d_eh )^2-1)) = lambda_deh;
lambda_ks * (1 + omega_sh/2 * (( k_sh/d_sh )^2-1)) = lambda_dsh;
lambda_ke * (1 + omega_el/2 * (( k_el/d_el )^2-1)) = lambda_del;
lambda_ks * (1 + omega_sl/2 * (( k_sl/d_sl )^2-1)) = lambda_dsl;

// 30
j_e = j_eh + j_el;

// 31
j_e^((Asigma_e-1)/Asigma_e) = g_e^((Asigma_e-1)/Asigma_e) * ( phi_eh *(i_eh/phi_eh)^((Asigma_e-1)/Asigma_e)
                + (1-phi_eh)*(i_el/(1-phi_eh))^((Asigma_e-1)/Asigma_e) );

// 32-33 
i_eh =g_e^(Asigma_e-1)*phi_eh*j_e*( lambda_je/lambda_yh)^Asigma_e;
i_el =g_e^(Asigma_e-1)*(1-phi_eh)*j_e*( lambda_je/lambda_yl)^Asigma_e;


// 34
j_s = j_sh + j_sl;

// 35
j_s^((Asigma_s-1)/Asigma_s) = g_s^((Asigma_s-1)/Asigma_s) * (phi_sh *(i_sh/phi_sh)^((Asigma_s-1)/Asigma_s)
                + (1-phi_sh)*(i_sl/(1-phi_sh))^((Asigma_s-1)/Asigma_s));
// 36 - 37
i_sh =g_s^(Asigma_s-1)*phi_sh*j_s*( lambda_js/lambda_yh)^Asigma_s;
i_sl =g_s^(Asigma_s-1)*(1-phi_sh)*j_s*( lambda_js/lambda_yl)^Asigma_s;


// 38 - 41
lambda_je = lambda_deh(1)* beta * z_eh *( 1- nu_eh/2* (j_eh*z_eh/j_eh(-1)/z_eh(-1)-1)^2 - j_eh*nu_eh*(j_eh*z_eh/j_eh(-1)/z_eh(-1)-1)*z_eh/j_eh(-1)/z_eh(-1) )
            + lambdaF1_deh(1)*beta^2 * z_eh(1)*j_eh(1)*nu_eh*(j_eh(1)*z_eh(1)/j_eh/z_eh-1)*j_eh(1)*z_eh(1)/j_eh^2/z_eh;

lambda_je = lambda_del(1)* beta * z_el *( 1- nu_el/2* (j_el*z_el/j_el(-1)/z_el(-1)-1)^2 - j_el*nu_el*(j_el*z_el/j_el(-1)/z_el(-1)-1)*z_el/j_el(-1)/z_el(-1) )
            + lambdaF1_del(1)*beta^2 * z_el(1)*j_el(1)*nu_el*(j_el(1)*z_el(1)/j_el/z_el-1)*j_el(1)*z_eh(1)/j_el^2/z_el;

lambda_js = lambda_dsh(1)* beta * z_sh *( 1- nu_sh/2* (j_sh/j_sh(-1)-1)^2 - j_sh*nu_sh*(j_sh/j_sh(-1)-1)/j_sh(-1) )
            + lambdaF1_dsh(1)*beta^2 * z_sh(1)*j_sh(1)*nu_sh*(j_sh(1)/j_sh-1)*j_sh(1)/j_sh^2;

lambda_js = lambda_dsl(1)* beta * z_sl *( 1- nu_sl/2* (j_sl/j_sl(-1)-1)^2 - j_sl*nu_sl*(j_sl/j_sl(-1)-1)/j_sl(-1) )
            + lambdaF1_dsl(1)*beta^2 * z_sl(1)*j_sl(1)*nu_sl*(j_sl(1)/j_sl-1)*j_sl(1)/j_sl^2;



// 42 - 45
lambdaF1_deh = lambda_deh(1);
lambdaF1_del = lambda_del(1);
lambdaF1_dsh = lambda_dsh(1);
lambdaF1_dsl = lambda_dsl(1);

// 46- 49
log(a_h) = rho_ah * log(a_h(-1)) + eps_ah;
log(a_l) = rho_al * log(a_l(-1)) + eps_al;
log(z_eh) = rho_zeh * log(z_eh(-1)) + eps_zeh;
log(z_el) = rho_zel * log(z_el(-1)) + eps_zel;
log(z_sh) = rho_zsh * log(z_sh(-1)) + eps_zsh;
log(z_sl) = rho_zsl * log(z_sl(-1)) + eps_zsl;
log(g_h) = rho_gh * log(g_h(-1)) + eps_gh;
log(g_l) = rho_gl * log(g_l(-1)) + eps_gl;



// Aggregate variables at constant  prices

// Quality-adjusted investment

  j_eh_cp = z_eh(1)*j_eh;

 j_el_cp = z_el(1)*j_el;

  j_sh_cp = z_sh(1)*j_sh;

  j_sl_cp = z_sl(1)*j_sl;

j_cp = z_eh(1)*j_eh + z_el(1)*j_el+ z_sh(1)*j_sh + z_sl(1)*j_sl;

// Consumption
c_cp  =  c_h+c_l;

// Total output
y_cp = j_cp+c_cp;

// Quality-adjusted sectoral output (only goes through with full specialization)
y_qa_h = z_eh(1)*j_eh + z_el(1)*j_el;
y_qa_l = z_sh(1)*j_sh + z_sl(1)*j_sl+c_l;

// Quality-adjusted sectoral output (only goes through with full specialization)
y_cp_share_h = y_qa_h/y_cp;
y_cp_share_l = y_qa_l/y_cp;

// Constant price share of investment
j_cp_share = j_cp/y_cp;

// Constant price share of consumption
c_cp_share = c_cp/y_cp;

//Consumption price of equipment
p_ce = lambda_je/(z_eh(1)*lambda_c);

// real interest rate
1/r = beta*lambda_c(1)/lambda_c;

// P_je/P_c
rpjepc = (lambda_je/lambda_c)/z_eh;

end;

////////////////////////////////////////////////////////////////////
// STEADY STATES
////////////////////////////////////////////////////////////////////
initval;

c = (1-s)*y;
c_h = phi_ch*c;
c_l = (1-phi_ch)*c;
y_h = y_h2y*y;
y_l = y_l2y*y;
n_h = n_h2y_h*y_h;
n_l = n_l2y_l*y_l;
k_eh = k_eh2y_h*y_h;
d_eh = k_eh;
j_eh = delta_eh*k_eh;
k_el = k_el2y_l*y_l;
d_el = k_el;
j_el = delta_el*k_el;
k_sh = k_sh2y_h*y_h;
d_sh = k_sh;
j_sh = delta_sh*k_sh;
k_sl = k_sl2y_l*y_l;
d_sl = k_sl;
j_sl = delta_sl*k_sl;
j_e = j_el + j_eh;
j_s = j_sl + j_sh;
i_eh = g_e^(-1)*phi_eh*j_e;
i_el = g_e^(-1)*(1-phi_eh)*j_e;
i_sh = g_s^(-1)*phi_sh*j_s;
i_sl = g_s^(-1)*(1-phi_sh)*j_s;
a_h = 1;
a_l = 1;
z_eh = 1;
z_el = 1;
z_sh = 1;
z_sl = 1;
g_h = 1;
g_l = 1;
jB1_eh = j_eh;
jB1_el = j_el;
jB1_sh = j_sh;
jB1_sl = j_sl;
lambda_c = c^(-gamma);
lambda_n = lambda_c;
lambda_yh = lambda_c;
lambda_yl = lambda_c;
lambda_je = lambda_c/g_e;
lambda_js = lambda_c/g_s;
rpjepc = lambda_je/lambda_c;
lambda_ke = lambda_je/beta;
lambda_ks = lambda_js/beta;
lambda_deh = lambda_je/beta;
lambda_del = lambda_je/beta;
lambda_dsh = lambda_js/beta;
lambda_dsl = lambda_js/beta;
lambdaF1_deh = lambda_deh;
lambdaF1_del = lambda_del;
lambdaF1_dsh = lambda_dsh;
lambdaF1_dsl = lambda_dsl;
 j_eh_cp = z_eh*j_eh;
 j_el_cp = z_el*j_el;
 j_sh_cp = z_sh*j_sh;
 j_sl_cp = z_sl*j_sl;
j_cp = z_eh*j_eh + z_el*j_el+ z_sh*j_sh + z_sl*j_sl;
c_cp = c_h+c_l;
y_cp = j_cp+c_cp;
j_cp_share = j_cp/y_cp;
c_cp_share = c_cp/y_cp;
p_ce = lambda_je/(z_eh(1)*lambda_c);
r = 1/beta;
y_qa_h = j_eh +j_el;
y_qa_l = j_sh + j_sl+c_l;
y_cp_share_h = y_qa_h/y_cp;
y_cp_share_l = y_qa_l/y_cp;
end;


shocks;
  var eps_ah; stderr 0.01;
  var eps_al; stderr 0.01;
  var eps_zeh; stderr 0.01;
  var eps_zel; stderr 0.01;
  var eps_zsh; stderr 0.01;
  var eps_zsl; stderr 0.01;
  var eps_gh; stderr 0.01;
  var eps_gl; stderr 0.01;
end;



stoch_simul(order=1,nocorr,nomoments,irf=0);
