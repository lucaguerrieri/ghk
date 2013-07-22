function dist = findss(p_nss)

global gamma_i gamma_c alpha_m alpha_n delta betap 


p_i2p_n= (1/gamma_i/p_nss)^gamma_i*(1/(1-gamma_i))^(1-gamma_i);
p_c2p_n= (1/gamma_c/p_nss)^gamma_c*(1/(1-gamma_c))^(1-gamma_c);

p_i2p_m=p_i2p_n*p_nss;
p_c2p_m=p_c2p_n*p_nss;

lss = 1;
r2p_i = (1-betap*(1-delta))/betap;


k_m2y_m = alpha_m/(r2p_i)/p_i2p_m;
k_n2y_n = alpha_n/(r2p_i)/p_i2p_n;

mymat = [1, 0, (-gamma_i*p_i2p_m+gamma_c*p_i2p_m)
         0, 1, -(1-gamma_i)*p_i2p_n+(1-gamma_c)*p_i2p_n
         delta*k_m2y_m, delta*k_n2y_n, -1];


myvars = mymat^-1*[gamma_c; (1-gamma_c)/p_nss; 0];

y_m2y = myvars(1);
y_n2y = myvars(2);


y_m2y_n = y_m2y/y_n2y;
l_m2y_m = k_m2y_m^(-alpha_m/(1-alpha_m));
l_n2y_n = k_n2y_n^(-alpha_n/(1-alpha_n));

l_m2l_n=l_m2y_m/l_n2y_n*y_m2y_n;


l_nss = lss*(1+l_m2l_n)^(-1);

k_nss = (l_n2y_n)^(-1/alpha_n)*l_nss;
y_nss = k_nss^(alpha_n)*l_nss^(1-alpha_n);

w_ss = (1-alpha_m)/l_m2y_m;

w_nss = (1-alpha_n)*p_nss*y_nss/l_nss;

dist = w_ss-w_nss;