function check_k(y_h2y,y_l2y,y, ...
                 k_eh2y_h, k_el2y_l, k_sh2y_h, k_sl2y_l)
             
y_h_ss = y_h2y*y;
y_l_ss = y_l2y*y;

k_eh = k_eh2y_h *y_h_ss;
k_el = k_el2y_l *y_l_ss;
k_sh = k_sh2y_h *y_h_ss;
k_sl = k_sl2y_l *y_l_ss;

check_vec = [k_eh,k_el,k_sh,k_sl];

if min(check_vec)<0
    error('Shares and elasticity combination imply negative capital')
end

