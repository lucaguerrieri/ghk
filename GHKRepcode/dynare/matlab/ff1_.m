% Copyright (C) 2001 Michel Juillard
%
function y=ff1_(x)
global it_ exo_nbr exo_det_nbr ex_ ex_det_ xkmin_ ykmin_ fname_ M_

n1 = size(x,1) - exo_nbr - exo_det_nbr;
ex_(it_+xkmin_-ykmin_,:) = x(n1+[1:exo_nbr])';
if M_.ex_det_length > 0
  ex_det_(it_,:) = x(end-exo_det_nbr+1:end)';
end
fh = str2func([fname_ '_ff']);
y=feval(fh,x(1:n1));



