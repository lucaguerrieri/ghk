% solves iteratively ax+bxc=d

function x=sylvester3a(x0,a,b,c,d)
a_1 = inv(a);
b = a_1*b;
d = a_1*d;
e = 1;
iter = 1;
while e > 1e-8 & iter < 500
  x = d-b*x0*c;
  e = max(max(abs(x-x0)));
  x0 = x;
  iter = iter + 1;
end
if iter == 500
  fprintf('Syvester3a : Only accuracy of %10.8f is achieved after 500 iterations. \n',e)
end