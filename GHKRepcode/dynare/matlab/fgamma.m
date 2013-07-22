function e = fgamma(p2,p,p1,perc)
  b = p2^2/p1;
  a = p1/b;
  e = p - pgamma(perc,a,b);