function e = figamm(p2,p,p1,perc)
  e = p - pgamma(1/perc,p2/2,2/p1);