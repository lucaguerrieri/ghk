function e = fbeta(p2,p,p1,perc)
% must restrict p2 such that a>0 and b>0 ....
  a = (1-p1)*p1^2/p2^2 - p1;
  b = a*(1/p1 - 1);

  e = p - pbeta(perc,a,b);