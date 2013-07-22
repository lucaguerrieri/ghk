% Copyright (C) 2001 Michel Juillard
%
% function used for dr_algo == 1
function ghs2=dr2(ys,dr)
  global fname_
  
  dr.ys = ys;
  fh = str2func([fname_ '_fff']);
  dr.fbias = 2*feval(fh,dr.ys);
  dr=dr1(dr,0);
  ghs2 = dr.ghs2;
