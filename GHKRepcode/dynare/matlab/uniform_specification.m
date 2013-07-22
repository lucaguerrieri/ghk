% Copyright (C) 2004 Michel Juillard
%
function [m,s,p1,p2] = uniform_specification(m,s,p3,p4)
    if ~(isnan(p3) | isnan(p4))
      p1 = p3;
      p2 = p4;
      m = (p3+p4)/2;
      s = (p4-p3)/(sqrt(12));
    else
      p1 = m-s*sqrt(3);
      p2 = m+s*sqrt(3);
    end
