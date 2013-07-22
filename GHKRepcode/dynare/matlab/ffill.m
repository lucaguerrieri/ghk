% Copyright (C) 2001 Michel Juillard
%
function [a,b] = ffill(x,ixc,y)

xc = size(x,1) ;

if isempty(y)
  b = [ixc; 0];
  a = [x zeros(size(x,1),1)];
else
  yc = size(y,1) ;
  b = [ixc;yc] ;

  if xc > yc
    a = [x [y;zeros(xc-yc,size(y,2))]] ;
  elseif yc > xc
    a = [[x;zeros(yc-xc,size(x,2))] y] ;
  else
    a = [x y] ;
  end

end

% 2001/09/1 MJ corrected for absent lags