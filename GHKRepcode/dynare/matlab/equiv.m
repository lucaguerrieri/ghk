% Copyright (C) 2001 Michel Juillard
%
function y=equiv(x)
ys_=x;
taylor;
y=newt('ff20_',ys_);
