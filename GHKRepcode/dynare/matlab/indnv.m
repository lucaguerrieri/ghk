% Copyright (C) 2001 Michel Juillard
%
function a=indnv(x,y)

a = zeros(size(x)) ;

for i = 1:size(x,1)
  j = find( x(i) == y );
  if isempty(j)
    a(i) = NaN;
  else
    a(i) = j;
  end
end



