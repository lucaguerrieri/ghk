function y=p2toperc(x)
  
  n = size(x,1);
  y = zeros(n,1);
  
  for i=1:n
%    if x(i,6) == 0  have to wait until the 2nd part works
    if 1
      y(i) = x(i,3);
    else
      if x(i,3) < x(i,2)
	p = 0.05;
      else
	p = 0.95;
      end
    
      if x(i,1) == 1
	y(i) = fsolve(@fbeta,1,p,x(i,2),x(i,3));
      elseif x(i,1) == 2
	y(i) = fsolve(@fgamma,1,p,x(i,2),x(i,3));
      elseif x(i,1) == 3
	y(i) = (x(i,3)-x(i,2))/qnorm(p,0,1);
      elseif x(i,1) == 4
	y(i) = fsolve(@figamm,1,p,x(i,2),x(i,3));
      elseif x(i,1) == 5
%	y(i) = fsolve(@fgamma,1,p,x(i,2));
      end
    end
  end