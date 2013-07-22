% solves a*x+b*x*c=d
function x=sylvester3(a,b,c,d)
  n = size(a,1);
  m = size(c,1);
  if n == 1
    x=d./(a*ones(1,m)+b*c);
    return
  end
  if m == 1
    x = (a+c*b)\d;
    return;
  end
  x=zeros(n,m);
  [u,t]=schur(c);
  [aa,bb,qq,zz]=qz(full(a),full(b),'real'); % available in Matlab version 6.0
  d=qq*d*u;
  i = 1;
  while i < m
    if t(i+1,i) == 0
      if i == 1
	c = zeros(n,1);
      else
	c = bb*(x(:,1:i-1)*t(1:i-1,i));
      end
      x(:,i)=(aa+bb*t(i,i))\(d(:,i)-c);
      i = i+1;
    else
      if i == n
	c = zeros(n,1);
	c1 = zeros(n,1);
      else
	c = bb*(x(:,1:i-1)*t(1:i-1,i));
	c1 = bb*(x(:,1:i-1)*t(1:i-1,i+1));
      end
      z = [aa+bb*t(i,i) bb*t(i+1,i); bb*t(i,i+1) aa+bb*t(i+1,i+1)]...
	  \[d(:,i)-c;d(:,i+1)-c1];
      x(:,i) = z(1:n);
      x(:,i+1) = z(n+1:end);
      i = i + 2;
    end
  end
  if i == m
    c = bb*(x(:,1:m-1)*t(1:m-1,m));
    x(:,m)=(aa+bb*t(m,m))\(d(:,m)-c);
  end
  x=zz*x*u';
  
% 01/25/03 MJ corrected bug for i==m (sign of c in x determination)
% 01/31/03 MJ added 'real' to qz call
