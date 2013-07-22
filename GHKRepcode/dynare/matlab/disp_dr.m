% Copyright (C) 2001 Michel Juillard
%
function disp_dr(dr,order,var_list)
global lgy_ lgx_ endo_nbr exo_nbr ykmin_
  
nx =size(dr.ghx,2);
nu =size(dr.ghu,2);
k = find(dr.kstate(:,2) <= ykmin_+1);
klag = dr.kstate(k,[1 2]);

k1 = dr.order_var;

nvar = size(var_list,1);
if nvar == 0
  nvar = length(k1);
  ivar = [1:nvar];
else
  ivar=zeros(nvar,1);
  for i=1:nvar
    i_tmp = strmatch(var_list(i,:),lgy_(k1,:),'exact');
    if isempty(i_tmp)
      disp(var_list(i,:));
      error (['One of the variable specified does not exist']) ;
    else
      ivar(i) = i_tmp;
    end
  end
end
disp('POLICY AND TRANSITION FUNCTIONS')
% variable names
str = '                        ';
for i=1:nvar
  str = [str sprintf('%16s',lgy_(k1(ivar(i)),:))];
end
disp(str);
%
% constant
%
str = 'Constant            ';
flag = 0;
for i=1:nvar
  x = dr.ys(k1(ivar(i)));
  if order > 1
    x = x + dr.ghs2(ivar(i))/2;
  end
  if abs(x) > 1e-6
    flag = 1;
    str = [str sprintf('%16.6f',x)];
  else
    str = [str '               0'];
  end
end
if flag
  disp(str)
end
if order > 1
  str = '(correction)        ';
  flag = 0;
  for i=1:nvar
    x = dr.ghs2(ivar(i))/2;
    if abs(x) > 1e-6
      flag = 1;
      str = [str sprintf('%16.6f',x)];
    else
      str = [str '               0'];
    end
  end
  if flag
    disp(str)
  end
end
%
% ghx
%
for k=1:nx
  flag = 0;
  str1 = sprintf('%s(%d)',lgy_(k1(klag(k,1)),:),klag(k,2)-ykmin_-2);
  str = sprintf('%-20s',str1);
  for i=1:nvar
    x = dr.ghx(ivar(i),k);
    if abs(x) > 1e-6
      flag = 1;
      str = [str sprintf('%16.6f',x)];
    else
      str = [str '               0'];
    end
  end
  if flag
    disp(str)
  end
end
%
% ghu
%
for k=1:nu
  flag = 0;
  str = sprintf('%-20s',lgx_(k,:));
  for i=1:nvar
    x = dr.ghu(ivar(i),k);
    if abs(x) > 1e-6
      flag = 1;
      str = [str sprintf('%16.6f',x)];
    else
      str = [str '               0'];
    end
  end
  if flag
    disp(str)
  end
end

if order > 1
  % ghxx
  for k = 1:nx
    for j = 1:k
      flag = 0;
      str1 = sprintf('%s(%d),%s(%d)',lgy_(k1(klag(k,1)),:),klag(k,2)-ykmin_-2, ...
		     lgy_(k1(klag(j,1)),:),klag(j,2)-ykmin_-2);
      str = sprintf('%-20s',str1);
      for i=1:nvar
	if k == j
	  x = dr.ghxx(ivar(i),(k-1)*nx+j)/2;
	else
	  x = dr.ghxx(ivar(i),(k-1)*nx+j);
	end
	if abs(x) > 1e-6
	  flag = 1;
	  str = [str sprintf('%16.6f',x)];
	else
	  str = [str '               0'];
	end
      end
      if flag
	disp(str)
      end
    end
  end
  %
  % ghuu
  %
  for k = 1:nu
    for j = 1:k
      flag = 0;
      str = sprintf('%-20s',[lgx_(k,:) ',' lgx_(j,:)] );
      for i=1:nvar
	if k == j
	  x = dr.ghuu(ivar(i),(k-1)*nu+j)/2;
	else
	  x = dr.ghuu(ivar(i),(k-1)*nu+j);
	end
	if abs(x) > 1e-6
	  flag = 1;
	  str = [str sprintf('%16.6f',x)];
	else
	  str = [str '               0'];
	end
      end
      if flag
	disp(str)
      end
    end
  end
  %
  % ghxu
  %
  for k = 1:nx
    for j = 1:nu
      flag = 0;
      str1 = sprintf('%s(%d),%s',lgy_(k1(klag(k,1)),:),klag(k,2)-ykmin_-2, ...
		     lgx_(j,:));
      str = sprintf('%-20s',str1);
      for i=1:nvar
	x = dr.ghxu(ivar(i),(k-1)*nu+j);
	if abs(x) > 1e-6
	  flag = 1;
	  str = [str sprintf('%16.6f',x)];
	else
	  str = [str '               0'];
	end
      end
      if flag
	disp(str)
      end
    end
  end
end

% $$$   dr.ghx
% $$$   dr.ghu
% $$$   dr.ghxx
% $$$   dr.ghuu
% $$$   dr.ghxu

% 01/08/2001 MJ  added test for order in printing quadratic terms
% 02/21/2001 MJ pass all variable names through deblank()
% 02/21/2001 MJ changed from f to g format to write numbers
% 10/09/2002 MJ corrected error on constant whith subset of variables 



