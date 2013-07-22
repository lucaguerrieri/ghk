% Copyright (C) 2001 Michel Juillard
%
function dynasave(s,var_list)
% DYNASAVE :	DYNASAVE ( [ 'filename' ] )	
%		This optional command saves the simulation results
%		in a .MAT file.

  global endo_nbr lgy_ y_

  n = size(var_list,1);
  if n == 0
    n = endo_nbr;
    ivar = [1:n]';
    var_list = lgy_;
  else
    ivar=zeros(n,1);
    for i=1:n
      i_tmp = strmatch(var_list(i,:),lgy_,'exact');
      if isempty(i_tmp)
	error (['One of the specified variables does not exist']) ;
      else
	ivar(i) = i_tmp;
      end
    end
  end


%  dyn2vec(var_list(1),var_list(1));
eval([var_list(1) '=y_(ivar(1),:)'';'])
eval(['save ' s ' ' var_list(1) ' -mat'])
  for i = 2:n
%    dyn2vec(var_list(i),var_list(i));
    eval([var_list(i) '=y_(ivar(i),:)'';'])
    eval(['save ' s ' ' var_list(i) ' -append -mat'])
  end


