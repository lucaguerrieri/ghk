% Copyright (C) 2005 Michel Juillard
%
% all local variables have complicated names in order to avoid name
% conflicts with possible user variable names

function dyn_data_01=read_variables(file_name_01,var_names_01,dyn_data_01)
  
  dyn_size_01 = size(dyn_data_01,1);
  
  if exist(file_name_01)
    dyn_instr_01 = file_name_01;
  else
    dyn_instr_01 = ['load ' file_name_01];
  end
  
  eval(dyn_instr_01);

  
  for dyn_i_01=1:size(var_names_01,1)
    dyn_tmp_01 = eval(var_names_01(dyn_i_01,:));
    if length(dyn_tmp_01) > dyn_size_01 & dyn_size_01 > 0
      error('data size is too large')
    end
    dyn_data_01(:,dyn_i_01) = dyn_tmp_01;
  end