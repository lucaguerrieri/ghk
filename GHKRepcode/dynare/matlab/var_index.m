% var_index(vname) returns the index of variable 'vname'
% in its respective variable list (endogeneous, exogeneous)
% input:
%        vname     a string of character
% output:
%        i         a positive integer
% an error message is given if vname isn't in the state vector
function k=var_index(vname)
  global lgy_
  
  k = strmatch(vname,lgy_,'exact');
