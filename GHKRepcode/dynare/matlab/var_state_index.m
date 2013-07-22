% var_state_index(vname) returns the index of variable 'vname'
% in the extended state space representation vector of the system
% input:
%        vname     a string of character
%        order     a vector containing the composition of the state
%                  vector. Normaly, dr_.order_var 

% output:
%        i         a positive integer
% an error message is given if vname isn't in the state vector
function i=var_state_index(vname,order)
  global lgy_
  
  k = strmatch(vname,lgy_,'exact');
  i = find(k==order);