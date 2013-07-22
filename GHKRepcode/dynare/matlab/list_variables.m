function list_variables()
  global lgy_ iy_ ykmin_ ykmax_ endo_nbr
  
  ipred = find(sum(iy_(1:ykmin_,:),1))';
  ifwrd = find(sum(iy_(ykmin_+2:end,:),1))';
  iboth = intersect(ipred,ifwrd);
  istatic = setdiff([1:endo_nbr]',union(ipred,ifwrd));
  disp('Purely predetermined variables')
  ipred1 = setdiff(ipred,iboth);
  if isempty(ipred1)
    disp('    none')
  else
    disp(lgy_(ipred1,:))
  end
  disp('Purely forward-looking variables')
  ifwrd1 = setdiff(ifwrd,iboth);
  if isempty(ifwrd1)
    disp('    none')
  else
    disp(lgy_(ifwrd1,:))
  end
  disp('Variables that are both forward-looking and predetermined')
  if isempty(iboth)
    disp('    none')
  else
    disp(lgy_(iboth,:))
  end
  disp('Static variables')
  if isempty(istatic)
    disp('    none')
  else
    disp(lgy_(istatic,:))
  end
  
