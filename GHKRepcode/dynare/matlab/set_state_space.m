% Copyright (C) 2001 Michel Juillard
%
function dr=set_state_space(dr)

global iy_ ykmin_ ykmax_ exo_nbr endo_nbr
global ex_ valf_ it_ exe_ xkmin_ xkmax_ ys_
global fname_ lgy_


xlen = xkmax_ + xkmin_ + 1;
klen = ykmin_ + ykmax_ + 1;

if ~ iy_(ykmin_+1,:) > 0
  error ('Error in model specification: some variables don"t appear as current') ;
end

fwrd_var = find(any(iy_(ykmin_+2:end,:),1))';
if ykmin_ > 0
  pred_var = find(any(iy_(1:ykmin_,:),1))';
  both_var = intersect(pred_var,fwrd_var);
  pred_var = setdiff(pred_var,both_var);
  fwrd_var = setdiff(fwrd_var,both_var);
  stat_var = setdiff([1:endo_nbr]',union(union(pred_var,both_var),fwrd_var));  % static variables
else
  pred_var = [];
  both_var = [];
  stat_var = setdiff([1:endo_nbr]',fwrd_var);
end
nboth = length(both_var);
npred = length(pred_var);
nfwrd = length(fwrd_var);
nstatic = length(stat_var);
order_var = [ stat_var; pred_var; both_var; fwrd_var];

% building kmask for z state vector in t+1
if ykmin_ > 0
  kmask = [];
  if ykmax_ > 0 
    kmask = [cumsum(flipud(iy_(ykmin_+2:end,order_var)),1)] ;
  end
  kmask = [kmask; flipud(cumsum(iy_(1:ykmin_,order_var),1))] ;
else
  kmask = cumsum(flipud(iy_(ykmin_+2:klen,order_var)),1) ;
end

kmask = kmask';
kmask = kmask(:);
i_kmask = find(kmask);          % index of nonzero entries in kmask
nd = size(i_kmask,1);           % size of the state vector
kmask(i_kmask) = [1:nd];

% auxiliary equations

% elements that are both in z(t+1) and z(t)
k1 = find([kmask(1:end-endo_nbr) & kmask(endo_nbr+1:end)] );
kad = [];
kae = [];
if ~isempty(k1)
  kad = kmask(k1+endo_nbr);
  kae = kmask(k1);
end

% composition of state vector
% col 1: variable;           col 2: lead/lag in z(t+1); 
% col 3: A cols for t+1 (D); col 4: A cols for t (E)
kstate = [ repmat([1:endo_nbr]',klen-1,1) kron([klen:-1:2]',ones(endo_nbr,1)) ...
	   zeros((klen-1)*endo_nbr,2)];
kiy = flipud(iy_(:,order_var))';
kiy = kiy(:);
kstate(1:ykmax_*endo_nbr,3) = kiy(1:ykmax_*endo_nbr)-endo_nbr;  
kstate(find(kstate(:,3) < 0),3) = 0;
kstate(ykmax_*endo_nbr+1:end,4) = kiy((ykmax_+1)*endo_nbr+1:end);  
% put in E only the current variables that are not already in D
kstate = kstate(i_kmask,:);

dr.order_var = order_var;
dr.nstatic = nstatic;
dr.npred = npred+nboth;
dr.kstate = kstate;
dr.kad = kad;
dr.kae = kae;
dr.nboth = nboth;
dr.nfwrd = nfwrd;
% number of forward variables in the state vector
dr.nsfwrd = sum(kstate(:,2) > ykmin_+1);
% number of predetermined variables in the state vector
dr.nspred = sum(kstate(:,2) <= ykmin_+1);

