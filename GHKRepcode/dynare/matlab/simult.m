% Copyright (C) 2001 Michel Juillard
%

function y_=simult(ys, dr)
global jacobia_ iy_ ykmin_ ykmax_ gstep_ exo_nbr endo_nbr
global ex_ valf_ it_ exe_ xkmin_ xkmax_ ys_
global fname_ means_ Sigma_e_ lgy_ options_

order = options_.order;
seed = options_.simul_seed;
iter_ = options_.periods;

it_ = ykmin_ + 1 ;

% eliminate shocks with 0 variance
i_exo_var = setdiff([1:exo_nbr],find(diag(Sigma_e_) == 0));
nxs = length(i_exo_var);
ex_ = zeros(xkmin_+xkmax_+iter_,exo_nbr);
chol_S = chol(Sigma_e_(i_exo_var,i_exo_var));

if isempty(seed)
  randn('state',sum(100*clock));
else
  randn('state',seed);
end
if ~isempty(Sigma_e_)
  ex_(:,i_exo_var) = randn(xkmin_+xkmax_+iter_,nxs)*chol_S;
end
y_ = simult_(ys,dr,ex_,order);



% 02/20/01 MJ replaced ys by dr.ys
% 02/22/01 MJ removed commented out lines
%             removed useless temps
%             stderr_ replaced by Sigma_e_
% 02/28/01 MJ changed expression for Sigma_e_
% 02/18/03 MJ added ys in the calling sequence for arbitrary initial values
%             suppressed useless calling parameter istoch
% 05/10/03 MJ removed repmat() in call to simult_() for lag > 1
% 05/29/03 MJ test for 0 variances
