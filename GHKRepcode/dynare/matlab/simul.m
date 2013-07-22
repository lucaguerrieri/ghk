% Copyright (C) 2001 Michel Juillard
%
function simul(dr)

global scalv_ ct_ endval_ endo_nbr ex_ exe_
global iy_ ykmin_ ykmax_ xkmin_ xkmax_
global valf_ ys_ ys0_ y_ lgy_ iter_ options_

if size(iy_,2)-nnz(iy_(ykmin_+1,:)) > 0
  mess = ['DYNARE: error in model specification : variable ' lgy_(find(iy_(ykmin_+1,:)==0),:)] ;
  mess = [mess ' doesn''t appear as current variable.'] ; 
  error (mess) ;
end

options_ = set_default_option(options_,'simul_algo',0);
if ~isfield(options_,'periods') & ~isempty(iter_)
  options_.periods = iter_
end
options_ = set_default_option(options_,'periods',0);
if options_.periods == 0
  error('SIMUL: number of periods for the simulation isn''t specified')
end
iter_ = options_.periods;
if options_.simul_algo == 0
  if ~ valf_
    make_y_;
    make_ex_;
  end

  if isempty(scalv_) | scalv_ == 0
    scalv_ = ys_ ;
  end

  scalv_= 1 ;

  if ykmin_ ==1 & ykmax_ <= 1
    sim1 ;
  else
    simk ;
  end
else
  set_default_option('replic',1);
  set_default_option('simul_seed',1);
  if isfield(dr,'ghxx')
    set_default_option('order',2);
  else
    set_defaut_option('order',1);
  end
  y_=simult(ys_,dr,options_);
end

dyn2vec;

% 6/18/01 MJ added dyn2vec if 40 variables or less
% 01/16/03 MJ use dyn2vec whatever the number of variables
% 02/18/03 MJ added ys_ for calling simult
% 05/24/03 MJ added options_ and options_.periods









