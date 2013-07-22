% Copyright (C) 2001 Michel Juillard
%
function y=ff2_(x)
global np_ endo_nbr nf_ ex_ it_ ykmin_ xkmin_ exo_nbr jp_ jf_ kf_ exe_
y=zeros(np_+endo_nbr+nf_,1);
y(1:np_)=x(1:np_);
ex_(it_-ykmin_+xkmin_,:)=x(np_+1:np_+exo_nbr)';
y(np_+jp_)=h1_(y(1:np_),ex_(it_-ykmin_+xkmin_,:)');
y(np_+jf_)=g1_(y(1:np_),ex_(it_-ykmin_+xkmin_,:)');
y(kf_)=g1_(y(np_+jp_),exe_);
y=ff_(y);



