% Copyright (C) 2001 Michel Juillard
%
function y=fx_(x)
global ys_ ykmin_ ykmax_ it_ iy_ xkmin_ ex_ fname_
y=repmat(ys_,ykmin_+ykmax_+1,1);
iyr0=find(reshape(iy_',size(iy_,1)*size(iy_,2),1));
y=y(iyr0);
ex_(it_+xkmin_-ykmin_,:)=x';
fh = str2func([fname_ '_ff']);
y=feval(fh,y);
