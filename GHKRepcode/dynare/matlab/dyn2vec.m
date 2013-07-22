% Copyright (C) 2001 Michel Juillard
%
function [z,zss]=dyn2vec(s1,s2);

  global lgy_ lgx_ y_ ys_ ykmin_ xkmin_ ex_ dsmpl_

  if dsmpl_ == 0
    k = [1:size(y_,2)];
  else
    k = [ykmin_+dsmpl_(1):ykmin_+dsmpl_(2)];
  end

  if nargin == 0
    if nargout > 0
      t = ['DYNARE dyn2vec error: the function doesn''t return values when' ...
	   ' used without input argument'];
      error(t);
    end
    for i=1:size(y_,1)
      assignin('base',deblank(lgy_(i,:)),y_(i,k)');
    end
    return
  else
    j = strmatch(s1,lgy_,'exact'); 
    if ~ isempty(j)
      z = y_(j,k)';
    else
      j = strmatch(s1,lgx_,'exact');
      if ~ isempty(j)
	if dsmpl_ == 0
	  z = ex_(:,j);
	else
	  z = ex_(xkmin_+dsmpl_(1):xkmin_+dsmpl_(2));
	end
      else
	t = ['DYNARE dyn2vec error: variable ' deblank(s1(i,:)) ' doesn''t' ...
	     ' exist.'] ;
	error (t) ;
      end
    end
  end

  if nargout == 0
    if nargin == 1
      assignin('base',s1,z);
    elseif nargin == 2
      assignin('base',s2,z);
    end
  else
    zss=ys_(j);
  end
  
% 02/23/01 MJ redone, incorporating FC's improvements
% 08/24/01 MJ replaced globlize by internal assignin
% 08/24/01 MJ added 'exact' to strmatch (thanks to David Vavra)
% 01/31/03 MJ added provision for alternative name of variable



