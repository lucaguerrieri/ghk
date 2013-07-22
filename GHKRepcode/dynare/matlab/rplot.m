% Copyright (C) 2001 Michel Juillard
%
function rplot(s1,rplottype)
% RPLOT :		RPLOT ( ['var1'; 'var2';  ...] , rplottype )	
%		This optionnal command creates the plot of the variable
%		trajectory. By default, the entire simulation period is 
%		ploted. The instruction DSAMPLE permits to reduce the 
%		number of periods in the plot.
		
global lgy_ lgx_ y_ ykmin_ ykmax_ iter_
global dsmpl_ ys_

if nargin == 1
	rplottype = 0;
end

col = ['y','c','r','g','b','w','m'] ;
ix = [1 - ykmin_:size(y_,2)-ykmin_]' ;

y = [];
for k=1:size(s1,1)
  if isempty(strmatch(s1(k,:),lgy_,'exact'))
    error (['One of the variable specified does not exist']) ;
  end

  y = [y; y_(strmatch(s1(k,:),lgy_,'exact'),:)] ;
end

if dsmpl_ == 0
        i = [ykmin_:size(y_,2)]' ;
else
	i = [dsmpl_(1)+ykmin_:dsmpl_(2)+ykmin_]' ;
end

t = ['Plot of '] ;
if rplottype == 0
	for j = 1:size(y,1)
		t = [t s1(j,:) ' '] ;
	end
        figure ;
        plot(ix(i),y(:,i)) ;
	title (t,'Interpreter','none') ;
	xlabel('Periods') ;
	if size(s1,1) > 1
	  legend(s1,0);
	end
elseif rplottype == 1
	for j = 1:size(y,1)
		figure ;
		plot(ix(i),y(j,i)) ;
		title(['Plot of ' s1(:,j)]) ;
		xlabel('Periods') ;
	end
elseif rplottype == 2
	figure ;
	nl = max(1,fix(size(y,1)/4)) ;
	nc = ceil(size(y,1)/nl) ;
	for j = 1:size(y,1)
		subplot(nl,nc,j) ;
		plot(ix(i),y(j,i)) ;
		hold on ;
		plot(ix(i),ys_(j)*ones(1,size(i,1)),'w:') ;
		xlabel('Periods') ;
		ylabel([s1(:,j)]) ;
		title(['Plot of ' s1(:,j)]) ;
	end
end

% 02/28/01 MJ replaced bseastr by MATLAB's strmatch
% 06/19/01 MJ added 'exact' to strmatch calls
% 06/25/03 MJ correction when dsmpl_ ~= 0









