% Copyright (C) 2001 Michel Juillard
%
function mcompare(s1,s2)
% MCOMPARE :	MCOMPARE ( [ 'file1' ; 'file2' ] , [ 'var1' ; 'var2' ...] )	
%		This optional command plots the relative differences between
%		two different simulations for a list of variables. One plot 
%		is drawn for each variable. The trajectories must have been
%		previously saved by the instruction DYNASAVE. The simulation
%		in file1 serves as the base simulation and the ploted quantity
%		is equal to the difference between the two simulation reported
%		to the first one. If, for a given variable, zero is one of the
%		value of the base simulation, the absolute difference is ploted
%		instead of the relative one.

global dsmpl_ iter_
global nvx nvy x y lag1

ftest(s1,s2) ;

ix = [1-lag1(1):size(x,2)-lag1(1)]' ;
i = [lag1(1):size(ix,1)-lag1(2)+1]' ;

if size(dsmpl_,1) == 1
	error(['DSAMPLE not specified.']) ;
end

if dsmpl_(3) > 0
	if dsmpl_(3) == 2
		if dsmpl_(1)<0 | dsmpl_(2)>size(x,2)-lag1(2)
			error ('Wrong sample.') ;
		end
		i = [dsmpl_(1)+lag1(1):dsmpl_(2)+lag1(1)]' ;
	elseif dsmpl_(3) == 1
		if dsmpl_(1)>size(x,2)-lag1(2)
			error ('Wrong sample.') ;
		end
		i = [lag1(1):dsmpl_(1)+lag1(1)]' ;
	end
end

for k = 1:size(x,1)
	figure ;
	x1 = x(k,i) ;
	y1 = y(k,i) ;
	if nnz(x1) < length(x1)
		plot(ix(i),(y1-x1)) ;
	else
		plot(ix(i),(y1-x1)./x1) ;
	end
	xlabel(['Periods']) ;
	title(['Variable ' s2(k)]) ;
end

return ;
