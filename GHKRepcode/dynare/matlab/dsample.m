% Copyright (C) 2001 Michel Juillard
%
function dsample(s1,s2)
% DSAMPLE :	DSAMPLE(d1,d2)
%		This optional command permits to reduce the number of
%		periods considered in following output commands. 
%               If only one argument is 
%		provided, output is from period 1 to the period 
%		specified in the DSAMPLE command. If two arguments are
%		present output is done for the interval between the 
%		two periods.
%               DSAMPLE without arguments reset the sample to the one 
%               specified by PERIODS 
		
global dsmpl_ iter_

dsmpl_ = zeros(2,1) ;

if s1 > iter_ | s2 > iter_
  t = ['DYNARE dsample error: one of the arguments is larger than the one' ...
       ' specified in PERIODS'];
  error(t);
end

if nargin == 0
	dsmpl_(1) = 1 ;
	dsmpl_(2) = iter_ ;
elseif nargin == 1
	dsmpl_(1) = 1 ;
	dsmpl_(2) = s1 ;
else
	dsmpl_(1) = s1 ;
	dsmpl_(2) = s2 ;
end

% 02/23/01 MJ added error checking