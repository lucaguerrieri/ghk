% Copyright (C) 2001 Michel Juillard
%
%	DYNARE ( 'Filename' )
%	This command runs the .MOD file (or .DYN) specified in Filename.
%	Filename could be enter with or without the .MOD extension.
function dynare(fname,varargin)

if ~ isstr(fname)
	error ('The argument in DYNARE must be a text string.') ;
end

pathfile = mfilename('fullpath');
command = [pathfile,'_m ',fname];
for i=2:nargin
  command = [command ' ' varargin{i-1}];
end
[status, result] = dos (command) ;
if status
  error(result)
end

if ~ isempty(find(abs(fname) == 46))
	fname = fname(:,1:find(abs(fname) == 46)-1) ;
end

evalin('base',fname) ;


% MJ 2/9/99: replace clear function by clear ff_
% MJ 4/7/00: change the path of dynare_m
% MJ 02/26/01: replaced local variable x by fname
% MJ 09/19/01: evaluates mod script in 'base' workspace
