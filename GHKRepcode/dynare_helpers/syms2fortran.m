function s=syms2fortran(sym_arg)

% test whether the symbolic expression contains symbols or only a number
is_number = 1;
try
    % if sym_arg contains symbols, matlab will not be able to evaluate it
    % and is_number will then be assigned the value 0
    eval(sym_arg);
catch
    is_number = 0;
end

if (is_number == 0 )
    % format for 
    s = eval(['strrep(strrep(deblank(fortran(sym_arg)),''t0 = '',''''),''     #'','''')']);
else
    % format for fortran double precision
    s = num2str(eval(sym_arg),'%.16fd0');
    
end