%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this program creates the analytical jacobian
% of a model specified using a dynare .mod file
%
% remember to run dynare before running this program
% in order to genrate the _ff file (remove the line near the top
% that initializes the size of the output vector) and 
% create the model solution .mat file needed below
%
% the program requires the model to contain
% only one lead and one lag of all variables
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

% contains parameter file
paramfile_name = 'param_simple_k.m';
composite_paramfile_name = 'param_composite2.m';

modfilename = 'ghk2_montecarlo';

% model solution needs to contain iy_ lgx_ lgy_ 
% create this mat file separately, after running dynare
load model_solution_ghk2_montecarlo
%load model_solution_rbc
%load model_solution_ghk3

if (size(iy_,1)>3)
    error('Can only process models with one lead and one lag')
end

% file creating steady states 
% as functions of fundamental parameters
steady_state_file ='steady_state_ghk2_montecarlo';
%steady_state_file = 'steady_state_rbc';
%steady_state_file = 'steady_state_ghk3';

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define model's parameter as sym objects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get names of parameters
paramfile = textread(paramfile_name,'%s', ...
                     'delimiter','\n','whitespace','','bufsize',40000);
paramfile_nlines = size(paramfile,1);
list_param = [];
for i=1:paramfile_nlines 
    tokens = tokenize(char(paramfile(i)),'=');
    if length(tokens)>1
        first_token = char(tokens(1));
        if first_token(1)~='%'
           list_param = [list_param,tokens(1)];
        end
    end
end
list_param= char(list_param);

% declare parameters as symbolic variables
nparams = size(list_param,1);
for i = 1:nparams
    eval(['global ',list_param(i,:),';'])
    eval(['syms ',list_param(i,:),';'])
end

% get names of composite parameters
paramfile = textread(composite_paramfile_name,'%s', ...
                     'delimiter','\n','whitespace','','bufsize',40000);
paramfile_nlines = size(paramfile,1);
list_comp_param = [];
indx = 0;
for i=1:paramfile_nlines 
    tokens = tokenize(char(paramfile(i)),'=');
    if length(tokens)>1
        first_token = char(tokens(1));
        if first_token(1)~='%'
           indx = indx+1;
           list_comp_param = [list_comp_param,tokens(1)];
           eval(['global ',char(list_comp_param(indx)),';'])
           eval(['syms comp_param',num2str(indx)]);
           eval(['syms ',char(list_comp_param(indx))]);
           eval(['comp_param',num2str(indx),'=',char(tokens(3)),';'])
        end
    end
end

list_comp_param = char(list_comp_param);
ncomp_params = size(list_comp_param,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define state vector as sym object
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eval(steady_state_file)

neq = size(lgy_,1);
nex = size(lgx_,1);

% initialize matrices
% leads
hl1 = zeros(neq);
% contemp.
h = hl1;
% lags
hm1 = hl1;
% innovations
d = zeros(neq,nex);

% find non-zero columns of hm1
lag_cols = find(iy_(1,:)~=0);
% find non-zero columns of h
con_cols = find(iy_(2,:)~=0);
% find non-zero columns of hl1
lea_cols = find(iy_(3,:)~=0);


% find number of entries for y vector used 
% by Dynare to stack the state vectors across time
ny = length(find(iy_~=0));
% build steady state y
y = ys_(lag_cols);
y = [y;ys_(con_cols)];
y = [y;ys_(lea_cols)];


global it_
it_ = 2;

for i=1:nex
    eval(['syms ',lgx_(i,:)]);
    eval(['ex_(1,i)=',lgx_(i,:),';']);
end
global ex_

for i=1:ny
    eval(['syms y',num2str(i)]);
    eval(['yvec(i) = y',num2str(i),';']);
end

% get dynamic equations from file
eval(['z=',modfilename,'_ff(yvec);'])

% create jacobian file
fid = fopen('get_hmat.f90','w+');
fprintf(fid,['SUBROUTINE get_hmat( hmat, hrows, hcols, psi_array, psin ) \n']);

fprintf(fid,['IMPLICIT NONE \n']);

fprintf(fid,['! \n']);
fprintf(fid,['! Input Variables \n']);
fprintf(fid,['! \n']);
fprintf(fid,['INTEGER, INTENT(IN) :: hrows, hcols, psin \n']);
fprintf(fid,['DOUBLE PRECISION, INTENT(IN), DIMENSION(psin) :: psi_array \n']);
fprintf(fid,['DOUBLE PRECISION, INTENT(OUT), DIMENSION(hrows,hcols) :: hmat \n']);
for k=1:ny
    fprintf(fid,['DOUBLE PRECISION :: y',num2str(k),' \n']);
end
for k = 1:nex
    fprintf(fid,['DOUBLE PRECISION :: ',char(sym(ex_(k))),' \n']);
end
for i=1:nparams
    fprintf(fid,['DOUBLE PRECISION :: ',list_param(i,:),' \n']);
end
ncparams = size(list_comp_param,1);
for i=1:ncparams
    fprintf(fid,['DOUBLE PRECISION :: ',list_comp_param(i,:),' \n']);
end

% include translation between parameter vector and actual parameters
for i=1:nparams
    fprintf(fid,[list_param(i,:),' = psi_array(',num2str(i),') \n']);
end

% include formulae for composite parameters
for i=1:ncparams
    fprintf(fid,[list_comp_param(i,:),'=',eval(['syms2fortran(comp_param',num2str(i),')']),' \n']);
end

% write down the steady state formulae,
% so that the jacobian will be evaluated at the steady state
for k=1:ny
    fprintf(fid,['y',num2str(k),' = ',syms2fortran(y(k)),' \n']);
end
for k = 1:nex
    fprintf(fid,[char(sym(ex_(k))),' = 0.d0 \n']);
end


fprintf(fid,'\n\n');
fprintf(fid,'! lags \n');
fprintf(fid,'hmat(:,:)=0.0d0 \n');

nlag_cols = length(lag_cols);
for i=1:nlag_cols
    display(['variable ',num2str(i),' of ',num2str(ny)])   
    for j=1:neq
        eval(['partial_deriv = diff(z(j),y',num2str(i),');'])
        if partial_deriv~=0         
         fprintf(fid,['hmat(',num2str(j),',',num2str(lag_cols(i)),')=',syms2fortran(partial_deriv),' \n']);
        end
    end
end


fprintf(fid,'\n\n');
fprintf(fid,'! current variables \n');
ncon_cols = length(con_cols);
offset = nlag_cols;
for i=1:ncon_cols
    display(['variable ',num2str(i+offset),' of ',num2str(ny)]) 
    for j=1:neq
        eval(['partial_deriv = diff(z(j),y',num2str(i+offset),');'])
        if partial_deriv~=0        
         fprintf(fid,['hmat(',num2str(j),',',num2str(con_cols(i)+neq+nex),')=',syms2fortran(partial_deriv),';\n']);
        end
    end
end

fprintf(fid,'\n\n');
fprintf(fid,'! leads \n');
nlea_cols = length(lea_cols);
offset = nlag_cols+ncon_cols;
for i=1:nlea_cols
    display(['variable ',num2str(i+offset),' of ',num2str(ny)]) 
    for j=1:neq
        eval(['partial_deriv = diff(z(j),y',num2str(i+offset),');'])
        if partial_deriv~=0         
            fprintf(fid,['hmat(',num2str(j),',',num2str(lea_cols(i)+2*(neq+nex)),')=',syms2fortran(partial_deriv),' \n']);
        end
    end
end


fprintf(fid,'\n\n');
fprintf(fid,'! innovations \n');
% deal with innovations
offset = neq;
for i = 1:nex
    display(['innovation ',num2str(i),' of ',num2str(nex)]) 
    for j = 1:neq
        partial_deriv = diff( z(j),ex_(i) );
        if partial_deriv~=0
            fprintf(fid,['hmat(',num2str(j),',',num2str(offset+i),')=',syms2fortran(partial_deriv),';\n']);
        end
    end
    % deal with extra definitions
    fprintf(fid,['hmat(',num2str(offset+i),',',num2str(offset+i+neq+nex),')=1.0d0 \n']);

end


fprintf( fid, '\n\n');

fprintf( fid, 'END SUBROUTINE' );
fclose(fid);


%% create file containing helper structures

% create jacobian file
fid = fopen([modfilename,'_data.m'],'w+');
fprintf(fid,['endog_=char( ... \n']);

for i=1:neq
    fprintf(fid,[char(39),lgy_(i,:),char(39),', ... \n']);
end
for i=1:nex-1
    fprintf(fid,[char(39),lgx_(i,:),char(39),', ... \n']);
end
fprintf(fid,[char(39),lgx_(nex,:),char(39),'); \n\n\n']);

fprintf(fid,['exog_=char( ... \n']);
for i=1:nex-1
    fprintf(fid,[char(39),lgx_(i,:),char(39),', ... \n']);
end
fprintf(fid,[char(39),lgx_(nex,:),char(39),'); \n\n\n']);

fprintf(fid,['param_=char( ... \n']);
for i=1:nparams-1
    fprintf(fid,[char(39),list_param(i,:),char(39),', ... \n']);
end
fprintf(fid,[char(39),list_param(nparams,:),char(39),'); \n\n\n']);

fprintf(fid,['comp_param_=char( ... \n']);

ncparams = size(list_comp_param,1);
for i=1:ncparams-1
    fprintf(fid,[char(39),list_comp_param(i,:),char(39),', ... \n']);
end
fprintf(fid,[char(39),list_comp_param(ncparams,:),char(39),'); \n\n\n']);

fclose(fid);


%% now reimport the file and split any line longer than 132 characters
fortranfile = textread('get_hmat.f90','%s', ...
                      'delimiter','\n','whitespace','','bufsize',40000);

nlines = size(fortranfile,1);


fid = fopen('get_hmat.f90','w+');

for i=1:nlines
    line = char(fortranfile(i));
    line_length = length(line);
    if line_length<=132
        %leave line unchanged
        fprintf(fid,[line,'\n']);
    elseif (line_length > 132 & line_length < 260)
        fprintf(fid,[line(1:130),'&\n']);
        fprintf(fid,['&',line(131:end),'\n']);
    else
        error('line too long for tapenade. Add cases in run_get_hsym_fortran program')
    end
end
fclose(fid);

%% helper file for fortran program
fid = fopen([modfilename,'_fortran.txt'],'w+');



% labels of endogenous variables
for i=1:neq
fprintf(fid,['endog(',num2str(i),')=',char(39),strtrim(lgy_(i,:)),char(39),' \n']);
end
indx = 0;
for i=neq+1:neq+nex;
    indx = indx+1;
    fprintf(fid,['endog(',num2str(i),')=',char(39),strtrim(lgx_(indx,:)),char(39),' \n']);
end
fprintf(fid,[' \n\n\n']);

param_simple_k
nparams = size(list_param,1);
% values for parameters
for i=1:nparams
       fprintf(fid,['psi(',num2str(i),')=',num2str(eval(list_param(i,:)),'%.16fd0 '),'    ! ',list_param(i,:),'\n']);
end


fclose(fid);