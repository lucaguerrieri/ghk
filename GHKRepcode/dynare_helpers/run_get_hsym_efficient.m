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
for i=1:paramfile_nlines 
    tokens = tokenize(char(paramfile(i)),'=');
    if length(tokens)>1
        first_token = char(tokens(1));
        if first_token(1)~='%'
           list_comp_param = [list_comp_param,tokens(1)];
        end
    end
end

% declare composite parameters as symbolic variables
list_comp_param = char(list_comp_param);
ncomp_params = size(list_comp_param,1);
for i = 1:ncomp_params
    eval(['global ',list_comp_param(i,:),';'])
    eval(['syms ',list_comp_param(i,:),';'])
end


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
fid = fopen([modfilename,'_jacobian.m'],'w+');

% fprintf(fid,'function [hm1 h hl1] =get_jacobian(params) \n');
% for i=1:nparams
%     fprintf(fid,['% params(',num2str(i),') = ',list_param(i,:),' \n']);
% end
% 
% 
% % declare the parameters as global variables
% for i=1:nparams
%     fprintf(fid,[list_param(i,:),'= params(',num2str(i),'); \n']);
% end
% 
% % insert definitions of composite parameters
% for i=1:paramfile_nlines
%     fprintf(fid,[char(paramfile(i)),'\n']);
% end


% write down the steady state formulae,
% so that the jacobian will be evaluated at the steady state
for k=1:ny
    fprintf(fid,['y',num2str(k),' = ',char(sym(y(k))),'; \n']);
end
for k = 1:nex
    fprintf(fid,[char(sym(ex_(k))),' = 0; \n']);
end


fprintf(fid,'\n\n');
fprintf(fid,'%% lags \n');
fprintf(fid,['hm1=zeros(',num2str(neq+nex),'); \n']);
nlag_cols = length(lag_cols);
for i=1:nlag_cols
    display(['variable ',num2str(i),' of ',num2str(ny)])   
    for j=1:neq
        eval(['partial_deriv = diff(z(j),y',num2str(i),');'])
        if partial_deriv~=0         
         fprintf(fid,['hm1(',num2str(j),',',num2str(lag_cols(i)),')=',char(sym(partial_deriv)),';\n']);
        end
    end
end


fprintf(fid,'\n\n');
fprintf(fid,'%% current variables \n');
fprintf(fid,['h=zeros(',num2str(neq+nex),'); \n']);
ncon_cols = length(con_cols);
offset = nlag_cols;
for i=1:ncon_cols
    display(['variable ',num2str(i+offset),' of ',num2str(ny)]) 
    for j=1:neq
        eval(['partial_deriv = diff(z(j),y',num2str(i+offset),');'])
        if partial_deriv~=0        
         fprintf(fid,['h(',num2str(j),',',num2str(con_cols(i)),')=',char(sym(partial_deriv)),';\n']);
        end
    end
end

fprintf(fid,'\n\n');
fprintf(fid,'%% leads \n');
fprintf(fid,['hl1=zeros(',num2str(neq+nex),'); \n']);
nlea_cols = length(lea_cols);
offset = nlag_cols+ncon_cols;
for i=1:nlea_cols
    display(['variable ',num2str(i+offset),' of ',num2str(ny)]) 
    for j=1:neq
        eval(['partial_deriv = diff(z(j),y',num2str(i+offset),');'])
        if partial_deriv~=0         
            fprintf(fid,['hl1(',num2str(j),',',num2str(lea_cols(i)),')=',char(sym(partial_deriv)),';\n']);
        end
    end
end


fprintf(fid,'\n\n');
fprintf(fid,'%% innovations \n');
% deal with innovations
offset = neq;
for i = 1:nex
    display(['innovation ',num2str(i),' of ',num2str(nex)]) 
    for j = 1:neq
        partial_deriv = diff( z(j),ex_(i) );
        if partial_deriv~=0
           fprintf(fid,['hm1(',num2str(j),',',num2str(offset+i),')=',char(sym(partial_deriv)),';\n']);
        end
    end
    % deal with extra definitions
    fprintf(fid,['h(',num2str(offset+i),',',num2str(offset+i),')=1;\n']);

end

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
for i=1:nparams
    fprintf(fid,[char(39),list_param(i,:),char(39),', ... \n']);
end
ncparams = size(list_comp_param,1);
for i=1:ncparams
    fprintf(fid,[char(39),list_comp_param(i,:),char(39),', ... \n']);
end
fprintf(fid,[char(39),lgx_(nex,:),char(39),'); \n\n\n']);



fclose(fid);