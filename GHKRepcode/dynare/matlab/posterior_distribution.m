function [borneinf,bornesup,x1,x2,f1,f2,top,nam,texnam] = ...
    posterior_distribution(indx,number_of_mh_files,first_mh_file,first_line,...
					number_of_blocks,number_of_simulations,...
					number_of_simulations_per_file,TeX);
% stephane.adjemian@cepremap.cnrs.fr [07-15-2004]

global fname_ bayestopt_ estim_params_ options_

number_of_grid_points = 2^9;      % 2^9 = 512 !... Must be a power of two.
bandwidth = 0;                    % Rule of thumb optimal bandwidth parameter.
kernel_function = 'gaussian';     % Gaussian kernel for Fast Fourrier Transform approximaton.  
             

mcsimulations = zeros(number_of_simulations,1);

if number_of_blocks == 1
  EndOfFile = number_of_simulations_per_file(first_mh_file+1)-first_line+1;
  instr = [fname_ '_mh' int2str(first_mh_file)];
  eval(['load ' instr]);
  clear post2 logpo2;
  mcsimulations(1:EndOfFile) = x2(first_line:end,indx);
  OldEndOfFile = EndOfFile;
  for f = first_mh_file+1:number_of_mh_files
    NewEndOfFile = number_of_simulations_per_file(f+1);
    instr = [fname_ '_mh' int2str(f)];
    eval(['load ' instr]);
    clear post2 logpo2;
    mcsimulations(OldEndOfFile+1:OldEndOfFile+NewEndOfFile) = x2(:,indx);
    OldEndOfFile = OldEndOfFile + NewEndOfFile;
  end
  clear x2;
else
  EndOfFile = number_of_simulations_per_file(first_mh_file+1)-first_line+1;
  NewStartLine = 0;
  inst = [fname_ '_mh' int2str(first_mh_file)];
  for b = 1:number_of_blocks
    instr = [inst '_blck' int2str(b)];
    eval(['load ' instr]);
    clear post2 logpo2;
    mcsimulations(NewStartLine+1:NewStartLine+EndOfFile,1) = x2(first_line:end,indx);
    NewStartLine = NewStartLine + EndOfFile;
  end
  for f = first_mh_file+1:number_of_mh_files
    EndOfFile = number_of_simulations_per_file(f+1);
    inst = [fname_ '_mh' int2str(f)];
    for B = 1:number_of_blocks
      instr = [inst '_blck' int2str(b)];
      eval(['load ' instr]);
      clear post2 logpo2;
      mcsimulations(NewStartLine+1:NewStartLine+EndOfFile,1) = x2(:,indx);
      NewStartLine = NewStartLine + EndOfFile;
    end
  end
  clear x2;
end

[nam,texnam] = get_the_name(indx,TeX);

%% Kernel estimator of the posterior density:
optimal_bandwidth = mh_optimal_bandwidth(mcsimulations,number_of_simulations,bandwidth,kernel_function); 
[x1,f1] = kernel_density_estimate(mcsimulations,number_of_grid_points,...
                                  optimal_bandwidth,kernel_function);
binf1 = x1(1);
bsup1 = x1(length(x1));

%% Prior density:
[x2,f2,abscissa,dens,binf2,bsup2] = draw_prior_density(indx);
clear abscissa dens;

borneinf = min(binf1,binf2);
bornesup = max(bsup1,bsup2);
top      = max([max(f1);max(f2)]);