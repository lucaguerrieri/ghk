function PosteriorOddsTable = model_comparison(ModelNames,ModelPriors)
% 05-30-2005
%
% type is a string  = Laplace
%                   = ModifiedHarmonicMean
% ModelPriors is a m*1 column vector
% ModelNames is m*1 cell array

global oo_ options_ fname_

tmp_oo = oo_;
if isfield(options_,'model_comparison_approximation')
  type = options_.model_comparison_approximation;
else
  type = '';
end

if strcmp(type,'Laplace')
  type = 'LaplaceApproximation';
end

NumberOfModels = size(ModelNames,1);
MarginalLogDensity = zeros(NumberOfModels,1);
if isempty(ModelPriors)
  ModelPriors = ones(NumberOfModels,1);
end

% Get the estimates of the (logged) marginal densities

init_loop = 1;
if isempty(type)
  if strcmp(ModelNames{1},fname_)
    oo_ = tmp_oo;
  else
    load([ModelNames{1} '_results.mat' ],'oo_');
  end
  try
    type = 'LaplaceApproximation';
    eval(['MarginalLogDensity(1) =' ...
	  ' oo_.MarginalDensity.LaplaceApproximation;']);
  catch
    try
      type = 'ModifiedHarmonicMean';
      eval(['MarginalLogDensity(1) = oo_.MarginalDensity.ModifiedHarmonicMean;']); 
    catch
      disp(['CompareModels :: I cant''t find any marginal density approximation associated to model ' ModelNames{1}])
      return
    end
  end
  init_loop = 2;
end
for i = init_loop:NumberOfModels
  if strcmp(ModelNames{i},fname_)
    oo_ = tmp_oo;
  else
    load([ModelNames{i} '_results.mat' ],'oo_');
  end
  try
    eval(['MarginalLogDensity(i) = oo_.MarginalDensity.' type ';']) 
  catch
    if strcmpi(type,'LaplaceApproximation')
      disp(['CompareModels :: I cant''t find the Laplace approximation associated to model ' ModelNames{i}])
      return
    elseif strcmpi(type,'ModifiedHarmonicMean')
      disp(['CompareModels :: I cant''t find the modified harmonic mean estimate associated to model ' ModelNames{i}])
      return
    end
  end
end

%MarginalDensity = exp(MarginalLogDensity);
LogConstantOfIntegration  = sum(log(ModelPriors)+MarginalLogDensity);
PosteriorProbabilities = log(ModelPriors) + MarginalLogDensity - ...
    LogConstantOfIntegration;

PosteriorOddsTable = exp(repmat(PosteriorProbabilities,1,NumberOfModels)- ...
			 repmat(PosteriorProbabilities',NumberOfModels,1));

% Now I display the posterior probabilities:
if NumberOfModels == 2
    disp(' ')
    disp(['Posterior odd (' ModelNames{1} '/' ModelNames{2} ') =  ' num2str(PosteriorOddsTable(1,2))])
    disp(' ')
else
    disp(' ')
    disp(' Posterior probabilities:')
    for i=1:NumberOfModels
            disp([ 'Model ' int2str(i) ' (' ModelNames{i} ') = ' num2str(PosteriorProbabilities(i))])
    end
    disp(' ')
end
oo_ = tmp_oo;