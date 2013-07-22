  function [MeanX,MedianX,StdX,DistribX,HPDX] = GetPosteriorStatistics(gend,B,infotype)
  % stephane.adjemian@gmail.com [05-26-2005]
  global options_ fname_ lgy_ lgx_ lgy_TeX_ lgx_TeX_ dr_ bayestopt_ oo_
  
  deciles = [round(0.1*B) ...
	     round(0.2*B)...
	     round(0.3*B)...
	     round(0.4*B)...
	     round(0.5*B)...
	     round(0.6*B)...
	     round(0.7*B)...
	     round(0.8*B)...
	     round(0.9*B)];
 
 tmp = zeros(B,1);
  
  if strcmpi(infotype,'SmoothedVariables')
    PrintOnScreen = 'Smooth variables';
    NumberOfX = size(lgy_,1);
    SaveName = 'SmoothedVariables';
    XNames = lgy_(dr_.order_var,:);
    NonGenericName = '_smooth';
    varname = 'smooth';
  elseif strcmpi(infotype,'SmoothedShocks')
    PrintOnScreen = 'Smooth structural shocks...';
    NumberOfX = size(lgx_,1);
    SaveName = 'SmoothedShocks';
    XNames = lgx_;
    NonGenericName = '_innovation';
    varname = 'innov';
  elseif strcmpi(infotype,'SmoothedObservationErrors');  
    PrintOnScreen = 'Smooth measurement error...';
    NumberOfX = size(options_.varobs,1);
    SaveName = 'SmoothedMeasurementErrors';
    XNames = options_.varobs;
    NonGenericName = '_error';
    varname = 'error';
  elseif strcmpi(infotype,'FilteredVariables');  
    PrintOnScreen = 'Filtered variables...';
    NumberOfX = size(lgy_,1);
    SaveName = 'FilteredVariables';
    XNames = lgy_(dr_.order_var,:);
    NonGenericName = '_filter';
    varname = 'filter';
  end
  eval(['sfile = size(dir(''' fname_ NonGenericName '*.mat''),1);'])
  disp(['MH: ' PrintOnScreen '...'])
  MeanX = zeros(NumberOfX,gend);
  MedianX = zeros(NumberOfX,gend);
  StdX = zeros(NumberOfX,gend);
  DistribX = zeros(NumberOfX,gend,9);
  HPDX = zeros(NumberOfX,gend,2);
  for i = 1:NumberOfX
     for t = 1:gend
	   StartLine = 0;
	   for file = 1:sfile;
	     instr = [fname_ NonGenericName int2str(file)];
	     eval(['load ' instr]);
	     eval(['X = stock_' varname ';'])
	     MeanX(i,t) = MeanX(i,t)+sum(X(i,t,:),3);
	     DeProfundis = size(X,3);
	     tmp(StartLine+1:StartLine+DeProfundis) = squeeze(X(i,t,:)); 
	     StartLine = StartLine+DeProfundis;
	   end
	   tmp = sort(tmp);
	   MedianX(i,t) = tmp(round(B*0.5));
	   StdX(i,t) = std(tmp);
	   DistribX(i,t,:) = reshape(tmp(deciles),1,1,9);
	   tt = floor(options_.mh_conf_sig*B);
	   a = 1; 
	   b = tt;
	   tmp2 = [1;tt;tmp(tt)-tmp(1)];
	   while b <= B
	     tmp1 = [a;b;tmp(b)-tmp(a)];
	     a = a + 1;
	     b = b + 1;
	     if tmp1(3,1) < tmp2(3,1)
	       tmp2 = tmp1;     
	     end    
	  end
	  HPDX(i,t,1) = tmp(tmp2(1,1));
	  HPDX(i,t,2) = tmp(tmp2(2,1));
    end
    disp(['    Variable: ' deblank(XNames(i,:))]);	
  end
  clear X;
  MeanX = MeanX/B;
  for i=1:NumberOfX
    eval(['oo_.Posterior' SaveName '.Mean.' deblank(XNames(i,:)) ' = MeanX(i,:)'';']);
    eval(['oo_.Posterior' SaveName '.Median.' deblank(XNames(i,:)) ' = MedianX(i,:)'';']);
    eval(['oo_.Posterior' SaveName '.Std.' deblank(XNames(i,:)) ' = StdX(i,:)'';']);
    eval(['oo_.Posterior' SaveName '.Distribution.' deblank(XNames(i,:)) ' = squeeze(DistribX(i,:,:))'';']);
    eval(['oo_.Posterior' SaveName '.HPDinf.' deblank(XNames(i,:)) ' = squeeze(HPDX(i,:,1))'';']);
    eval(['oo_.Posterior' SaveName '.HPDsup.' deblank(XNames(i,:)) ' = squeeze(HPDX(i,:,2))'';']);
  end
  disp(['MH: ' PrintOnScreen ', done!'])
  disp(' ')