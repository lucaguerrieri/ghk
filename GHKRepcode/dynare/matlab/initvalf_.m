function initvalf_(fname,period,varargin)
  global lgy_ lgx_ lgr_ y_ ex_ valf_ ykmin_ xkmin_ ykmax_ xkmax_ iter_ ...
      y_start_date ex_start_date freq_ start_date

  if ~isempty(strfind(upper(fname),'.XLS'))
    [data,names_v]=xlsread(fname);
  else
    load(fname);
  end
  
  if length(period) == 2
    period = dy_date(period(1),period(2));
  end
  
  if period - max(ykmin_,xkmin_) < 0
    error(['INITVALF_: not enough data points in database for number of' ...
	   ' lags. Start later!'])
  end
  
  if nargin > 2
    if strcmp(upper(varargin{1}),'SERIES')
      series = 1 ;
    elseif strcmp(upper(varargin{1}),'MAT')
      series = 0 ;
    else
      error(['INITVALF: unknown option ' varargin{1}])
    end
  else
    series = 0 ;
  end
  
  y1 = floor((period-ykmin_)/freq_);
  p1 = period-ykmin_-freq_*y1;
  y_start_date(2) = start_date(2) + p1-1;
  if y_start_date(2) > freq_
    y_start_date(2) = y_start_date(2) - freq_;
    y1 = y1 + 1;
  end
  y_start_date(1) = start_date(1)+y1;
  
  y1 = floor((period-xkmin_)/freq_);
  p1 = period-xkmin_-freq_*y1;
  ex_start_date(2) = start_date(2) + p1-1;
  if y_start_date(2) > freq_
    ex_start_date(2) = ex_start_date(2) - freq_;
    y1 = y1 + 1;
  end
  ex_start_date(1) = start_date(1)+y1;
  
  clear y1, p1;
  
  valf_ = 1;
  y_ = [];
  ex_ = [];
  
  for i=1:size(lgy_,1)
    if series == 1
      x = eval([lgy_(i,:) '(period-ykmin_:period+iter_+ykmax_-1);']);
      y_ = [y_; x'];
    else
      k = strmatch(upper(lgy_(i,:)),names_v,'exact');
      if isempty(k)
	error(['INITVALF: ' lgy_(i,:) ' not found'])
      end
      x = data(:,k);
      y_ = [y_; x(period-ykmin_:period+iter_+ykmax_-1)']; 
    end
  end
  
  for i=1:size(lgx_,1)
    if series == 1
      x = eval([lgx_(i,:) '(period-xkmin_:period+iter_+xkmax_-1);']);
      ex_ = [ex_ x];
    else
      k = strmatch(upper(lgx_(i,:)),names_v,'exact');
      if isempty(k)
	error(['INITVALF: ' lgx_(i,:) ' not found'])
      end
      x = data(:,k);
      ex_ = [ex_ x(period-xkmin_:period+iter_+xkmax_-1)]; 
    end
  end
    
% $$$   if any(isnan(y_(:,1))) | any(isnan(ex_(1,:)))
% $$$     error('INITVALF: missing value in first period')
% $$$   end
% $$$   
% $$$   if any(isnan(y_(:,end))) | any(isnan(ex_(end,:)))
% $$$     error('INITVALF: missing value in last period')
% $$$   end
  
% 8/23/01 MJ changed argument 'FILE' to 'MAT'











