function y=dy_date(year,period)
  global start_date freq_
  
  y = freq_*(year-start_date(1))+period-start_date(2)+1;
  
  