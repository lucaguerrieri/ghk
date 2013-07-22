function z=ff_simul1(x,y1,m,state,o1,o2,o3,k)
  global fname_ dr_

  fh = str2func([fname_ '_ff']);
  tempx = [x(dr_.order_var(o1:o2))-dr_.ys(dr_.order_var(o1:o2)); state];
  y2 = m+dr_.ghx(o3:end,:)*tempx;
  z=feval(fh,[y1; x; y2(k)]);
