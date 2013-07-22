function z=ff_simul2(x,y1,m,state,o1,o2,o3,k)
  global fname_ dr_

  fh = str2func([fname_ '_ff']);
  tempx = [x(dr_.order_var(o1:o2))-dr_.ys(dr_.order_var(o1:o2)); state];
  tempxx = kron(tempx,tempx);
%  y2 = dr.ys(dr.order_var)+dr.ghs2/2+dr.ghx*tempx+dr.ghu*tempu+0.5*(dr.ghxx*tempxx+dr.ghuu*tempuu)+dr.ghxu*tempxu;
  y2 = m+dr_.ghx(o3:end,:)*tempx+0.5*dr_.ghxx(o3:end,:)*tempxx;
  z=feval(fh,[y1; x; y2(k)])+dr_.fuu;
