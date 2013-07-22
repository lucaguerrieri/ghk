function z=ff_simul(x,y1,m,state,tempu,o1,o2,o3,delta)
  global fname_ dr_

  fh = str2func([fname_ '_ff']);
  tempx = [x(dr.order_var(o1:o2))-dr.ys(dr.order_var(o1:o2)); state];
  tempxx = kron(tempx,tempx);
  tempxu = kron(tempx,tempu);
%  y2 = dr.ys(dr.order_var)+dr.ghs2/2+dr.ghx*tempx+dr.ghu*tempu+0.5*(dr.ghxx*tempxx+dr.ghuu*tempuu)+dr.ghxu*tempxu;
  y2(dr.order_var(o3:end)) = m+dr.ghx(o3:end,:)*tempx+0.5*dr.ghxx(o3:end,:)*tempxx+dr.ghxu(o3:end,:)*tempxu;
  z=feval(fh,(y1; x; y2))+delta;
