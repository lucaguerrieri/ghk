% makes transition matrices out of ghx and ghu for Kalman filter
% still needs to eliminate unobserved static variables
% order of variables s p b f p_1

function [A,B] = kalman_transition_matrix(dr)
  global ykmin_ exo_nbr 
  nx = size(dr.ghx,2)+dr.nfwrd+dr.nstatic;
  kstate = dr.kstate;
  ikx = [dr.nstatic+1:dr.nstatic+dr.npred];
  
  A = zeros(nx,nx);
  k0 = kstate(find(kstate(:,2) <= ykmin_+1),:);
  i0 = find(k0(:,2) == ykmin_+1);
  n0 = size(dr.ghx,1);
  A(1:n0,dr.nstatic+1:dr.nstatic+dr.npred) = dr.ghx(:,1:dr.npred);
  A(1:n0,dr.nstatic+dr.npred+dr.nfwrd+1:end) = dr.ghx(:,dr.npred+1:end);
  B = zeros(nx,exo_nbr);
  B(1:n0,:) = dr.ghu;
  offset_col = dr.nstatic;
  for i=ykmin_:-1:2
    i1 = find(k0(:,2) == i);
    n1 = size(i1,1);
    j = zeros(n1,1);
    for j1 = 1:n1
      j(j1) = find(k0(i0,1)==k0(i1(j1),1));
    end
    if i == ykmin_-1
      offset_col = dr.nstatic+dr.nfwrd;
    end
    A(n0+i1-dr.npred,offset_col+i0(j))=eye(n1);
    i0 = i1;
  end
