% makes transition matrices out of ghx and ghu

function [A,B] = transition_matrix(dr)
  global ykmin_ exo_nbr 
  nx = size(dr.ghx,2);
  kstate = dr.kstate;
  ikx = [dr.nstatic+1:dr.nstatic+dr.npred];
  
  A = zeros(nx,nx);
  k0 = kstate(find(kstate(:,2) <= ykmin_+1),:);
  i0 = find(k0(:,2) == ykmin_+1);
  A(i0,:) = dr.ghx(ikx,:);
  B = zeros(nx,exo_nbr);
  B(i0,:) = dr.ghu(ikx,:);
  for i=ykmin_:-1:2
    i1 = find(k0(:,2) == i);
    n1 = size(i1,1);
    j = zeros(n1,1);
    for j1 = 1:n1
      j(j1) = find(k0(i0,1)==k0(i1(j1),1));
    end
    A(i1,i0(j))=eye(n1);
    i0 = i1;
  end
