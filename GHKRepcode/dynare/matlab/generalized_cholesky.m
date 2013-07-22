% proc gmchol(A);
%    /* calculates the Gill-Murray generalized choleski decomposition */
%    /* input matrix A must be non-singular and symmetric */
%    /* Author: Jeff Gill. Part of the Hessian Project. */
%    local i,j,k,n,sum,R,theta_j,norm_A,phi_j,delta,xi_j,gamm,E,beta_j;

function AA = generalized_cholesky(A);

n = rows(A);
R = eye(n);
E = zeros(n,n);
norm_A = max(sum(abs(A))');
gamm   = max(abs(diag(A))); 
delta = max([eps*norm_A;eps]);

for j = 1:n; 
    theta_j = 0;
    for i=1:n;
        somme = 0;
        for k=1:i-1;	
	        somme = somme + R(k,i)*R(k,j);
        end;
        R(i,j) = (A(i,j) - somme)/R(i,i);
	    if (A(i,j) -somme) > theta_j;
	        theta_j = A(i,j) - somme;
        end;
	    if i > j;
	        R(i,j) = 0;
        end;
    end;
    somme = 0;
    for k=1:j-1;	
        somme = somme + R(k,j)^2;
    end;
    phi_j = A(j,j) - somme;
    if j+1 <= n;
        xi_j = max(abs(A((j+1):n,j)));
    else;
        xi_j = abs(A(n,j));
    end;
    beta_j = sqrt(max([gamm ; (xi_j/n) ; eps]));
    if all(delta >= [abs(phi_j);((theta_j^2)/(beta_j^2))]);
        E(j,j) = delta - phi_j;
    elseif all(abs(phi_j) >= [((delta^2)/(beta_j^2));delta]);
        E(j,j) = abs(phi_j) - phi_j;
    elseif all(((theta_j^2)/(beta_j^2)) >= [delta;abs(phi_j)]);
        E(j,j) = ((theta_j^2)/(beta_j^2)) - phi_j;
    end;
    R(j,j) = sqrt(A(j,j) - somme + E(j,j));
end;
AA = R'*R;



