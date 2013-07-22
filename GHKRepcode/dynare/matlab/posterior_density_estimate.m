function [abscissa,f,h] = posterior_density_estimate(data,number_of_grid_points,bandwidth,kernel_function) 
%%   
%%  This function aims at estimating posterior univariate densities from realisations of a Metropolis-Hastings 
%%  algorithm. A kernel density estimator is used (see Silverman [1986]) and the main task of this function is 
%%  to obtain an optimal bandwidth parameter. 
%% 
%%  * Silverman [1986], "Density estimation for statistics and data analysis". 
%%  * M. Skold and G.O. Roberts [2003], "Density estimation for the Metropolis-Hastings algorithm". 
%%
%%  The last section is adapted from Anders Holtsberg's matlab toolbox (stixbox).
%%
%%  stephane.adjemian@cepremap.cnrs.fr [01/16/2004]. 

if size(data,2) > 1 & size(data,1) == 1; 
    data = data'; 
elseif size(data,2)>1 & size(data,1)>1; 
    error('density_estimate: data must be a one dimensional array.'); 
end;
test = log(number_of_grid_points)/log(2);
if ( abs(test-round(test)) > 10^(-12));
    error('The number of grid points must be a power of 2.');
end;

n = size(data,1); 


%% KERNEL SPECIFICATION...

if strcmp(kernel_function,'gaussian'); 
    k    = inline('inv(sqrt(2*pi))*exp(-0.5*x.^2)'); 
    k2   = inline('inv(sqrt(2*pi))*(-exp(-0.5*x.^2)+(x.^2).*exp(-0.5*x.^2))'); % second derivate of the gaussian kernel 
    k4   = inline('inv(sqrt(2*pi))*(3*exp(-0.5*x.^2)-6*(x.^2).*exp(-0.5*x.^2)+(x.^4).*exp(-0.5*x.^2))'); % fourth derivate... 
    k6   = inline('inv(sqrt(2*pi))*(-15*exp(-0.5*x.^2)+45*(x.^2).*exp(-0.5*x.^2)-15*(x.^4).*exp(-0.5*x.^2)+(x.^6).*exp(-0.5*x.^2))'); % sixth derivate... 
    mu02 = inv(2*sqrt(pi)); 
    mu21 = 1; 
elseif strcmp(kernel_function,'uniform'); 
    k    = inline('0.5*(abs(x) <= 1)'); 
    mu02 = 0.5; 
    mu21 = 1/3; 
elseif strcmp(kernel_function,'triangle'); 
    k    = inline('(1-abs(x)).*(abs(x) <= 1)'); 
    mu02 = 2/3; 
    mu21 = 1/6; 
elseif strcmp(kernel_function,'epanechnikov'); 
    k    = inline('0.75*(1-x.^2).*(abs(x) <= 1)'); 
    mu02 = 3/5; 
    mu21 = 1/5;     
elseif strcmp(kernel_function,'quartic'); 
    k    = inline('0.9375*((1-x.^2).^2).*(abs(x) <= 1)'); 
    mu02 = 15/21; 
    mu21 = 1/7; 
elseif strcmp(kernel_function,'triweight'); 
    k    = inline('1.09375*((1-x.^2).^3).*(abs(x) <= 1)'); 
    k2   = inline('(105/4*(1-x.^2).*x.^2-105/16*(1-x.^2).^2).*(abs(x) <= 1)'); 
    k4   = inline('(-1575/4*x.^2+315/4).*(abs(x) <= 1)'); 
    k6   = inline('(-1575/2).*(abs(x) <= 1)'); 
    mu02 = 350/429; 
    mu21 = 1/9;     
elseif strcmp(kernel_function,'cosinus'); 
    k    = inline('(pi/4)*cos((pi/2)*x).*(abs(x) <= 1)'); 
    k2   = inline('(-1/16*cos(pi*x/2)*pi^3).*(abs(x) <= 1)'); 
    k4   = inline('(1/64*cos(pi*x/2)*pi^5).*(abs(x) <= 1)'); 
    k6   = inline('(-1/256*cos(pi*x/2)*pi^7).*(abs(x) <= 1)'); 
    mu02 = (pi^2)/16; 
    mu21 = (pi^2-8)/pi^2;     
end;     


%% OPTIMAL BANDWIDTH PARAMETER....

if bandwidth == 0;  %  Rule of thumb bandwidth parameter (Silverman [1986] corrected by 
                    %  Skold and Roberts [2003] for Metropolis-Hastings). 
    sigma = std(data); 
    h = 2*sigma*(sqrt(pi)*mu02/(12*(mu21^2)*n))^(1/5); % Silverman's optimal bandwidth parameter. 
    A = 0; 
    for i=1:n; 
        j = i; 
        while j<= n & data(j,1)==data(i,1); 
            j = j+1; 
        end;     
        A = A + 2*(j-i) - 1; 
    end; 
    A = A/n; 
    h = h*A^(1/5); % correction 
elseif bandwidth == -1;     % Adaptation of the Sheather and Jones [1991] plug-in estimation of the optimal bandwidth 
                            % parameter for metropolis hastings algorithm. 
    if strcmp(kernel_function,'uniform')      | ... 
       strcmp(kernel_function,'triangle')     | ... 
       strcmp(kernel_function,'epanechnikov') | ... 
       strcmp(kernel_function,'quartic'); 
       error('I can''t compute the optimal bandwidth with this kernel... Try the gaussian, triweight or cosinus kernels.'); 
    end;         
    sigma = std(data); 
    A = 0; 
    for i=1:n; 
        j = i; 
        while j<= n & data(j,1)==data(i,1); 
            j = j+1; 
        end;     
        A = A + 2*(j-i) - 1; 
    end; 
    A = A/n; 
    Itilda4 = 8*7*6*5/(((2*sigma)^9)*sqrt(pi)); 
    g3      = abs(2*A*k6(0)/(mu21*Itilda4*n))^(1/9); 
    Ihat3 = 0; 
    for i=1:n; 
        Ihat3 = Ihat3 + sum(k6((data(i,1)-data)/g3)); 
    end;     
    Ihat3 = -Ihat3/((n^2)*g3^7); 
    g2      = abs(2*A*k4(0)/(mu21*Ihat3*n))^(1/7); 
    Ihat2 = 0; 
    for i=1:n; 
        Ihat2 = Ihat2 + sum(k4((data(i)-data)/g2)); 
    end;     
    Ihat2 = Ihat2/((n^2)*g2^5); 
    h       = (A*mu02/(n*Ihat2*mu21^2))^(1/5);    % equation (22) in Skold and Roberts [2003] --> h_{MH} 
elseif bandwidth == -2;     % Bump killing... We construct local bandwith parameters in order to remove 
                            % spurious bumps introduced by long rejecting periods.   
    if strcmp(kernel_function,'uniform')      | ... 
       strcmp(kernel_function,'triangle')     | ... 
       strcmp(kernel_function,'epanechnikov') | ... 
       strcmp(kernel_function,'quartic'); 
        error('I can''t compute the optimal bandwidth with this kernel... Try the gaussian, triweight or cosinus kernels.'); 
    end;         
    sigma = std(data); 
    A = 0; 
    T = zeros(n,1); 
    for i=1:n; 
        j = i; 
        while j<= n & data(j,1)==data(i,1); 
            j = j+1; 
        end;     
        T(i) = (j-i); 
        A = A + 2*T(i) - 1; 
    end; 
    A = A/n; 
    Itilda4 = 8*7*6*5/(((2*sigma)^9)*sqrt(pi)); 
    g3      = abs(2*A*k6(0)/(mu21*Itilda4*n))^(1/9); 
    Ihat3 = 0; 
    for i=1:n; 
        Ihat3 = Ihat3 + sum(k6((data(i,1)-data)/g3)); 
    end;     
    Ihat3 = -Ihat3/((n^2)*g3^7); 
    g2      = abs(2*A*k4(0)/(mu21*Ihat3*n))^(1/7); 
    Ihat2 = 0; 
    for i=1:n; 
        Ihat2 = Ihat2 + sum(k4((data(i)-data)/g2)); 
    end;     
    Ihat2 = Ihat2/((n^2)*g2^5); 
    h = ((2*T-1)*mu02/(n*Ihat2*mu21^2)).^(1/5); % Note that h is a column vector (local banwidth parameters). 
elseif bandwidth > 0; 
    h = bandwidth; 
else; 
    error('density_estimate: bandwidth must be positive or equal to 0,-1 or -2.'); 
end; 

%% COMPUTE DENSITY ESTIMATE, using the optimal bandwidth parameter.
%%
%% This section is adapted from Anders Holtsberg's matlab toolbox
%% (stixbox --> plotdens.m).


a  = min(data) - (max(data)-min(data))/3;
b  = max(data) + (max(data)-min(data))/3;
abscissa = linspace(a,b,number_of_grid_points)';
d  = abscissa(2)-abscissa(1); 
xi = zeros(number_of_grid_points,1);
xa = (data-a)/(b-a)*number_of_grid_points; 
for i = 1:n;
    indx = floor(xa(i));
    temp = xa(i)-indx;
    xi(indx+[1 2]) = xi(indx+[1 2]) + [1-temp,temp]';
end;    
xk = [-number_of_grid_points:number_of_grid_points-1]'*d;
kk = k(xk/h);
kk = kk / (sum(kk)*d*n);
f = ifft(fft(fftshift(kk)).*fft([xi ;zeros(size(xi))]));
f = real(f(1:number_of_grid_points));