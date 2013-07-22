function [R_interval,R_c1,R_c2] = mcmcdiags(x2,alpha,c1,c2,origin,step_size);
% This function computes univariate mcmc diagnostics. See S. Brooks and A.
% Gelman "General Methods for Monitoring Convergence of Iterative
% Simulations".
% 
% stephane.adjemian@cepremap.cnrs.fr [01/05/2004] 

global fname_


f1 = [fname_ '_MH_diagnostics_1'];
f2 = [fname_ '_MH_diagnostics_2'];


npar  = size(x2,2);  
nruns = size(x2,1);
nblck = size(x2,3);
ligne = 0;
LIGNE = 0;

R_interval = zeros(nruns,npar,2);
MU1   = zeros(1,npar,nblck);

BigArray  = cat(2,x2,zeros(nruns,npar,nblck));
for i = 1:npar;
    for j = 1:nblck;
        [tmp,indx] = sort(BigArray(:,i,j));
        BigArray(:,i,j) = tmp;
        BigArray(:,i+npar,j) = indx;
    end;
end;
save(f1,'BigArray');
clear BigArray;
BigMatrix = zeros(nruns*nblck,2*npar);
tmp = (1:nruns)';
for i=1:npar;
    for j=1:nblck;
        BigMatrix((j-1)*nruns+1:j*nruns,i) = x2(:,i,j);
        BigMatrix((j-1)*nruns+1:j*nruns,npar+i) = tmp;
    end;
    pmt = sortrows([BigMatrix(:,i),BigMatrix(:,npar+i)],1);
    BigMatrix(:,i) = pmt(:,1);
    BigMatrix(:,npar+i) = pmt(:,2);    
end;   
save(f2,'BigMatrix');
clear BigMatrix;

load(f2);
for iter  =  origin:step_size:nruns;
    ligne = ligne + 1;
    linea = ceil(0.5*iter);
    n = iter-linea+1;
    cinf = round(nblck*n*alpha/2);
    csup = round(nblck*n*(1-alpha/2));
    for i=1:npar;
        linie = find(BigMatrix(:,npar+i) >= linea & ...
            BigMatrix(:,npar+i)<=iter);
        tmp = BigMatrix(linie,i); 
        R_interval(ligne,i,1) = tmp(csup)-tmp(cinf);
    end;
end;
clear BigMatrix;LIGNE = ligne;ligne = 0;
R_interval = R_interval(1:LIGNE,:,:);
load(f1);
for iter  = origin:step_size:nruns;
    ligne = ligne+1;
    linea = ceil(0.5*iter);
    n = iter-linea+1;
    cinf = round(n*alpha/2);
    csup = round(n*(1-alpha/2));
    for i=1:npar;
        for j=1:nblck;
            linie = find(BigArray(:,npar+i,j) >= linea & ...
                BigArray(:,npar+i,j)<=iter);
            tmp = BigArray(linie,i,j); 
            R_interval(ligne,i,2) = tmp(csup)-tmp(cinf);
        end;    
    end;
end;
clear BigArray;ligne = 0;

R_c1 = zeros(LIGNE,npar,2);
R_c2 = zeros(LIGNE,npar,2);

for iter  = origin:step_size:nruns;
    ligne  = ligne + 1; 
    linea  = ceil(0.5*iter);
    n      = iter-linea+1;
    selvec = linea:iter;
    MEAN1  = zeros(n,npar,nblck);
    MEAN2  = zeros(n,npar,nblck);
    tmp = x2(selvec,:,:);                    
    mu1 = mean(mean(tmp,3),1);
    mu2 = mean(tmp,1);
    for j=1:nblck;
        MU1(1,:,j) = mu1;
    end;    
    for t = 1:n;
        MEAN1(t,:,:) = MU1;
        MEAN2(t,:,:) = mu2;
    end;    
    R_c1(ligne,:,1) = sum(sum(abs(tmp-MEAN1).^c1,3),1)/(nblck*n-1);
    R_c2(ligne,:,1) = sum(sum(abs(tmp-MEAN1).^c2,3),1)/(nblck*n-1);
    R_c1(ligne,:,2) = sum(sum(abs(tmp-MEAN2).^c1,3),1)/(nblck*(n-1));
    R_c2(ligne,:,2) = sum(sum(abs(tmp-MEAN2).^c2,3),1)/(nblck*(n-1));
end;    

R_c1 = R_c1.^(1/c1);
R_c2 = R_c2.^(1/c2);

delete([f1 '.mat'])
delete([f2 '.mat'])

% 2004/1/12 SA corrected bug line 59