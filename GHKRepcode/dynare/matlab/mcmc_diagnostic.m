function mcmc_diagnostic(x2)
% 
% See S. Brooks and Gelman [1998, Journal ]
%

global estim_params_ lgx_ options_

nblck = size(x2,3);
npar  = size(x2,2);
nruns = size(x2,1)

origin = 1000;
step_size = 100;

det_W = zeros(nruns,1);
det_B = zeros(nruns,1);
det_V = zeros(nruns,1);
R     = zeros(nruns,1);
ligne = 0;

for iter = origin:step_size:nruns;
    ligne = ligne + 1; 
    W = zeros(npar,npar);
    B = zeros(npar,npar);
    linea = ceil(0.5*iter);
    n = iter-linea+1;
    muB = mean(mean(x2(linea:iter,:,:),3),1)';  
    for j = 1:nblck;
        muW  = mean(x2(linea:iter,:,j))';
        B = B + (muW-muB)*(muW-muB)';
        for t = linea:iter;
            W = W + (x2(t,:,j)'-muW)*(x2(t,:,j)-muW');     
        end;    
    end;
    W = inv(nblck*(n-1))*W;
    B = n*inv(nblck-1)*B;
    V = inv(n)*(n-1)*W + (1+inv(nblck))*B/n;
    det_W(ligne,1) = det(W);
    det_B(ligne,1) = det(B);
    det_V(ligne,1) = det(V);
    lambda = max(eig(inv(n*W)*B));
    R(ligne,1) = (n-1)/n + lambda*(nblck+1)/nblck;
end;    

det_W = det_W(1:ligne,1);
det_V = det_V(1:ligne,1);
R     = R(1:ligne,1);

figure('Name','MCMC multivariate diagnostic');
subplot(2,1,1);
plot(origin:step_size:nruns,R,'-k','linewidth',2);
subplot(2,1,2);
plot(origin:step_size:nruns,det_W,'--r','linewidth',2);
hold on;
plot(origin:step_size:nruns,det_V,'--b','linewidth',2);
hold off;


[R_interval,R_c1,R_c2] = mcmcdiags(x2,0.2,3,4,origin,step_size);
pages = floor(npar/4);
k = 0;

% pages
for i = 1:pages;
    figure('Name','MCMC univariate diagnostic');
    boxplot = 1;
    for j = 1:4;
        k = k+1;
        if k <= estim_params_.nvx;
            vname = deblank(lgx_(estim_params_.var_exo(k,1),:));
            nam=['SE_{',vname,'}'];
        elseif  k <= (estim_params_.nvx+estim_params_.nvn);
            vname = deblank(options_.varobs(estim_params_.var_endo(k-estim_params_.nvx,1),:));
            nam=['SE_{EOBS_',vname,'}'];
        elseif  k <= (estim_params_.nvx+estim_params_.nvn+estim_params_.ncx);
            jj = k - (estim_params_.nvx+estim_params_.nvn);
            k1 = estim_params_.corrx(jj,1);
            k2 = estim_params_.corrx(jj,2);
            vname = [deblank(lgx_(k1,:)) ',' deblank(lgx_(k2,:))];
            nam=['CC_{',vname,'}'];
        elseif  k <= (estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+ ...
		estim_params_.ncn)
            jj = k - (estim_params_.nvx+estim_params_.nvn+estim_params_.ncx);
            k1 = estim_params_.corrn(jj,1);
            k2 = estim_params_.corrn(jj,2);
            vname = [deblank(lgx_(k1,:)) ',' deblank(lgx_(k2,:))];
            nam=['CC_{EOBS_',vname,'}'];
        else
            jj = k - (estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+estim_params_.ncn);
            nam = deblank(estim_params_.param_names(jj,:));
        end
        xx = origin:step_size:nruns;
        for crit = 1:3;
            if crit == 1;
                plt1 = R_interval(:,k,1);
                plt2 = R_interval(:,k,2);
                namnam  = [nam , ', R_{interval}']; 
            elseif crit == 2;
                plt1 = R_c1(:,k,1);
                plt2 = R_c1(:,k,2);
                namnam  = [nam , ', R_{c1}'];
            elseif crit == 3;    
                plt1 = R_c2(:,k,1);
                plt2 = R_c2(:,k,2);
                namnam  = [nam , ', R_{c2}'];
            end;
            subplot(4,3,boxplot);
            plot(xx,plt1,'-b');
            hold on;
            plot(xx,plt2,'-r');
            hold off;
            title(namnam,'Interpreter','none');
            boxplot = boxplot + 1;
        end;
    end;
end;
reste = npar-k;
figure('Name','MCMC univariate diagnostic');
boxplot = 1;
for j = 1:reste;
    k = k+1;
    if k <= estim_params_.nvx;
        vname = deblank(lgx_(estim_params_.var_exo(k,1),:));
        nam=['SE_{',vname,'}'];
    elseif  k <= (estim_params_.nvx+estim_params_.nvn);
        vname = deblank(options_.varobs(estim_params_.var_endo(k-estim_params_.nvx,1),:));
        nam=['SE_{EOBS_',vname,'}'];
    elseif  k <= (estim_params_.nvx+estim_params_.nvn+estim_params_.ncx);
        jj = k - (estim_params_.nvx+estim_params_.nvn);
        k1 = estim_params_.corrx(jj,1);
        k2 = estim_params_.corrx(jj,2);
        vname = [deblank(lgx_(k1,:)) ',' deblank(lgx_(k2,:))];
        nam=['CC_{',vname,'}'];
    elseif  k <= (estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+ ...
		estim_params_.ncn)
        jj = k - (estim_params_.nvx+estim_params_.nvn+estim_params_.ncx);
        k1 = estim_params_.corrn(jj,1);
        k2 = estim_params_.corrn(jj,2);
        vname = [deblank(lgx_(k1,:)) ',' deblank(lgx_(k2,:))];
        nam = ['CC_{EOBS_',vname,'}'];
    else
        jj = k - (estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+estim_params_.ncn);
        nam = deblank(estim_params_.param_names(jj,:));
    end
    xx = origin:step_size:nruns;
    for crit = 1:3;
        if crit == 1;
            plt1 = R_interval(:,k,1);
            plt2 = R_interval(:,k,2);
            namnam  = [nam , ', R_{interval}']; 
        elseif crit == 2;
            plt1 = R_c1(:,k,1);
            plt2 = R_c1(:,k,2);
            namnam  = [nam , ', R_{c1}'];
        elseif crit == 3;    
            plt1 = R_c2(:,k,1);
            plt2 = R_c2(:,k,2);
            namnam  = [nam , ', R_{c2}'];
        end;
        subplot(reste,3,boxplot);
        plot(xx,plt1,'-b');
        hold on;
        plot(xx,plt2,'-r');
        hold off;
        title(namnam,'Interpreter','none');
        boxplot = boxplot + 1;
    end;
end;    
    
    
    

