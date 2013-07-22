% Copyright (C) 2004 Marco Ratto
% adapted from Michel Juillard original rutine hessian.m
%
%  [hessian_mat, gg, htol1, ihh, hh_mat0] = mr_hessian(func,x,hflag,htol0,varargin)
%
%  numerical gradient and Hessian, with 'automatic' check of numerical
%  error 
%
%  func =  name of the function: func must give two outputs: 
%    - the log-likelihood AND the single contributions at times t=1,...,T
%    of the log-likelihood to compute outer product gradient
%  x = parameter values
%  hflag = 0, Hessian computed with outer product gradient, one point
%  increments for partial derivatives in  gradients
%  hflag = 1, 'mixed' Hessian: diagonal elements computed with numerical second order derivatives
%             with correlation structure as from outer product gradient;
%             two point evaluation of derivatives for partial derivatives
%             in gradients
%  hflag = 2, full numerical Hessian, computes second order partial derivatives
%          uses Abramowitz and Stegun (1965) formulas 25.3.24 and 25.3.27
%          p. 884.
%  htol0 = 'precision' of increment of function values for numerical
%  derivatives
%
%  varargin: other parameters of func
%

function [hessian_mat, gg, htol1, ihh, hh_mat0] = mr_hessian(func,x,hflag,htol0,varargin)
global gstep_ bayestopt_
persistent h1 htol

if isempty(htol), htol = 1.e-4; end
func = str2func(func);
[f0, ff0]=feval(func,x,varargin{:});
n=size(x,1);
h2=bayestopt_.ub-bayestopt_.lb;
%h1=max(abs(x),gstep_*ones(n,1))*eps^(1/3);
%h1=max(abs(x),sqrt(gstep_)*ones(n,1))*eps^(1/6);
if isempty(h1),
    h1=max(abs(x),sqrt(gstep_)*ones(n,1))*eps^(1/4);
end

if htol0<htol, 
    htol=htol0;
end
xh1=x;
f1=zeros(size(f0,1),n);
f_1=f1;
ff1=zeros(size(ff0));
ff_1=ff1;

%for i=1:n,
i=0;
while i<n,
    i=i+1;
    h10=h1(i);
    hcheck=0;
    dx=[];
    xh1(i)=x(i)+h1(i);
    [fx, ffx]=feval(func,xh1,varargin{:});
    it=1;
    dx=(fx-f0);
    ic=0;
    if abs(dx)>(2*htol),
        c=mr_nlincon(xh1,varargin{:});
        while c
            h1(i)=h1(i)*0.9;
            xh1(i)=x(i)+h1(i);
            c=mr_nlincon(xh1,varargin{:});        
            ic=1;
        end   
        if ic,
            [fx, ffx]=feval(func,xh1,varargin{:});
            dx=(fx-f0);
        end
    end
    
    icount = 0;
    h0=h1(i);
    while (abs(dx(it))<0.5*htol | abs(dx(it))>(2*htol)) & icount<10 & ic==0,
        %while abs(dx(it))<0.5*htol & icount< 10 & ic==0,
        icount=icount+1;
        %if abs(dx(it)) ~= 0,
        if abs(dx(it))<0.5*htol
            if abs(dx(it)) ~= 0,
                h1(i)=min(0.3*abs(x(i)), 0.9*htol/abs(dx(it))*h1(i));
            else
                h1(i)=2.1*h1(i);
            end
            xh1(i)=x(i)+h1(i);
            c=mr_nlincon(xh1,varargin{:});
            while c
                h1(i)=h1(i)*0.9;
                xh1(i)=x(i)+h1(i);
                c=mr_nlincon(xh1,varargin{:});        
                ic=1;
            end  
            [fx, ffx]=feval(func,xh1,varargin{:});
        end
        if abs(dx(it))>(2*htol),
            h1(i)= htol/abs(dx(it))*h1(i);
            xh1(i)=x(i)+h1(i);
            [fx, ffx]=feval(func,xh1,varargin{:});
            while (fx-f0)==0,
                h1(i)= h1(i)*2;
                xh1(i)=x(i)+h1(i);
                [fx, ffx]=feval(func,xh1,varargin{:});
                ic=1;
            end
        end
        it=it+1;
        dx(it)=(fx-f0);
        h0(it)=h1(i);
        if h1(i)<1.e-12*min(1,h2(i)),
            ic=1;
            hcheck=1;
        end
        %else
        % h1(i)=1;
        % ic=1;
        %end
    end
    %     if (it>2 & dx(1)<10^(log10(htol)/2)) ,        
    %         [dum, is]=sort(h0); 
    %         if find(diff(sign(diff(dx(is)))));
    %             hcheck=1;
    %         end           
    %     elseif (it>3 & dx(1)>10^(log10(htol)/2)) ,        
    %         [dum, is]=sort(h0); 
    %         if find(diff(sign(diff(dx(is(1:end-1))))));
    %             hcheck=1;
    %         end           
    %     end
    f1(:,i)=fx;
    ff1=ffx;
    if hflag,  % two point based derivatives
        xh1(i)=x(i)-h1(i);
        c=mr_nlincon(xh1,varargin{:});
        ic=0;
        while c
            h1(i)=h1(i)*0.9;
            xh1(i)=x(i)-h1(i);
            c=mr_nlincon(xh1,varargin{:});  
            ic = 1;
        end    
        [fx, ffx]=feval(func,xh1,varargin{:});
        f_1(:,i)=fx;
        ff_1=ffx;
        if ic,
            xh1(i)=x(i)+h1(i);
            [f1(:,i), ff1]=feval(func,xh1,varargin{:});
        end
        ggh(:,i)=(ff1-ff_1)./(2.*h1(i));
    else
        ggh(:,i)=(ff1-ff0)./h1(i);
    end
    xh1(i)=x(i);
    if hcheck & htol<1,
        htol=min(1,max(min(abs(dx))*2,htol*10));
        h1(i)=h10;
        i=0;
    end
    save hess
end

h_1=h1;
xh1=x;
xh_1=xh1;

if hflag,
    gg=(f1'-f_1')./(2.*h1);
else
    gg=(f1'-f0)./h1;
end

if hflag==2,
    gg=(f1'-f_1')./(2.*h1);
    hessian_mat = zeros(size(f0,1),n*n);
    for i=1:n
        if i > 1
            k=[i:n:n*(i-1)];
            hessian_mat(:,(i-1)*n+1:(i-1)*n+i-1)=hessian_mat(:,k);
        end 
        hessian_mat(:,(i-1)*n+i)=(f1(:,i)+f_1(:,i)-2*f0)./(h1(i)*h_1(i));
        temp=f1+f_1-f0*ones(1,n);
        for j=i+1:n
            xh1(i)=x(i)+h1(i);
            xh1(j)=x(j)+h_1(j);
            xh_1(i)=x(i)-h1(i);
            xh_1(j)=x(j)-h_1(j);
            %hessian_mat(:,(i-1)*n+j)=-(-feval(func,xh1,varargin{:})-feval(func,xh_1,varargin{:})+temp(:,i)+temp(:,j))./(2*h1(i)*h_1(j));
            %temp1 = feval(func,xh1,varargin{:});
            c=mr_nlincon(xh1,varargin{:});
            lam=1;
            while c, 
                lam=lam*0.9;
                xh1(i)=x(i)+h1(i)*lam;
                xh1(j)=x(j)+h_1(j)*lam;
                %disp( ['hessian warning cross ', num2str(c) ]), 
                c=mr_nlincon(xh1,varargin{:});
            end
            temp1 = f0+(feval(func,xh1,varargin{:})-f0)/lam;
            
            %temp2 = feval(func,xh_1,varargin{:});
            c=mr_nlincon(xh_1,varargin{:});
            while c, 
                lam=lam*0.9;
                xh_1(i)=x(i)-h1(i)*lam;
                xh_1(j)=x(j)-h_1(j)*lam;
                %disp( ['hessian warning cross ', num2str(c) ]), 
                c=mr_nlincon(xh_1,varargin{:});
            end
            temp2 = f0+(feval(func,xh_1,varargin{:})-f0)/lam;
            
            hessian_mat(:,(i-1)*n+j)=-(-temp1 -temp2+temp(:,i)+temp(:,j))./(2*h1(i)*h_1(j));
            xh1(i)=x(i);
            xh1(j)=x(j);
            xh_1(i)=x(i);
            xh_1(j)=x(j);
            j=j+1;
            save hess
        end
        i=i+1;
    end
    
elseif hflag==1,
    hessian_mat = zeros(size(f0,1),n*n);
    for i=1:n,
        dum = (f1(:,i)+f_1(:,i)-2*f0)./(h1(i)*h_1(i));
        if dum>eps,
            hessian_mat(:,(i-1)*n+i)=dum;
        else
            hessian_mat(:,(i-1)*n+i)=max(eps, gg(i)^2);
        end                        
    end
    %hessian_mat2=hh_mat(:)';
end

gga=ggh.*kron(ones(size(ff1)),2.*h1');  % re-scaled gradient
hh_mat=gga'*gga;  % rescaled outer product hessian 
hh_mat0=ggh'*ggh;  % outer product hessian
A=diag(2.*h1);  % rescaling matrix
if hflag>0 & min(eig(reshape(hessian_mat,n,n)))>0,
    hh0 = A*reshape(hessian_mat,n,n)*A';  %rescaled second order derivatives
    sd0=sqrt(diag(hh0));   %rescaled 'standard errors' using second order derivatives
    sd=sqrt(diag(hh_mat));  %rescaled 'standard errors' using outer product
    hh_mat=hh_mat./(sd*sd').*(sd0*sd0');  %rescaled outer product with 'true' std's
    hh0 = reshape(hessian_mat,n,n);  % second order derivatives
    sd0=sqrt(diag(hh0));   % 'standard errors' using second order derivatives
    sd=sqrt(diag(hh_mat0));  % 'standard errors' using outer product
    hh_mat0=hh_mat0./(sd*sd').*(sd0*sd0');  % rescaled outer product with 'true' std's
end
igg=inv(hh_mat);  % inverted rescaled outer product hessian
ihh=A'*igg*A;  % inverted outer product hessian
if hflag==0,
    hessian_mat=hh_mat0(:);
end

if isnan(hessian_mat),
    hh_mat0=eye(length(hh_mat0));
    ihh=hh_mat0;
    hessian_mat=hh_mat0(:);    
end
hh1=h1;
htol1=htol;
save hess
% 11/25/03 SA Created from Hessian_sparse (removed sparse)


