function [xparam1, hh, gg, fval] = newrat(func0,x,hh,gg,flag,varargin)
%
%  [xparam1, hh, gg, fval] = newrat(func0,x,hh,gg,flag,varargin)
%
%  Standard Newton search
%
%  flag = 0, to begin with a pure gradient search (the program shifts to
%  pure Newton when improvement of pure gradient is below ftol)
%
%  flag = 1, to start with the compelte Newton search

  global bayestopt_
icount=0;
nx=length(x);
xparam1=x;
lamtol=1.e-7;
ftol=1.e-5;
options=optimset('fminunc');
options.MaxFunEvals=200;
options.TolFun= 1.0000e-005;
options.MaxIter=1;
options.LargeScale='off';

optim_options=optimset('fminsearch');
optim_options.display='iter';
optim_options.MaxFunEvals=1000;
optim_options.TolFun= 1.0000e-003;
optim_options.TolX= 1.0000e-006;


func = str2func(func0);
fval0=feval(func,x,varargin{:});
if isempty(hh)
    [dum, gg]=mr_hessian(func0,x,flag,varargin{:});
    hh = reshape(dum,nx,nx);
end
disp(['Gradient norm ',num2str(norm(gg))])
disp(['Minimum Hessian eigenvalue ',num2str(min(eig(hh)))])
disp(['Maximum Hessian eigenvalue ',num2str(max(eig(hh)))])
g=gg;
h{1}=hh;
check=0;
if max(eig(hh))<0, disp('Negative definite Hessian! Local maximum!'), pause, end,
while norm(gg)>1.e-3 & check==0,
    icount=icount+1;
    bayestopt_.penalty = fval0(icount);
    disp([' '])
    disp(['Iteration ',num2str(icount)])
    x0=xparam1-inv(hh)*gg;
    c=mr_nlincon(x0,varargin{:},1);
    lam=1;
    while c
        lam=lam*0.9;
        x0=xparam1-inv(hh)*gg.*lam;
        c=mr_nlincon(x0,varargin{:},1);
    end        
    fval=feval(func,x0,varargin{:});
%     if (fval0(icount)-fval)<ftol & flag==0,
%         fvala=fval;
%         x0a=x0;
%         disp('Try to modify Hessian')
%         x0=xparam1-inv(gg*gg')*gg;
%         c=mr_nlincon(x0,varargin{:},1);
%         lam=1;
%         while c
%             lam=lam*0.9;
%             x0=xparam1-inv(gg*gg')*gg.*lam;
%             c=mr_nlincon(x0,varargin{:},1);
%         end        
%         fval=feval(func,x0,varargin{:});
%         if fvala<=fval, 
%             x0=x0a;
%             fval=fvala;
%         end            
%     end            
    if (fval0(icount)-fval)<ftol,
        disp('Try line search')
        [lam,fval,EXITFLAG,OUTPUT,GRAD,HESSIAN]=fminunc(@lsearch, 0, options, func, xparam1, inv(hh)*gg , varargin{:});
        x0=xparam1-inv(hh)*gg.*lam;
    end
    if (fval0(icount)-fval)<ftol & flag==1,
        fvala=fval;
        x0a=x0;
        disp('Try gradient direction')
        [lam,fval,EXITFLAG,OUTPUT,GRAD,HESSIAN]=fminunc(@lsearch, 0, options, func, xparam1, gg , varargin{:});
        if fvala<=fval, 
            x0=x0a;
            fval=fvala;
        else
            x0=xparam1-gg*lam;
            if (fval0(icount)-fval)>ftol,
                flag=0;
            end
        end
    end
    if (fval0(icount)-fval)<ftol,
        fvala=fval;
        x0a=x0;
        disp('Try some simplex iterations')
        [x0,fval,EXITFLAG,OUTPUT] = fminsearch(func, xparam1, optim_options, varargin{:});
        if fvala<fval, 
            x0=x0a;
            fval=fvala;
        else
            lam = NaN;
        end
    end
    if (fval0(icount)-fval)<ftol*ftol & flag==1;,
        %         if fvala<fval,
        %             fval=fvala;
        %             x0=x0a;
        %         end
        disp('No further improvement is possible!')
        check=1;
    else
        
        if (fval0(icount)-fval)<ftol & flag==0,
            flag=1;
        end
        
        xparam1=x0;
        x(:,icount+1)=xparam1;
        fval0(icount+1)=fval;
        disp(['LAMBDA        ',num2str(lam)])
        %disp(['DX norm       ',num2str(norm(inv(hh)*gg.*lam))])
        disp(['DX norm       ',num2str(norm(x(:,end)-x(:,end-1)))])
        disp(['FVAL          ',num2str(fval)])
        disp(['Improvement   ',num2str(fval0(icount)-fval)])
        
        if norm(x(:,icount)-xparam1)>1.e-12,
            %[dum, gg]=hessian('mj_optmumlik',xparam1,gend,data,1);
            [dum, gg]=mr_hessian(func0,xparam1,flag,varargin{:});
            hh = reshape(dum,nx,nx);
        end
        disp(['Gradient norm  ',num2str(norm(gg))])
        disp(['Minimum Hessian eigenvalue ',num2str(min(eig(hh)))])
        disp(['Maximum Hessian eigenvalue ',num2str(max(eig(hh)))])
        if max(eig(hh))<0, disp('Negative definite Hessian! Local maximum!'), pause, end,
        
        h{icount+1}=hh;
        g(:,icount+1)=gg;
        save m1 x h g fval0
    end
end

return

%  
function f00 = lsearch(lam,func,x,dx,varargin)


x0=x-dx*lam;
f00=feval(func,x0,varargin{:});






