function hessian2(xparam1,gend,data)
% source Schorfheide
%/* filename:    lbdhess.g
%** description: The program computes the hessianfit at the posterior mode 
%** created:     05/05/00
%/******************************************************** 
%/* Compute hessianfit, element by element, fine tune with dxscale
%   Compute hessianfit for two step-seize (dx) and take average to prevent singularity  
%/******************************************************** 
npara = length(xparam1);
para = xparam1;
ndx = 1;%6
%dx =  exp(-seqa(6,2,ndx));
%dx =  exp([-6:-2:-16]);
dx =  exp([-10]);
hessianfit = zeros( npara, npara );
gradx = zeros(ndx,1);
grady = zeros(ndx,1);
gradxy = zeros(ndx,1);
hessdiag = zeros(ndx,1);
dxscale = ones(npara,1);
%  dxscale(5,1)=10;
%  dxscale(13,1)=10;
%  dxscale(17,1)=10;
%  dxscale(8,1)=.10;
%  dxscale(11,1)=.10;
%  dxscale(31,1)=.10;
  

%/* Compute Diagonal elements first
%*/
seli = 1;
fx  = mj_optmumlik(para,gend,data,1);
%do until seli > npara;
while seli <= npara;
%     locate 1,1;
%     "hessianfit Element    (" seli seli ")";
     i=1;
     while i <= ndx;
      paradx = para;
      parady = para;
      paradx(seli) = paradx(seli) + dx(i)*dxscale(seli);
      parady(seli) = parady(seli) - dx(i)*dxscale(seli);
      paradxdy = paradx;
      paradxdy(seli) = paradxdy(seli) - dx(i)*dxscale(seli);
%     fx  = optmumlik20(para,q,h2,aa,mf,gitno,gst,gobs,gend,gstart,nlags,rawdata,lpd);
fdx  = mj_optmumlik(paradx,gend,data,1);
%      fdx  = optmumlik20(paradx,q,h2,aa,mf,gitno,gst,gobs,gend,gstart,nlags,rawdata,lpd);
fdy  = mj_optmumlik(parady,gend,data,1);
   %   fdy  = optmumlik20(parady,q,h2,aa,mf,gitno,gst,gobs,gend,gstart,nlags,rawdata,lpd);
%     fdxdy  = optmumlik20(paradxdy,q,h2,aa,mf,gitno,gst,gobs,gend,gstart,nlags,rawdata,lpd);
      fdxdy  = fx;
      gradx(i) = -( fx - fdx )/ (dx(i)*dxscale(seli));
      grady(i) = ( fx - fdy )/ (dx(i)*dxscale(seli));
      gradxy(i) = -(fx -fdxdy)/ sqrt( (dx(i)*dxscale(seli))^2 + (dx(i)*dxscale(seli))^2 );
      hessdiag(i) = -( 2*fx - fdx - fdy)/(dx(i)*dxscale(seli))^2; 
      hessdiag(i) = -( fx - fdx - fdy + fdxdy )/(dx(i)*dx(i)*dxscale(seli)*dxscale(seli));
      i = i+1;
  end;
%     "Values";
%     -hessdiag';
     hessianfit(seli,seli) = -1*(hessdiag(1));
%     hessianfit(seli,seli) = -0.5*(hessdiag(3)+hessdiag(4));
%     locate 6,1;
%     "Value Used:" hessianfit[seli,seli];
   seli = seli+1
end;

diag(hessianfit)

%/* Now compute off-diagonal elements
%** Make sure that correlations are between -1 and 1
%** errorij contains the index of elements that are invalid
%*/
errorij = [ 0 0 0];

seli = 1;
for seli = 1:npara;
%   selj = seli+1;
   for selj =seli+1:npara;
       disp([seli selj]);
%     locate 1,1;
%     "hessianfit Element    (" seli selj ")";
     i=1;
     while i <= ndx;
      paradx = para;
      parady = para;
      paradx(seli) = paradx(seli) + dx(i)*dxscale(seli);
      parady(selj) = parady(selj) - dx(i)*dxscale(selj);
      paradxdy = paradx;
      paradxdy(selj) = paradxdy(selj) - dx(i)*dxscale(selj);
%      fx  = optmumlik20(para,q,h2,aa,mf,gitno,gst,gobs,gend,gstart,nlags,rawdata,lpd);
fdx  = mj_optmumlik(paradx,gend,data,1);
%      fdx  = optmumlik20(paradx,q,h2,aa,mf,gitno,gst,gobs,gend,gstart,nlags,rawdata,lpd);
fdy  = mj_optmumlik(parady,gend,data,1);
%      fdy  = optmumlik20(parady,q,h2,aa,mf,gitno,gst,gobs,gend,gstart,nlags,rawdata,lpd);
fdy  = mj_optmumlik(paradxdy,gend,data,1);
%      fdxdy  = optmumlik20(paradxdy,q,h2,aa,mf,gitno,gst,gobs,gend,gstart,nlags,rawdata,lpd);
      gradx(i) = -( fx - fdx )/ (dx(i)*dxscale(seli));
      grady(i) = ( fx - fdy )/ (dx(i)*dxscale(selj));
      gradxy(i) = -(fx -fdxdy)/ sqrt( (dx(i)*dxscale(selj))^2 + (dx(i)*dxscale(seli))^2 );
      hessdiag(i) = -( 2*fx - fdx - fdy)/(dx(i)*dxscale(seli))^2; 
      hessdiag(i) = -( fx - fdx - fdy + fdxdy )/(dx(i)*dx(i)*dxscale(seli)*dxscale(selj));
      i = i+1;
     end;
%     "Values";
%     -hessdiag';

%     hessianfit(seli,selj) = -0.5*(hessdiag(3)+hessdiag(4));
     hessianfit(seli,selj) = -1*(hessdiag(1));
     
     if ( hessianfit(seli,seli) == 0 ) | (hessianfit(selj,selj) == 0);
        corrij = 0;
     else;
        corrij = hessianfit(seli,selj)/sqrt(hessianfit(seli,seli)*hessianfit(selj,selj));
    end;

     if (corrij < -1) | (corrij > 1);
        hessianfit(seli,selj)=0;
        errorij = [ errorij [seli selj corrij] ];
    end;   
     hessianfit(selj,seli) = hessianfit(seli,selj);

%     locate 6,1;
%     "Value Used: " hessianfit[seli,selj];
%     "Correlation:" corrij;
%     "Number of Errors:" rows(errorij)-1;
%     selj=selj+1;
 end;
%   seli = seli+1;
end;

%cls;
disp('Errors')
disp(errorij);


%/*******************************************************************************

bbbb=xparam1;
%  func =fval;
%  grad=grad;
%  retcode=exitflag
  opfhessfit = (-hessianfit);
  invhess=inv(opfhessfit);
  stdh=sqrt(diag(invhess));
  pr =length(xparam1);
  tstath=zeros(pr,1);
  i = 1; while i <= pr ; %do until i>pr;
    tstath(i)=abs(bbbb(i))/stdh(i);
  i=i+1; end ; %endo;
%tstath
%  print "t-stats. from the Hessian";
  disp('print "t-stats. from the Hessian" ') ;
 disp([xparam1 stdh tstath]);   

bbbb=xparam1;
%  func =fval;
%  grad=grad;
%  retcode=exitflag
 % opfhessfit = (-hessianfit);
 % invhess=inv(opfhessfit*.5+opfhess*.5);
 % stdh=sqrt(diag(invhess));
 % pr =length(xparam1);
 % tstath=zeros(pr,1);
 % i = 1; while i <= pr ; %do until i>pr;
 %   tstath(i)=abs(bbbb(i))/stdh(i);
 % i=i+1; end ; %endo;
%tstath
%  print "t-stats. from the Hessian";
 % disp('print "t-stats. from the Hessian" ') ;
% disp([xparam1 stdh tstath]);   

%opfhessfit = -hessianfit*.5+opfhess*.5;
hessian=opfhessfit;



return;




























npara = length(xparam1);
para = xparam1;
ndx = 1;%6
%dx =  exp(-seqa(6,2,ndx));
%dx =  exp([-6:-2:-16]);
dx =  exp([-10]);
hessianfit = zeros( npara, npara );
gradx = zeros(ndx,1);
grady = zeros(ndx,1);
gradxy = zeros(ndx,1);
hessdiag = zeros(ndx,1);
dxscale = ones(npara,1);


%/* Compute Diagonal elements first
%*/
seli = 1;
%do until seli > npara;
while seli <= npara;
%     locate 1,1;
%     "hessianfit Element    (" seli seli ")";
     i=1;
     while i <= ndx;
      paradx = para;
      parady = para;
      paradx(seli) = paradx(seli) + dx(i)*dxscale(seli);
      parady(seli) = parady(seli) - dx(i)*dxscale(seli);
      paradxdy = paradx;
      paradxdy(seli) = paradxdy(seli) - dx(i)*dxscale(seli);
      fx  = optmumlik20(para,q,h2,aa,mf,gitno,gst,gobs,gend,gstart,nlags,rawdata,lpd);
      fdx  = optmumlik20(paradx,q,h2,aa,mf,gitno,gst,gobs,gend,gstart,nlags,rawdata,lpd);
      fdy  = optmumlik20(parady,q,h2,aa,mf,gitno,gst,gobs,gend,gstart,nlags,rawdata,lpd);
      fdxdy  = optmumlik20(paradxdy,q,h2,aa,mf,gitno,gst,gobs,gend,gstart,nlags,rawdata,lpd);
      gradx(i) = -( fx - fdx )/ (dx(i)*dxscale(seli));
      grady(i) = ( fx - fdy )/ (dx(i)*dxscale(seli));
      gradxy(i) = -(fx -fdxdy)/ sqrt( (dx(i)*dxscale(seli))^2 + (dx(i)*dxscale(seli))^2 );
      hessdiag(i) = -( 2*fx - fdx - fdy)/(dx(i)*dxscale(seli))^2; 
      hessdiag(i) = -( fx - fdx - fdy + fdxdy )/(dx(i)*dx(i)*dxscale(seli)*dxscale(seli));
      i = i+1;
  end;
%     "Values";
%     -hessdiag';
     hessianfit(seli,seli) = -1*(hessdiag(1));
%     hessianfit(seli,seli) = -0.5*(hessdiag(3)+hessdiag(4));
%     locate 6,1;
%     "Value Used:" hessianfit[seli,seli];
   seli = seli+1;
end;

%/* Now compute off-diagonal elements
%** Make sure that correlations are between -1 and 1
%** errorij contains the index of elements that are invalid
%*/
errorij = [ 0 0 0];

seli = 1;
for seli = 1:npara;
   selj = seli+1
   while selj <= npara;
%     locate 1,1;
%     "hessianfit Element    (" seli selj ")";
     i=1;
     while i <= ndx;
      paradx = para;
      parady = para;
      paradx(seli) = paradx(seli) + dx(i)*dxscale(seli);
      parady(selj) = parady(selj) - dx(i)*dxscale(selj);
      paradxdy = paradx;
      paradxdy(selj) = paradxdy(selj) - dx(i)*dxscale(selj);
      fx  = optmumlik20(para,q,h2,aa,mf,gitno,gst,gobs,gend,gstart,nlags,rawdata,lpd);
      fdx  = optmumlik20(paradx,q,h2,aa,mf,gitno,gst,gobs,gend,gstart,nlags,rawdata,lpd);
      fdy  = optmumlik20(parady,q,h2,aa,mf,gitno,gst,gobs,gend,gstart,nlags,rawdata,lpd);
      fdxdy  = optmumlik20(paradxdy,q,h2,aa,mf,gitno,gst,gobs,gend,gstart,nlags,rawdata,lpd);
      gradx(i) = -( fx - fdx )/ (dx(i)*dxscale(seli));
      grady(i) = ( fx - fdy )/ (dx(i)*dxscale(selj));
      gradxy(i) = -(fx -fdxdy)/ sqrt( (dx(i)*dxscale(selj))^2 + (dx(i)*dxscale(seli))^2 );
      hessdiag(i) = -( 2*fx - fdx - fdy)/(dx(i)*dxscale(seli))^2; 
      hessdiag(i) = -( fx - fdx - fdy + fdxdy )/(dx(i)*dx(i)*dxscale(seli)*dxscale(selj));
      i = i+1
  end;
%     "Values";
%     -hessdiag';

%     hessianfit(seli,selj) = -0.5*(hessdiag(3)+hessdiag(4));
     hessianfit(seli,selj) = -1*(hessdiag(1));
     
     if ( hessianfit(seli,seli) == 0 ) or (hessianfit(selj,selj) == 0);
        corrij = 0;
     else;
        corrij = hessianfit(seli,selj)/sqrt(hessianfit(seli,seli)*hessianfit(selj,selj));
    end;

     if (corrij < -1) or (corrij > 1);
        hessianfit(seli,selj)=0;
        errorij = [ errorij [seli selj corrij] ];
    end;   
     hessianfit(selj,seli) = hessianfit(seli,selj);

%     locate 6,1;
%     "Value Used: " hessianfit[seli,selj];
%     "Correlation:" corrij;
%     "Number of Errors:" rows(errorij)-1;
     selj=selj+1
 end;
%   seli = seli+1;
end;

%cls;
disp('Errors')
disp(errorij);


%/*******************************************************************************























%new;
%closeall;
%library user, pgraph, lbdlib;
%format /mb1 /ros 16,8;
%cls;

%/******************************************************** 
%**      Estimate the DSGE Model
%**      Models: 1 = RBC
%**              2 = LBD
%**              3 = LBD + Effort
 
%mspec  = 3;
%mprior = 2;
%npara  = 12;

%/******************************************************** 
%** Import data on output growth and inflation: series (nobs,2)
%** observations from 1954:III to 1997:IV: 
%**
%** YY is composed of gdpq_cld and blsh_cl
%*/
%#include c:\projects\active\persist\Gauss\prog_t03\loaddata.g
%loadm path = ^lpath para_names;
%
%/******************************************************** 
%** Load Posterior Mode Estimates

%lpara = lpath $+ "\\" $+ lmodel $+ lprior $+ "mode";
%open fhpara = ^lpara for read;
%para = readr(fhpara,npara);
%closeall fhpara;

%"Parameter   | Estimate ";
%outmat = para_names'~para;
%let mask[1,2] = 0 1;
%%let fmt[2,3] =
%   "-*.*s " 10 4
%   "*.*lf " 10 4;
%d = printfm(outmat,(0 ~ 1),fmt);

%"";
%"Prior*Likelihood at Mode";
%fcn(para);
%"Press Key to Continue";
%k = keyw;
%cls;
% $$$  
% $$$ /*
% $$$ goto evalhess;
% $$$ */
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ /* Initialize Output files
% $$$ */
% $$$ ohess = lpath $+ "\\" $+ lmodel $+ lprior $+ "hess";
% $$$ create fhhess=^ohess with hess, npara, 8;
% $$$ wr    = writer(fhhess,hessianfit[1:npara,1:npara]);
% $$$ closeall fhhess;
% $$$    
% $$$ 
% $$$ 
% $$$ 
% $$$ /* Load hessianfit, compute penalty
% $$$ */
% $$$ evalhess:
% $$$ 
% $$$ lhess = lpath $+ "\\" $+ lmodel $+ lprior $+ "hess";
% $$$ open fhhess = ^lhess for read;
% $$$ HHm   = readr( fhhess,npara );
% $$$ closeall fhhess;
% $$$ 
% $$$ /* hessianfit is of reduced rank. Keep zero rows and columns and do SVD
% $$$ */
% $$$ if mspec == 1;
% $$$    rankHHm = 9;
% $$$ elseif mspec == 2;
% $$$    rankHHm = 11;
% $$$ else;
% $$$    rankHHm = 12;
% $$$ endif;
% $$$ 
% $$$ /* Create Inverse by Singular Value Decomposition
% $$$ */
% $$$ {u , s, v} = svd1(HHM);
% $$$ invHHMdet = 1;
% $$$ 
% $$$ i = 1;
% $$$ do until i > npara;
% $$$    if i > rankHHM;
% $$$       s[i,i] = 0;
% $$$    else;
% $$$       s[i,i]    = 1/s[i,i];
% $$$       invHHMdet = invHHMdet*s[i,i];
% $$$    endif;
% $$$    i = i+1;
% $$$ endo;
% $$$ 
% $$$ invHHM  = u*s*u';
% $$$ sigmult = u*sqrt(s);
% $$$ 
% $$$ "Determinant of minus hessianfit";
% $$$ invHHMdet;
% $$$ "sqrt(Diagonal of Inverse hessianfit)";
% $$$ sqrt(diag(invHHM));
% $$$ 
% $$$ "Post Mode Penalty";
% $$$ penalt = (rankHHM/2)*ln(2*pi) + 0.5*ln(invHHMdet);
% $$$ penalt;
% $$$ 
% $$$ /* Initialize Output files
% $$$ */
% $$$ omult = lpath $+ "\\" $+ lmodel $+ lprior $+ "mult";
% $$$ create fhmult=^omult with MULT, npara, 8;
% $$$ wr = writer(fhmult,sigmult);
% $$$ 
% $$$ closeall fhmult;
% $$$ end;
% $$$ 
% $$$ 
% $$$ %/****************************************************/
% $$$ %/*                 PROCEDURES                       */
% $$$ %/****************************************************/
% $$$ 
% $$$ 
% $$$ %proc (1) = fcn(para);
% $$$ %local lnpY, lnprio1, lnprio2, obsmean, obsvar;
% $$$ 
% $$$ %{lnpy, obsmean, obsvar} = evallbd( para,mspec,T0,YY);
% $$$ 
% $$$ %/* Evaluate the Prior density
% $$$       
% $$$ %  lnprio1 = priodens( para, pmean, pstdd, pshape);
% $$$ %  lnprio2 = priomuphi( para );
% $$$ 
% $$$ %retp(real(lnpY+lnprio1+lnprio2));  
% $$$ %endp;
% $$$ 
% $$$ /***************************************************************************
% $$$ */
% $$$ 
% $$$ proc (1) = priomuphi(para);
% $$$ local muphi, lnprio;
% $$$ muphi = para[7:8];
% $$$ if mspec > 1;
% $$$    lnprio = -ln(2*pi) - 0.5*ln(det(muphi_v0))
% $$$             - 0.5*(muphi - muphi_m0)'*inv(muphi_v0)*(muphi - muphi_m0);
% $$$ else;
% $$$    lnprio = 0;
% $$$ endif;
% $$$ retp(lnprio);
% $$$ endp;
% $$$ 
% $$$ 
% $$$ 
