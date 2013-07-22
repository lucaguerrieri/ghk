% Program calculating the non-l;inear constraint for posterior mximisation
% M. Ratto
% adapted from mj_optmumlik.m
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
function [c, ceq] = mr_nlincon(xparam1,gend,rawdata,algo);


global bayestopt_ exo_nbr dr_ estim_params_ Sigma_e_ options_ xparam_test ...
    trend_coeff_ 


c=0;
ceq=0;

xparam_test = xparam1;
cost_flag = 1;
if options_.mode_compute ~= 1 && any(xparam1 < bayestopt_.lb)
  k = find(xparam1 < bayestopt_.lb);
  c = sum(bayestopt_.lb(k)-xparam1(k));
  cost_flag = 0;
  return;
end
if options_.mode_compute ~= 1 && any(xparam1 > bayestopt_.ub)
  k = find(xparam1 > bayestopt_.ub);
  c = sum(xparam1(k)-bayestopt_.ub(k));
  cost_flag = 0;
  return;
end

nobs = size(options_.varobs,1);

Q = Sigma_e_;
for i=1:estim_params_.nvx
	k =estim_params_.var_exo(i,1);
	Q(k,k) = xparam1(i)*xparam1(i);
end
offset = estim_params_.nvx;
if estim_params_.nvn
	H = zeros(nobs,nobs);
	for i=1:estim_params_.nvn
		k = estim_params_.var_endo(i,1);
		H(k,k) = xparam1(i+offset)*xparam1(i+offset);
	end
	offset = offset+estim_params_.nvn;
end	
if estim_params_.ncx
	for i=1:estim_params_.ncx
		k1 =estim_params_.corrx(i,1);
		k2 =estim_params_.corrx(i,2);
		Q(k1,k2) = xparam1(i+offset)*sqrt(Q(k1,k1)*Q(k2,k2));
		Q(k2,k1) = Q(k1,k2);
	end
	[CholQ,testQ] = chol(Q);
	if testQ 	%% The variance-covariance matrix of the structural innovations is not definite positive.
    			%% We have to compute the eigenvalues of this matrix in order to build the penalty.
    	a = diag(eig(Q));
		%fval = bayestopt_.penalty*min(1e3,exp(sum(-a(a<=0))));
    	c =  sum(-a(a<=0));
    	cost_flag = 0;
		return
	end
	offset = offset+estim_params_.ncx;
end
if estim_params_.nvn & estim_params_.ncn 
	for i=1:estim_params_.ncn
		k1 = options_.lgyidx2varobs(estim_params_.corrn(i,1));
		k2 = options_.lgyidx2varobs(estim_params_.corrn(i,2));
		H(k1,k2) = xparam1(i+offset)*sqrt(H(k1,k1)*H(k2,k2));
		H(k2,k1) = H(k1,k2);
	end
	[CholH,testH] = chol(H);
	if testH
		a = diag(eig(H));
		if nobs == estim_params_.nvn
			% fval = bayestopt_.penalty*min(1e3,exp(sum(-a(a<=0))));
   			c = sum(-a(a<=0));
   			cost_flag = 0;
			return
		else
			if sum(abs(a)<crit) == nobs-estim_params_.nvn
				if any(a<0)
					% fval = bayestopt_.penalty*min(1e3,exp(sum(-a(a<0))));
   					c = sum(-a(a<0));
   					cost_flag = 0;
					return					
				else
					% All is fine, there's nothing to do here...
				end 					
			else
				fval = bayestopt_.penalty*min(1e3,exp(sum(-a(a<=0))));
   				c = sum(-a(a<=0));
   				cost_flag = 0;
				return			
			end 
		end
	end
	offset = offset+estim_params_.ncn;
end
for i=1:estim_params_.np
	assignin('base',deblank(estim_params_.param_names(i,:)),xparam1(i+offset));
end
Sigma_e_ = Q;


%------------------------------------------------------------------------------
% 2. call model setup & reduction program
%------------------------------------------------------------------------------

[A,B,ys,info] = dynare_resolve;
if info(1) == 1 | info(1) == 2 | info(1) == 5
  fval = bayestopt_.penalty;
  cost_flag = 0;
  return
elseif info(1) == 3 | info(1) == 4 | info(1) == 20
  fval = bayestopt_.penalty*min(1e3,exp(info(2)));
  cost_flag = 0;
  return
end



% 11/18/03 MJ changed input parameters for priordens()
