function [nam,texnam] = get_the_name(k,TeX)
% stephane.adjemian@cepremap.cnrs.fr [07-13-2004]
global lgx_ lgy_ lgx_TeX_ estim_params_ options_

nam = [];
texnam = [];

if k <= estim_params_.nvx
  vname = deblank(lgx_(estim_params_.var_exo(k,1),:));
  nam=['SE_',vname];
  if TeX
    tname  = deblank(lgx_TeX_(estim_params_.var_exo(k,1),:));
    % tname = vname;
    texnam = ['$ SE_{' tname '} $'];
  end
elseif  k <= (estim_params_.nvx+estim_params_.nvn)
  vname = deblank(options_.varobs(estim_params_.var_endo(k-estim_params_.nvx,1),:));
  nam=['SE_EOBS_',vname];
  if TeX
    tname  = deblank(options_.TeX_varobs(estim_params_.var_endo(k-estim_params_.nvx,1),:));
    % tname = vname;
    texnam = ['$ SE_{' tname '} $'];
  end
elseif  k <= (estim_params_.nvx+estim_params_.nvn+estim_params_.ncx)
  jj = k - (estim_params_.nvx+estim_params_.nvn);
  k1 = estim_params_.corrx(jj,1);
  k2 = estim_params_.corrx(jj,2);
  vname = [deblank(lgx_(k1,:)) '_' deblank(lgx_(k2,:))];
  nam=['CC_',vname];
  if TeX
    tname  = [deblank(lgx_TeX_(k1,:)) ',' deblank(lgx_TeX_(k2,:))];
    % tname = vname;
    texnam = ['$ CC_{' tname '} $'];
  end
elseif  k <= (estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+ ...
              estim_params_.ncn)
  jj = k - (estim_params_.nvx+estim_params_.nvn+estim_params_.ncx);
  k1 = estim_params_.corrn(jj,1);
  k2 = estim_params_.corrn(jj,2);
  vname = [deblank(lgy_(k1,:)) '_' deblank(lgy_(k2,:))];
  nam=['CC_EOBS_' vname];
  if TeX
    tname  = [deblank(lgy_TeX_(k1,:)) ',' deblank(lgy_TeX_(k2,:))];
    % tname = vname;
    texnam =['$ CC_{' tname '} $'];
  end
else
  jj = k - (estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+estim_params_.ncn); 
  nam = deblank(estim_params_.param_names(jj,:));
  if TeX
    texnam = ['$ '  deblank(estim_params_.tex(jj,:))  ' $'];
    % texnam = nam;
  end    
end


% SA 07-15-2004 Added TeX names.
% SA 12-02-2004 Changed non-TeX names format.
