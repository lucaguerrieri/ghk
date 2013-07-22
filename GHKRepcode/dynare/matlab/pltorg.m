function [nbplt,nr,nc,lr,lc,nstar] = pltorg(number)
% stephane.adjemian@cepremap.cnrs.fr [06-07-2004]
nrstar = 3;
ncstar = 3;
nstar  = nrstar*ncstar;
nbplt  = 0;
nr     = 0;
nc     = 0;
lr     = 0;
lc     = 0;
if number == 1
    nbplt = 1;
    nr    = 1;
    nc    = 1;
    nstar = number;
elseif number == 2
    nbplt = 1;
    nr    = 2;
    nc    = 1;
    nstar = number;
elseif number == 3
    nbplt = 1;
    nr    = 3;
    nc    = 1;
    nstar = number;
elseif number == 4
    nbplt = 1;
    nr    = 2;
    nc    = 2;
    nstar = number;
elseif number == 5
    nbplt = 1;
    nr    = 3;
    nc    = 2;
    nstar = number;
elseif number == 6
    nbplt = 1;
    nr    = 3;
    nc    = 2;
    nstar = number;    
elseif number == 7
    nbplt = 1;
    nr    = 3;
    nc    = 3;
    nstar = number;    
elseif number == 8
    nbplt = 1;
    nr    = 3;
    nc    = 3;
    nstar = number;
elseif number == 9
    nbplt = 1;
    nr    = 3;
    nc    = 3;
else
    if number/nstar == round(number/nstar)
        nbplt = number/nstar;
        nr    = nrstar;
        nc    = ncstar;
        lr    = nr;
        lc    = nc; 
    else
        nbplt = ceil(number/nstar);
        nr    = nrstar;
        nc    = ncstar;
        reste = number-(nbplt-1)*nstar;
        if reste == 1
            lr    = 1;
            lc    = 1;
        elseif reste == 2
            lr    = 2;
            lc    = 1;
        elseif reste == 3
            lr    = 3;
            lc    = 1;
        elseif reste == 4
            lr    = 2;
            lc    = 2;
        elseif reste == 5
            lr    = 3;
            lc    = 2;
        elseif reste == 6
            lr    = 3;
            lc    = 2;    
        elseif reste == 7
            lr    = 3;
            lc    = 3;    
        elseif reste == 8
            lr    = 3;
            lc    = 3;
        end
    end
end
