function my_subplot(i,imax,irow,icol,fig_title)
% spreads subplots on several figures according to a maximum number of
% subplots per figure
%
% INPUTS
% i: subplot number
% imax: total number of subplots
% irow: maximum number of rows in a figure
% icol: maximum number of columns in a figure
% fig_title: title to be repeated on each figure
  nfig_max = irow*icol;
  if imax < nfig_max
    icol = ceil(sqrt(imax));
    irow=icol;
    if (icol-1)*(icol-2) >= imax
      irow = icol-2;
      icol = icol-1;
    elseif (icol)*(icol-2) >= imax
      irow = icol-2;
    elseif icol*(icol-1) >= imax
      irow = icol-1;
    end
  end

  i1 = mod(i-1,nfig_max);
  if i1 == 0
    figure('Name',fig_title);
  end
  
  subplot(irow,icol,i1+1);