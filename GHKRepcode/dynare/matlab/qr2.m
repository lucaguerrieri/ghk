function [Q,R] = qr2(X)
% stephane.adjemian@ens.fr [12-07-2005]
%
% This routine performs a qr decomposition of matrix X such that the 
% diagonal scalars of the upper-triangular matrix R are positive.
[Q,R] = qr(X);
indx = find(diag(R)<0);
if ~isempty(indx)
    Q(:,indx) = -Q(:,indx);
    R(indx,:) = -R(indx,:);
end