% Copyright (C) 2001 Michel Juillard
%
% apply fast Givens rotation to a pair of rows

function [a] = rfrot(a, alph, bet, type)
     if type == 1
         tau = bet*a(1,:) + a(2,:);
         a(2,:) = a(1,:) + alph*a(2,:);
         a(1,:) = tau;
     else
         tau = a(1,:) + bet*a(2,:);
         a(2,:) = alph*a(1,:) + a(2,:);
         a(1,:) = tau;
     end
