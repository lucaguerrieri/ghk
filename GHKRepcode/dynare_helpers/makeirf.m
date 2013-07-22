%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get IRFS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
history=mymkirf(dr_,ys_,nperiods,shock,0);  % returns matrix containing IRFs for all variables

reordered = lgy_(dr_.order_var,:);          % the vector reordered contains names of variables
                                            % corresponding to rows of
                                            % history

% this loop assigns each row of history to its variable
% it also assigns steady state values
% IRFs are stored under varname_irf  where varname is the variable name
% SS values are stored under varname
% NB: SS values are scalar, IRFs are column vectors.
nvars = size(reordered,1);
for indxi = 1:nvars
    eval([dset,'_',deblank(reordered(indxi,:)),'_irf=transpose(history(indxi,:));']);
    eval([dset,'_',deblank(lgy_(indxi,:)),'= ys_(indxi);']);
end




