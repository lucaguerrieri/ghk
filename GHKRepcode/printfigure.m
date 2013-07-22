function printfigure(figurename,nfigures)

for i=1:nfigures
    figure(i)
    if (i==1)
       eval(['print -dpsc2 ',figurename,'.ps']);
    else
       eval(['print -dpsc2 -append ',figurename,'.ps']);
    end
end