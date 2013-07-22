function printpref

% sets preferences for printing

% chooses landscape orientation
% and enlarges charts to fill page
orient('portrait');

set(gcf, 'PaperUnits', 'inches');

% gcf returns handle to current active figure
papersize = get(gcf, 'PaperSize');
width = .85*papersize(1);         % Initialize a variable for width.
height = .85*papersize(2);          % Initialize a varible for height.
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf, 'PaperPosition', myfiguresize);