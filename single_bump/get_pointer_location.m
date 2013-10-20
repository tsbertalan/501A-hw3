% [x, y] = get_pointer_location(ax)
%
% Given a Matlab graphics handle to existing Axes, returns the x,y location
% of the current position of the mouse in the units of the passed axis.
%

function [x, y] = get_pointer_location(ax)

pt = get(0, 'PointerLocation');

fig = get(ax, 'Parent'); 
set(fig, 'Units', 'pixels');
figpos = get(fig, 'Position');

set(ax, 'Units', 'pixels'); axpos  = get(ax, 'Position'); 

xl = get(ax, 'XLim');
x = round(xl(1) + diff(xl) * (pt(1) - figpos(1) - axpos(1))/axpos(3));

yl = get(ax, 'YLim');
y = round(yl(1) + diff(yl) * (pt(2) - figpos(2) - axpos(2))/axpos(4));

