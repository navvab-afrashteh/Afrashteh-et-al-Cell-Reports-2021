function CranialWinOrient(midX, midY, hL, sp, orientation, FontSize, col)

set(gca,'unit','inches');
pos = get(gca,'pos'); r = min(pos(3),pos(4))/2;
r = sqrt(r);
linewidth = 2*r;

if ~nargin
    % cranial window orientation
    hL = 12*r; % half length of the cross
    hL = 2*round(hL/2);
    sp = 2*r; % space between letters and the ends of cross
    FontSize = 12*r;
    dim1 = 256; dim2 = 256; a = 0; b = 25; SB = 7.8;
    midX = dim2-SB/2-a;
    midY = dim1-b-hL-sp-FontSize*2/r;
    col = 'r';
    orientation.top = 'A'; orientation.bottom = 'P'; orientation.left = 'M'; orientation.right = 'L';
end

hold on;
plot(midX+[-hL hL], midY*[1 1], 'color', col,'linewidth',linewidth)
plot(midX*[1 1], midY+[-hL hL], 'color', col,'linewidth',linewidth)
text(midX,midY+hL+sp*1.33,orientation.top,'color',col,'FontSize',FontSize,'HorizontalAlignment','center','VerticalAlignment','baseline');
text(midX,midY-hL-sp*1.33,orientation.bottom,'color',col,'FontSize',FontSize,'HorizontalAlignment','center','VerticalAlignment','cap');
text(midX+hL+sp*0.67,midY,orientation.right,'color',col,'FontSize',FontSize,'HorizontalAlignment','left','VerticalAlignment','middle');
text(midX-hL-sp*1.33,midY,orientation.left,'color',col,'FontSize',FontSize,'HorizontalAlignment','right','VerticalAlignment','middle');
n=0;
