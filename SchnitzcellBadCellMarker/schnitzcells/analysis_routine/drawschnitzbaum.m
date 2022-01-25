<<<<<<< HEAD

function drawschnitzbaum(schnitzcells,me,x,level);

% draws a tree from a schnitzcell structure, starting with "me" as the
% root.
% when you run this yourself, omit "x" and "level" (those parameters are 
% used by drawschnitzbaum when it calls itself recursively.

miny = 0;
minc = 0;

ccolor = [1 1 1];
ycolor = [1 1 1];

cedgecolor = 0.5*[0 1 0];
yedgecolor = 0.5*[1 0 0];

% maxy = 3/2;
% maxc = 211/2;

maxc = 0.9*max([schnitzcells.MC]);
maxy = 0.9*max([schnitzcells.MY]);

if nargin == 2,
    x = 0;
    level = 1;
    figure;
end;

epsilon = 0.01;
epsilon = 0.02;

% first draw ourself
frames = schnitzcells(me).frames;
plot(x*ones(length(frames),1),frames,'.-');
hold on;
c = schnitzcells(me).MC;
y = schnitzcells(me).MY;
c(c<minc) = minc;
y(y<miny) = miny;

c = (c - minc)./(maxc-minc);
y = (y - miny)./(maxy-miny);

c(c>1) = 1;
y(y>1) = 1;

for i = 1:length(frames),
    if ~isnan(c(i)),
        h1 = plot(x+epsilon,frames(i), 's');
        set(h1,'markerfacecolor', (c(i)*ccolor));
        set(h1,'markeredgecolor', cedgecolor);
        set(h1,'markersize', 10);
    end;
    if ~isnan(y(i)),
        h2 = plot(x-epsilon,frames(i), 's');
        set(h2,'markerfacecolor', (y(i)*ycolor));
        set(h2,'markeredgecolor', yedgecolor);
        set(h2,'markersize', 10);
    end;
end;
% text((x+0.05)*ones(length(frames),1),frames,num2str(schnitzcells(me).cellno'),'color','r','fontsize',7);
%  text((x+0.05)*ones(length(frames),1),frames,num2str(schnitzcells(me).ancestor'),'color','r','fontsize',7);
% now draw our daughters:
if schnitzcells(me).D>0,
    drawschnitzbaum(schnitzcells,schnitzcells(me).D,x+1/2^level,level+1);
    line([x x+1/2^level],[frames(end) frames(end)+1])
end;
if schnitzcells(me).E>0,
    drawschnitzbaum(schnitzcells,schnitzcells(me).E,x-1/2^level,level+1);
    line([x x-1/2^level],[frames(end) frames(end)+1])
end;

if nargin == 2,
%     set(gca,'color','k');
    a = axis;
    a(1:2) = [-1 1];
    axis(a);
end;
=======

function drawschnitzbaum(schnitzcells,me,x,level);

% draws a tree from a schnitzcell structure, starting with "me" as the
% root.
% when you run this yourself, omit "x" and "level" (those parameters are 
% used by drawschnitzbaum when it calls itself recursively.

miny = 0;
minc = 0;

ccolor = [1 1 1];
ycolor = [1 1 1];

cedgecolor = 0.5*[0 1 0];
yedgecolor = 0.5*[1 0 0];

% maxy = 3/2;
% maxc = 211/2;

maxc = 0.9*max([schnitzcells.MC]);
maxy = 0.9*max([schnitzcells.MY]);

if nargin == 2,
    x = 0;
    level = 1;
    figure;
end;

epsilon = 0.01;
epsilon = 0.02;

% first draw ourself
frames = schnitzcells(me).frames;
plot(x*ones(length(frames),1),frames,'.-');
hold on;
c = schnitzcells(me).MC;
y = schnitzcells(me).MY;
c(c<minc) = minc;
y(y<miny) = miny;

c = (c - minc)./(maxc-minc);
y = (y - miny)./(maxy-miny);

c(c>1) = 1;
y(y>1) = 1;

for i = 1:length(frames),
    if ~isnan(c(i)),
        h1 = plot(x+epsilon,frames(i), 's');
        set(h1,'markerfacecolor', (c(i)*ccolor));
        set(h1,'markeredgecolor', cedgecolor);
        set(h1,'markersize', 10);
    end;
    if ~isnan(y(i)),
        h2 = plot(x-epsilon,frames(i), 's');
        set(h2,'markerfacecolor', (y(i)*ycolor));
        set(h2,'markeredgecolor', yedgecolor);
        set(h2,'markersize', 10);
    end;
end;
% text((x+0.05)*ones(length(frames),1),frames,num2str(schnitzcells(me).cellno'),'color','r','fontsize',7);
%  text((x+0.05)*ones(length(frames),1),frames,num2str(schnitzcells(me).ancestor'),'color','r','fontsize',7);
% now draw our daughters:
if schnitzcells(me).D>0,
    drawschnitzbaum(schnitzcells,schnitzcells(me).D,x+1/2^level,level+1);
    line([x x+1/2^level],[frames(end) frames(end)+1])
end;
if schnitzcells(me).E>0,
    drawschnitzbaum(schnitzcells,schnitzcells(me).E,x-1/2^level,level+1);
    line([x x-1/2^level],[frames(end) frames(end)+1])
end;

if nargin == 2,
%     set(gca,'color','k');
    a = axis;
    a(1:2) = [-1 1];
    axis(a);
end;
>>>>>>> 88dd4f4ac7eb6f5ae758552b04c47df241263be0
