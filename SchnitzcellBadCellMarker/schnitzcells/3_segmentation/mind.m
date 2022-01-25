<<<<<<< HEAD

function dist = mind(pt_x, pt_y, img, cellno);
% function dist = mind(pt_x, pt_y, img, cellno);
%
% returns minimum distance of (pt_x, pt_y) from cell denoted by cellno
%  in image img

[fy,fx] = find(img==cellno);
dy = fy - pt_y;
dx = fx - pt_x;
d = sqrt(dx.^2 + dy.^2);
dist = min(d);

=======

function dist = mind(pt_x, pt_y, img, cellno);
% function dist = mind(pt_x, pt_y, img, cellno);
%
% returns minimum distance of (pt_x, pt_y) from cell denoted by cellno
%  in image img

[fy,fx] = find(img==cellno);
dy = fy - pt_y;
dx = fx - pt_x;
d = sqrt(dx.^2 + dy.^2);
dist = min(d);

>>>>>>> 88dd4f4ac7eb6f5ae758552b04c47df241263be0
