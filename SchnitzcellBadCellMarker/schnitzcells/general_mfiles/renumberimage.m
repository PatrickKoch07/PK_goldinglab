<<<<<<< HEAD

function imnew = renumberimage(imold);

imnew = imold;
u = unique(imold(:));
% pk 11/26/2021, only had u_p before
u_p = u(u >= 0);
u_n = u(u <= 0);
%
for i = 2:length(u_p),
    imnew(imold==u_p(i))=i-1;
end;
for i = 1:(length(u_n)-1),
    imnew(imold==u_n(i))=i*-1;
end;

if islogical(imnew),
    imnew = +imnew;
end;
=======

function imnew = renumberimage(imold);

imnew = imold;
u = unique(imold(:));
% pk 11/26/2021, only had u_p before
u_p = u(u >= 0);
u_n = u(u <= 0);
%
for i = 2:length(u_p),
    imnew(imold==u_p(i))=i-1;
end;
for i = 1:(length(u_n)-1),
    imnew(imold==u_n(i))=i*-1;
end;

if islogical(imnew),
    imnew = +imnew;
end;
>>>>>>> 88dd4f4ac7eb6f5ae758552b04c47df241263be0
