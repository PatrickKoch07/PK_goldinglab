function xpand_segimage(p,imnamenum,amount);
% function X_expand_segimage(imnamenum,amount);
% Incase the segmented image cuts off some cells. 
% Resaves the images of frame number imnamenum
% with bigger size expanded by amount.

% NOTES on variable name changes
%   tifdir    -> p.imageDir
%   imprefix  -> p.movieName
%   regsize   -> p.regsize
%   outprefix -> p.outprefix

mynum = str3(imnamenum);
if p.numphaseslices == 1
    pname = [p.imageDir,p.movieName,'-p-',mynum,'.tif'];
else
    pname = [p.imageDir,p.movieName,'-p-2-',mynum,'.tif'];
end
phim = imread(pname);
Lname= [p.outprefix,mynum];   
cname= [p.imageDir,p.movieName,'-c-',mynum,'.tif'];
yname= [p.imageDir,p.movieName,'-y-',mynum,'.tif'];
gname= [p.imageDir,p.movieName,'-g-',mynum,'.tif'];
rname= [p.imageDir,p.movieName,'-r-',mynum,'.tif'];

load([p.segmentationDir,Lname])
newrect = [ max(1+p.regsize,rect(1)-amount), ...
	    max(1+p.regsize,rect(2)-amount), ...
	    min(size(phim,1)-p.regsize,rect(3)+amount), ...
	    min(size(phim,2)-p.regsize,rect(4)+amount)      ];
newLNsub = zeros(size(phim));
newLNsub(rect(1):rect(3), rect(2):rect(4)) = LNsub;
rect = newrect;
LNsub = newLNsub(rect(1):rect(3), rect(2):rect(4));
phsub = phim(rect(1):rect(3), rect(2):rect(4));
savelist=['''phsub'',''LNsub'',''rect'''];
if exist(cname)==2
    disp('found CFP image');
    [creg, cshift, cback]= quicknoreg(LNsub, cname, rect, p.regsize, p.fullsize);
    [exptcstr, gainc, exptc]= imsettings(cname);
    savelist=[savelist,',''creg'',''cshift'',''exptc'',''gainc'',''cback'''];
end 
if exist(yname)==2
    disp('found YFP image');

    [yreg, yshift, yback]= quicknoreg(LNsub, yname, rect, p.regsize, p.fullsize);
    [exptystr, gainy, expty]= imsettings(yname);
    savelist=[savelist,',''yreg'',''yshift'',''expty'',''gainy'',''yback'''];
end 
if exist(gname)==2
    disp('found GFP image');
    [greg, gshift, gback]= quicknoreg(LNsub, gname, rect, p.regsize, p.fullsize);
    [exptystr, gaing, exptg]= imsettings(gname);
    savelist=[savelist,',''greg'',''gshift'',''exptg'',''gaing'',''gback'''];
end 
if exist(rname)==2
    disp('found RFP image');
    [greg, gshift, gback]= quicknoreg(LNsub, rname, rect, p.regsize, p.fullsize);
    [exptystr, gaing, exptg]= imsettings(gname);
    savelist=[savelist,',''rreg'',''rshift'',''exptr'',''gainr'',''rback'''];
end 

eval(['save(''',p.segmentationDir,Lname,''',',savelist,');']);
disp(['saved file ',p.segmentationDir,Lname]);
