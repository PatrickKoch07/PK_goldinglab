function p = manualcheckseg (p, varargin);
% function P = manualcheckseg (p, varargin);
%
%   MANUALCHECKSEG allows users to review the results of image segmentation
%   in order to manually correct any merged cells or delete false positives.
%
%   MANUALCHECKSEG(P,'Field1',Value1,'Field2',Value2,...) also performs
%   manual/interactive correction of segmentation results, but permits users
%   to adjust any parameters describing the movie or parameters controlling
%   the manual segmentation checking process.  The routine first sets all of
%   the movie analysis parameters to their default values defined for the
%   given
%   movieKind, and then overrides any specific parameters provided by setting
%   P.Field1 = Value1, P.Field2 = Value2, etc.  Thus any/all schnitzcells
%   parameter values can be defined in the function call via these optional
%   field/value pairs.  (This is in the style of setting MATLAB properties
%   using optional property/value pairs.)
%
%   MANUALCHECKSEG returns a struct (1x1 struct array) referred to as the
%   schnitzcells parameter structure that contains fields and values
%   contained in P, including unchanged/original parameters plus any of those
%   added or overridden in the list of properties & values provided in the
%   optional arguments.
%
%-------------------------------------------------------------------------------
% Schnitzcells Fields / Parameters that you can adjust to control manualcheckseg
%
%   outprefix     overrides [p.movieName 'seg']
%   manualRange   specify frame range to check
%   override      if 1, program will redo frames already having Lc; default=0
%   expandvalue   increases each segmentation image border by this amount
%                 (in case the segmentation cut off some cells); defaults to 30
%   frnum         check one frame only, also sets p.override = 1
%   Dskip         (what does this do? probably frame skipping?), defaults to 1
%   upend         frame size / location; default is 680
%   leftend       frame size / location; default is 516
%   min_size      [height width] of figure; default [upend-60 leftend-10]
%   finetuneimage if 1, program will not renumber the image; default = 0
%   regsize       maximum size in pixels of any translation between phase and
%                 fluorescent images; default is 3
%-------------------------------------------------------------------------------
%

%-------------------------------------------------------------------------------
% variable renaming during rewrite:
%   outprefix -> p.outprefix
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Parse the input arguments, input error checking
%-------------------------------------------------------------------------------

numRequiredArgs = 1;
if (nargin < 1) | ...
        (mod(nargin,2) == 0) | ...
        (~isSchnitzParamStruct(p))
    errorMessage = sprintf ('%s\n%s\n%s\n',...
        'Error using ==> manualcheckseg:',...
        '    Invalid input arguments.',...
        '    Try "help manualcheckseg".');
    error(errorMessage);
end

%-------------------------------------------------------------------------------
% Override any schnitzcells parameters/defaults given optional fields/values
%-------------------------------------------------------------------------------

% Loop over pairs of optional input arguments and save the given fields/values
% to the schnitzcells parameter structure
numExtraArgs = nargin - numRequiredArgs;
if numExtraArgs > 0
    for i=1:2:(numExtraArgs-1)
        if (~isstr(varargin{i}))
            errorMessage = sprintf ('%s\n%s%s%s\n%s\n',...
                'Error using ==> manualcheckseg:',...
                '    Invalid property ', num2str(varargin{i}), ...
                ' is not (needs to be) a string.',...
                '    Try "help manualcheckseg".');
            error(errorMessage);
        end
        fieldName = schnitzfield(varargin{i});
        p.(fieldName) = varargin{i+1};
    end
end



disp(' ')
disp('Instructions: press left mouse button on two consecutive cells to join them.')
disp('              Press right mouse button to cut at that point.')
disp('              Shift + left mouse button erases current cell.')
disp('              When you''re done, press <space> to go to the next frame.')
disp('              press <escape> to redo this frame from the original autoseg file.')
disp('              to save progress during work on a frame, press:')
disp('                    ''s'' to save work to memory without writing to the file.')
disp('                    ''w'' to write a temporary partial correction to the file.')
disp('              press ''x'' to black out an area.')
disp('              press ''t'' to mark terraced area.')
disp('              press ''c'' to crop out only populated area.')
disp('              press ''r'' to renumber the cell you are pointing to.')
disp('              press ''p'' to show a square around the position, on phase image.')
disp('              press ''b'' to mark the cell you are pointing to, on phase image.')
disp('              press ''o'' to obliterate all but the cell you''re pointing to.')
disp('              press ''l'' to add frame to list of badly segmented frames.')
disp('              press ''.'' to skip frame (without saving).')
disp('              press '','' to go back a frame (without saving).')
disp('              press ''a'' to add a new cell where the mouse is pointing.')
disp('              press ''q'' to quit.')
disp(' ')
disp('              press ''e'' to expand image.')
disp('              press ''f'' for "fine-tuning" (to avoid renumbering the image).')
disp('              press ''g'' to goto indexnum = ... .')
disp('              press ''R'' to renumber all cells.')
disp(' ')
disp('              PK ADD: PRESS ''`'' TO MARK CELL AS BAD (makes dark and flips cellnum to negative).')
disp(' ')

quit_now=0;
global pos Limage ourfig res pp phfig

ourfig = figure;
phfig  = figure;

outl=length(p.segmentationDir);

% set some defaults if they don't exist
if ~existfield(p,'Dskip')
    p.Dskip=1;
end

if ~existfield(p,'outprefix')
    p.outprefix = [p.movieName 'seg'];
end

D = dir([p.segmentationDir, p.outprefix, '*.mat']);
[s,I] = sort({D.name}');
D = D(I);
numpos= findstr(D(1).name, '.mat')-3;
if ~existfield(p,'manualRange')
    % p.manualRange=1:length(dir([p.segmentationDir, p.outprefix, '*.mat']));
    segNameStrings = char(s);
    p.manualRange = str2num(segNameStrings(:,numpos:numpos+2))';
end

if existfield(p,'override')~=1
    p.override=0;
end

if ~existfield(p,'expandvalue'),
    p.expandvalue = 30;
end;

if existfield(p,'frnum')
    %  p.manualRange=1:length(dir([p.segmentationDir, p.outprefix,
    %  '*.mat']));
    newloopindex = find(p.manualRange==frnum);
    if isempty(newloopindex)
        disp(['Could not check frame ', num2str(gotoframe), ...
            ' because it is not in the manualRange of segmented files.']);
    else
        p.manualRange=frnum;
        p.override=1;
    end
end

if ~existfield(p,'upend')
    p.upend   = 680;
end
if ~existfield(p,'leftend')
    p.leftend = 516;
end
if ~existfield(p,'min_size')
    p.min_size = [p.upend-60 p.leftend-10];
end

if ~existfield(p,'finetuneimage')
    p.finetuneimage = 0;
end

if ~existfield(p,'regsize')
    p.regsize = 3;
end

backwards = 0;
gotoframenum=0;
loopindex = 1;
while loopindex <= length(p.manualRange);
    i = p.manualRange(loopindex);

    clear Lc creg yreg rreg greg savelist rect newrect oldrect;
    name= [p.segmentationDir,p.movieName,'seg',str3(i)];
    tempsegcorrect=0;
    load(name);
    % nitzan's changes June24th
    if exist('phaseFullSize')==1
        p.fullsize = phaseFullSize;
    end
    if ~isfield(p,'fullsize')
        p.fullsize = [1024 , 1344];
        disp('setting size of images arbitrarily to [1024 , 1344].');
    end
    % nitzan's changes June24th
    if (exist('Lc')~=1 | backwards | p.override==1 | (exist('Lc')==1 & tempsegcorrect==1)) & ~mod(i-1,p.Dskip)

        if exist('Lc')==1
            LNsub=Lc;
        end

        % size(g,1) is the height of the fig
        % size(g,2) is the width  of the fig
        % min_size(1) is the height of the desired figure
        % min_size(2) is the width  of the desired figure
        % pos11(3) is the width  of the actual fig
        % pos11(4) is the height of the actual fig

        %show phase image...
        g = double(phsub);
        g = (g-minmin(g)) / (maxmax(g)-minmin(g));
        g = g+0.2;
        gones = ones(size(g));
        g(g>1)= gones(g>1);
        g = uint8(g * 255);

        rescale=(p.min_size./size(LNsub));
        res=min(rescale);
        if res>2 res=2;end
        figure(phfig);
        clf reset;
        
        %--------------------------------------------%
        % mod MW Nov 26,2019
%         % show spot channel
%         cx = imread([pwd '\c3images\' p.movieName ...
%                 'seg' str3(i) 'c3' 'MaxP.tif']) ;
%         cx = cx(rect(1):rect(3), rect(2):rect(4)) ; 
%         cx = double(cx) ; 
%         I = zeros(size(cx)) ;
%         G = g ; 
%         R = double2rgb(cx, [900 1500] , 0) ; 
%         merged = cat(3, R, G, I); 
%         imshow(imresize(merged,res)) 
%         --------------------------------------------%

        imshow(imresize(g,res));
%         colormap((1:100)'*[0 1 0]/100)
        pos11=get(phfig,'position');
%         set(phfig, 'position', ...
%             [p.leftend-pos11(3)-8, p.upend-pos11(4), pos11(3), pos11(4)]);
%         set(phfig,'position',...
%             [p.leftend-pos11(3)-8,p.upend-pos11(4),1.5*pos11(3),1.5*pos11(4)]); % MW Sep 04, 2013
%         set(phfig,'position',...
%             [14 114 779 779]); % MW Jul,12,2017
         set(phfig,'position',...
             [19   175   778   778]);
%             1e3*[0.0140    0.1140    1.25    1.25]); % MW Mar,11,2020

        %    set(phfig,'name',[name(outl+1:end-4)]);
        set(phfig,'name',['Frame ',str3(i),' phase']);
        hold on;

        is_done=0;
        savelist=['''Lc'''];
        crop_pop=0;
        while ~is_done
            [Lc,is_done,quit_now,dontsave,addtolist,crop_pop,newrect,savetemp,backwards,p.finetuneimage,gotoframenum] = ...
                manual_kant(p, LNsub, p.leftend, p.upend, g, ...
                p.finetuneimage, rect, p.min_size, i, p.expandvalue);
            if backwards==1 & (loopindex>1),
                % JCR: Note this is kinda ugly- later loopindex is *incremented*
                % by Dskip, so we have to back up 2 Dskips now.
                loopindex = loopindex-2*p.Dskip;
                % disp(['backing up to loopindex = ',num2str(loopindex)]);
                % because of above ugliness must correct for 1 dskip offset here
                properloopindex = loopindex + p.Dskip;
                disp(['backing up to frame = ',...
                     str3(p.manualRange(properloopindex)),...
                     '(loopindex = ', num2str(properloopindex), ')']);
                % JCR: Note this code is still kinda bad because when 
                % loopindex = 1 and user presses backwards, loopindex skips 
                % forward anyway.  Probably should simplify and get rid of 
                % Dskip all together, now that segRange and manualRange should 
                % handle arbitrary frames
            end;
            if quit_now
                close(phfig);close(ourfig);
                return;
            end;
            if crop_pop
                LNfull=zeros(p.fullsize(1),p.fullsize(2));
                LNfull(rect(1):rect(3), rect(2):rect(4))=LNsub;
                oldrect=rect;
                rect=newrect;
                LNsub = LNfull(rect(1):rect(3), rect(2):rect(4));
                LNfull=zeros(p.fullsize(1),p.fullsize(2));
                LNfull(oldrect(1):oldrect(3), oldrect(2):oldrect(4))=phsub;
                phsub = LNfull(rect(1):rect(3), rect(2):rect(4));
%                 savelist=[savelist,',''phsub'',''LNsub'',''rect'''];
                LcFull = zeros(phaseFullSize) ; %Tommy Sep 23, 2008
                LcFull(rect(1):rect(3),rect(2):rect(4)) = Lc ; %Tommy Sep 23, 2008
                savelist=[savelist,',''phsub'',''LNsub'',''rect'',''phaseFullSize'',''LcFull''']; %Tommy Sep 23, 2008
                if exist('creg')==1
                    savelist=[savelist,',''creg'''];
                    LNfull(oldrect(1):oldrect(3), oldrect(2):oldrect(4))=creg;
                    creg = LNfull(rect(1):rect(3), rect(2):rect(4));
                end
                if exist('yreg')==1
                    savelist=[savelist,',''yreg'''];
                    LNfull(oldrect(1):oldrect(3), oldrect(2):oldrect(4))=yreg;
                    yreg = LNfull(rect(1):rect(3), rect(2):rect(4));
                end
                if exist('greg')==1
                    savelist=[savelist,',''greg'''];
                    LNfull(oldrect(1):oldrect(3), oldrect(2):oldrect(4))=greg;
                    greg = LNfull(rect(1):rect(3), rect(2):rect(4));
                end
                if exist('rreg')==1
                    savelist=[savelist,',''rreg'''];
                    LNfull(oldrect(1):oldrect(3), oldrect(2):oldrect(4))=rreg;
                    rreg = LNfull(rect(1):rect(3), rect(2):rect(4));
                end
                
                crop_pop=0;

                g = double(phsub);
                g = (g-minmin(g)) / (maxmax(g)-minmin(g));
                g = g+0.2;
                gones = ones(size(g));
                g(g>1)= gones(g>1);
                g = uint8(g * 255);

                rescale=(p.min_size./size(LNsub));
                res=min(rescale);
                if res>2 res=2;end
                figure(phfig);
                clf reset;

                imshow(imresize(g,res));
                colormap((1:100)'*[0 1 0]/100)
                pos11=get(phfig,'position');
                set(phfig, 'position', ...
                    [p.leftend-pos11(3)-8, p.upend-pos11(4), ...
                     pos11(3), pos11(4)]);
                % set(phfig,'name',[name(outl+1:end-4)]);
                set(phfig,'name',['Frame ',str3(i),' phase']);
            end
            if savetemp
                LNsub=Lc;
                LcFull = zeros(phaseFullSize) ; %Tommy Sep 23, 2008
                LcFull(rect(1):rect(3),rect(2):rect(4)) = Lc ; %Tommy Sep 23, 2008
                savelist=[savelist,',''phsub'',''LNsub'',''rect'',''phaseFullSize'',''LcFull''']; %Tommy Sep 23, 2008
                if savetemp==2
                    tempsegcorrect=1;
                    eval(['save(''',name,''',',savelist,...
                          ',''tempsegcorrect'',''-append'');']);
                    disp(['Saved partial file ',name(outl:end),...
                          '   # of cells: ',num2str(max2(Lc))]);
                end
            end
        end
        if ~dontsave,
            tempsegcorrect=0;
            LcFull = zeros(phaseFullSize) ; %Tommy Sep 23, 2008
            LcFull(rect(1):rect(3),rect(2):rect(4)) = Lc ; %Tommy Sep 23, 2008
            savelist=[savelist,',''phsub'',''LNsub'',''rect'',''phaseFullSize'',''LcFull''']; %Tommy Sep 23, 2008
            eval(['save(''',name,''',',savelist,...
                  ',''tempsegcorrect'',''-append'');']);
            disp(['Updated file ',name(outl:end),'   # of cells: ',...
                  num2str(max2(Lc))]);
        elseif addtolist
            if exist(badseglistname)==2
                load(badseglistname);
            end
            if exist('badseglist')~=1
                badseglist=[];
            end
            badseglist=[badseglist,i];
            save(badseglistname,'badseglist');
            disp(['Added frame ',name(outl:end),' to bad segmentation list.']);
        else
            disp(['Skipped file ',name(outl:end),' <--']);
        end;
    end
    if backwards~=2
        loopindex = loopindex + p.Dskip;
    end
    if gotoframenum
        newloopindex = find(p.manualRange==gotoframenum);
        if isempty(newloopindex)
            disp(['Could not goto frame ', num2str(gotoframenum), ...
                ' because it is not in the manualRange.']);
        else
            loopindex = newloopindex;
        end
    end
end;
close(phfig);close(ourfig);
