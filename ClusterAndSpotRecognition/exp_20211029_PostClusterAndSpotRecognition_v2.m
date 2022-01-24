%%

clear all
clc

tic;

ip = exp_20211029_InitializeExp(); %Update the file name here
channels = {'Cy3', 'Cy5', 'GFP'};
channel_threshType = {'-nonatrous', '-atrous100', 'dualthresh'};
channel_threshNum = {'75', '0'};
channel_date = {'08-Nov-2021', '03-Dec-2021'};
sample_name = {'W3110 rpoc-GFP 1mM IPTG', 'W3110 rpoc-GFP 10uM IPTG', 'W3110 rpoc-GFP 0mM IPTG'};
gfp_datafolder = ['\\igolding-01.engr.illinois.edu\igolding\Shares\Patrick\' ... 
    '20211029_RNAPgfp_ribosomeldrFISH_lacZFISH\cluster_quantify\Run-GFP-0.4low0.53high-cellnormed-thresh-30-Nov-2021'];

rng(0);

mini_default = 1000;

%new structures
distancelist_new = {{zeros(0,8+3), zeros(0,8+3), zeros(0,8+3)},...
                    {zeros(0,3+3), zeros(0,3+3), zeros(0,3+3)}};

%{
    distancelist_new{1}{}(:,1) = LOW CLUSTER #
    distancelist_new{1}{}(:,2) = CLOSEST DIST TO LOW CLUSTER BOUNDARY
    distancelist_new{1}{}(:,3) = HIGH CLUSTER #
    distancelist_new{1}{}(:,4) = CLOSEST DIST TO high CLUSTER BOUNDARY
                    
    distancelist_new{1}{}(:,5) = CELL #
    distancelist_new{1}{}(:,6) = FRAME #
    distancelist_new{1}{}(:,7) = SPOT #
    distancelist_new{1}{}(:,8) = Cy3 spot intensity                    
    distancelist_new{1}{}(:,9) = Total GFP intensity within rad2 #1
    distancelist_new{1}{}(:,10) = Avg GFP int/pix within rad2 #1
    distancelist_new{1}{}(:,11) = Total GFP intensity within rad2 #2
    distancelist_new{1}{}(:,12) = Avg GFP int/pix within rad2 #2

    distancelist_new{1}{}(:,13) = CLOSEST CY5 SPOT #
    distancelist_new{1}{}(:,14) = DIST TO THAT CY5 SPOT
    distancelist_new{1}{}(:,15) = INT OF THAT CY5 SPOT
    distancelist_new{1}{}(:,16) = TOTAL # OF ALL CY5 SPOTS IN x RAD 1
    distancelist_new{1}{}(:,17) = TOTAL INT OF ALL CY5 SPOTS IN x RAD 1
    distancelist_new{1}{}(:,18) = TOTAL # OF ALL CY5 SPOTS IN x RAD 2
    distancelist_new{1}{}(:,19) = TOTAL INT OF ALL CY5 SPOTS IN x RAD 2
    distancelist_new{1}{}(:,20) = TOTAL # OF ALL CY5 SPOTS IN x RAD 3
    distancelist_new{1}{}(:,21) = TOTAL INT OF ALL CY5 SPOTS IN x RAD 3
                    
    -----------------------------------------------------------------------

    distancelist_new{2}{}(:,1) = LOW CLUSTER # of fake 
    distancelist_new{2}{}(:,2) = CLOSEST DIST TO LOW CLUSTER BOUNDARY of fake 
    distancelist_new{2}{}(:,3) = HIGH CLUSTER # of fake 
    distancelist_new{2}{}(:,4) = CLOSEST DIST TO high CLUSTER BOUNDARY of fake 
                    
    distancelist_new{2}{}(:,5) = CELL #
    distancelist_new{2}{}(:,6) = FRAME #
    distancelist_new{2}{}(:,7) = SPOT #
    distancelist_new{2}{}(:,8) = Cy3 spot intensity 
    distancelist_new{2}{}(:,9) = Total GFP intensity within rad2 #1 of fake 
    distancelist_new{2}{}(:,10) = Avg GFP int/pix within rad2 #1 of fake 
    distancelist_new{2}{}(:,11) = Total GFP intensity within rad2 #2 of fake 
    distancelist_new{2}{}(:,12) = Avg GFP int/pix within rad2 #2 of fake 
                    
    distancelist_new{2}{}(:,13) = CLOSEST CY5 SPOT # to fake
    distancelist_new{2}{}(:,14) = DIST TO THAT CY5 SPOT to fake
    distancelist_new{2}{}(:,15) = INT OF THAT CY5 SPOT
    distancelist_new{2}{}(:,16) = TOTAL # OF ALL CY5 SPOTS IN x RAD 1 of fake
    distancelist_new{2}{}(:,17) = TOTAL INT OF ALL CY5 SPOTS IN x RAD 1 of fake
    distancelist_new{2}{}(:,18) = TOTAL # OF ALL CY5 SPOTS IN x RAD 2 of fake
    distancelist_new{2}{}(:,19) = TOTAL INT OF ALL CY5 SPOTS IN x RAD 2 of fake
    distancelist_new{2}{}(:,20) = TOTAL # OF ALL CY5 SPOTS IN x RAD 3 of fake
    distancelist_new{2}{}(:,21) = TOTAL INT OF ALL CY5 SPOTS IN x RAD 3 of fake                 
                    
    ======================================================================= 

    distancelist_new{3}{}(:,1) = LOW CLUSTER #
    distancelist_new{3}{}(:,2) = CLOSEST DIST TO LOW CLUSTER BOUNDARY
    distancelist_new{3}{}(:,3) = HIGH CLUSTER #
    distancelist_new{3}{}(:,4) = CLOSEST DIST TO high CLUSTER BOUNDARY     
               
    distancelist_new{3}{}(:,5) = CELL #
    distancelist_new{3}{}(:,6) = FRAME #
    distancelist_new{3}{}(:,7) = SPOT #
    distancelist_new{3}{}(:,8) = Cy5 spot intensity
    distancelist_new{3}{}(:,9) = Total GFP intensity within rad2 #1
    distancelist_new{3}{}(:,10) = Avg GFP int/pix within rad2 #1
    distancelist_new{3}{}(:,11) = Total GFP intensity within rad2 #2
    distancelist_new{3}{}(:,12) = Avg GFP int/pix within rad2 #2
                    
    distancelist_new{3}{}(:,13) = CLOSEST CY5 SPOT #
    distancelist_new{3}{}(:,14) = DIST TO THAT CY5 SPOT
    distancelist_new{3}{}(:,15) = INT OF THAT CY5 SPOT
    distancelist_new{3}{}(:,16) = TOTAL # OF ALL CY5 SPOTS IN x RAD 1
    distancelist_new{3}{}(:,17) = TOTAL INT OF ALL CY5 SPOTS IN x RAD 1
    distancelist_new{3}{}(:,18) = TOTAL # OF ALL CY5 SPOTS IN x RAD 2
    distancelist_new{3}{}(:,19) = TOTAL INT OF ALL CY5 SPOTS IN x RAD 2
    distancelist_new{3}{}(:,20) = TOTAL # OF ALL CY5 SPOTS IN x RAD 3
    distancelist_new{3}{}(:,21) = TOTAL INT OF ALL CY5 SPOTS IN x RAD 3
                 
    -----------------------------------------------------------------------

    distancelist_new{4}{}(:,1) = LOW CLUSTER # of fake 
    distancelist_new{4}{}(:,2) = CLOSEST DIST TO LOW CLUSTER BOUNDARY of fake 
    distancelist_new{4}{}(:,3) = HIGH CLUSTER # of fake 
    distancelist_new{4}{}(:,4) = CLOSEST DIST TO high CLUSTER BOUNDARY of fake 
                    
    distancelist_new{4}{}(:,5) = CELL #
    distancelist_new{4}{}(:,6) = FRAME #
    distancelist_new{4}{}(:,7) = SPOT #
    distancelist_new{4}{}(:,8) = Cy5 spot intensity 
    distancelist_new{4}{}(:,9) = Total GFP intensity within rad2 #1 of fake 
    distancelist_new{4}{}(:,10) = Avg GFP int/pix within rad2 #1 of fake 
    distancelist_new{4}{}(:,11) = Total GFP intensity within rad2 #2 of fake 
    distancelist_new{4}{}(:,12) = Avg GFP int/pix within rad2 #2 of fake
                    
    distancelist_new{4}{}(:,13) = CLOSEST CY5 SPOT # to fake
    distancelist_new{4}{}(:,14) = DIST TO THAT CY5 SPOT to fake
    distancelist_new{4}{}(:,15) = INT OF THAT CY5 SPOT
    distancelist_new{4}{}(:,16) = TOTAL # OF ALL CY5 SPOTS IN x RAD 1 of fake
    distancelist_new{4}{}(:,17) = TOTAL INT OF ALL CY5 SPOTS IN x RAD 1 of fake
    distancelist_new{4}{}(:,18) = TOTAL # OF ALL CY5 SPOTS IN x RAD 2 of fake
    distancelist_new{4}{}(:,19) = TOTAL INT OF ALL CY5 SPOTS IN x RAD 2 of fake
    distancelist_new{4}{}(:,20) = TOTAL # OF ALL CY5 SPOTS IN x RAD 3 of fake
    distancelist_new{4}{}(:,21) = TOTAL INT OF ALL CY5 SPOTS IN x RAD 3 of fake      
%}
                    
rad(1) = 5;
rad(2) = 10;
rad(3) = 15;
% rad(4) = 20;

rad2(1) = 2;
rad2(2) = 4;

%Loop over samples
for sample_num = 1:length(sample_name)
    counter = 0;
    counter2 = 0;
    %Loop over frames
    for n_image = ip.sample{1,sample_num}.idx
        n_frame = ip.exp.splimg2frm(sample_num, n_image);
        
        %Progress update
        time = toc;
        hours = floor(time/3600);
        minutes = floor(time/60) - 60*hours;
        seconds = floor(time - 60*minutes - 3600*hours);
        
        disp([newline 'Loading Frame ', num2str(n_frame), ' of ', ...
            num2str(ip.exp.splimg2frm(end,end)) '. Elapsed time = ' ...
            num2str(hours,'%02d') 'hr:' ...
            num2str(minutes,'%02d') 'min:' ...
            num2str(seconds,'%02d') 'sec']);
        
        %Cell masks to iterate thru objects in the same cell only
        datafolder = ['\\igolding-01.engr.illinois.edu\igolding\Shares\Patrick\20211029_RNAPgfp_ribosomeldrFISH_lacZFISH\segmentation\masks'];
        cellmasks = load([datafolder '\sampleseg0' num2str(n_frame, '%02d') '.mat'], 'LcFull');
        
        %peakdata for spots in current frame and inside a cell
        datafolder = [pwd '\spots_quantify\Run-Cy3' channel_threshType{1} ...
            '-thresh' channel_threshNum{1} '-' channel_date{1} '-gated\data'];
        peakdata_cy3 = load([datafolder '\peakdata0' num2str(n_frame, '%02d')]);
        cy3_spotnum = 1:length(peakdata_cy3.peakdata);
        cy3_spotnum = cy3_spotnum((peakdata_cy3.peakdata(:,12) > 0) & (imag(peakdata_cy3.peakdata(:,14)) == 0));
        peakdata_cy3 = peakdata_cy3.peakdata(...
            (peakdata_cy3.peakdata(:,12) > 0) & (imag(peakdata_cy3.peakdata(:,14)) == 0),:);
        
        datafolder = [pwd '\spots_quantify\Run-Cy5' channel_threshType{2} ...
            '-thresh' channel_threshNum{2} '-' channel_date{2} '-gated\data'];
        peakdata_cy5 = load([datafolder '\peakdata0' num2str(n_frame, '%02d')]);
        cy5_spotnum = 1:length(peakdata_cy5.peakdata);
        cy5_spotnum = cy5_spotnum((peakdata_cy5.peakdata(:,12) > 0) & (imag(peakdata_cy5.peakdata(:,14)) == 0));
        peakdata_cy5 = peakdata_cy5.peakdata(...
            (peakdata_cy5.peakdata(:,12) > 0) & (imag(peakdata_cy5.peakdata(:,14)) == 0),:);
        
        %GFP fluor image
        all_low_cluster_mask = load([gfp_datafolder '\data\clusterdata_low0' ...
            num2str(n_frame, '%02d') '.mat'], 'Lcluster_low');
        all_high_cluster_mask = load([gfp_datafolder '\data\clusterdata_high0' ...
            num2str(n_frame, '%02d') '.mat'], 'Lcluster_high');        
        
        spot_folder = [ip.exp.path '\cluster_quantify\Run-GFP-0.4low0.53high-cellnormed-thresh-07-Dec-2021\'];
        sr = InitializeClusterRecognitionParameters_c234(ip,n_frame,4,spot_folder);
        z_range = sr.image.zrange;
        N = max(double(cellmasks.LcFull(:)));
        all_cv = zeros(numel(z_range),N);
        fprintf(1,['Finding best z for each cell.' sprintf('\n')]);
        for i = 1:numel(z_range)
            slice_name = [sr.image.fullname 'z' num2str(z_range(i), '%01d') 'c' num2str(sr.image.channel, '%01d') '.tif'];
            FluorImage = imread(slice_name);
            
            for ii = 1:N
                cell = immultiply(cellmasks.LcFull == ii, FluorImage);
                all_cv(i, ii) = sqrt(var(double(cell(cellmasks.LcFull(:) == ii)))) / mean(double(cell(cellmasks.LcFull(:) == ii)));
            end
        end
        [~, max_cv_z] = max(all_cv);
        max_FluorImage = zeros(size(cellmasks.LcFull));
        fprintf(1,['Applying best z for each cell.' sprintf('\n')]);
        for i = 1:numel(z_range)
            slice_name = [sr.image.fullname 'z' num2str(z_range(i), '%01d') 'c' num2str(sr.image.channel, '%01d') '.tif'];
            FluorImage = imread(slice_name);          
            for ii = 1:N
                if max_cv_z(ii) == i
                    max_FluorImage = max_FluorImage + double(immultiply(cellmasks.LcFull == ii, FluorImage));
                end
            end
        end
        
        %Look at Cy3 spots
        progress_ = '';
        for spot_num = 1:height(peakdata_cy3)
            counter = counter + 1;
            
            for d_=1:1:size(progress_,2) ; fprintf(1,'\b') ; end
            progress_ = ['    Looking at Cy3 spot #' num2str(spot_num) ' of ' num2str(height(peakdata_cy3)) '.'] ;
            fprintf(1,progress_);
            
            for fake_spot = [0 1]
                cellnum = peakdata_cy3(spot_num, 12);
                
                if ~fake_spot
                    x = peakdata_cy3(spot_num, 2);
                    y = peakdata_cy3(spot_num, 3);
                else
                    cell_coor = cellmasks.LcFull == cellnum;
                    
                    which_coor = randi([1 sum(cell_coor(:))], 1,1);
                    coor = find(cell_coor', which_coor);
                    coor = coor(end);
                    
                    x = mod(coor, ip.image.size);
                    y = ceil(coor / ip.image.size) + 1;
                end

                %is inside a low cluster?
                low_cluster_mask = all_low_cluster_mask.Lcluster_low .* (cellmasks.LcFull == cellnum);

                inside = low_cluster_mask(round(y),round(x));
                if isempty(inside)
                    inside = 0;
                end
                %{ 
                    visualize the overlap:
                    % red
                    img(:,:,1) = uint8(255 * logical(low_cluster_mask));
                    banana = zeros(size(low_cluster_mask));
                    banana(round(y), round(x)) = 1;
                    banana(round(y)+1, round(x)) = 1;
                    banana(round(y)-1, round(x)) = 1;
                    banana(round(y), round(x)+1) = 1;
                    banana(round(y), round(x)-1) = 1;
                    % green
                    perim = bwperim(cellmasks.LcFull == cellnum);
                    img(:,:,2) = uint8(255 * perim);
                    % blue
                    img(:,:,3) = uint8(255 * logical(banana));
                    imshow(img);
                %}

                %YES: distance to edge of cluster, cluster# & dist
                if inside
                    perim = bwperim(low_cluster_mask == inside);
                    perim_ind = find(perim ~= 0);
                    perim_ind = [ floor(perim_ind / 1608), mod(perim_ind,1608) ];

                    dist = sqrt((x - perim_ind(:,1)).^2 + (y - perim_ind(:,2)).^2);
                    [mini, ~] = min(dist);
                %NO: distance to edge of closest cluster, cluster# * -1 & dist
                else
                    perim = bwperim(low_cluster_mask);
                    perim_ind = find(perim ~= 0);
                    perim_ind = [ floor(perim_ind / 1608), mod(perim_ind,1608) ];

                    dist = sqrt((x - perim_ind(:,1)).^2 + (y - perim_ind(:,2)).^2);
                    [mini, ind] = min(dist);
                    inside = -1 * low_cluster_mask( ...
                            perim_ind(ind,1)*1608 + perim_ind(ind,2) );
                end
                if isempty(inside)
                    inside = 0;
                    mini = mini_default;
                end
                distancelist_new{1 + fake_spot}{sample_num}(counter,1) = inside;
                distancelist_new{1 + fake_spot}{sample_num}(counter,2) = mini;
                %end if in low cluster

                %is inside high cluster?, 0 or 1
                high_cluster_mask = all_high_cluster_mask.Lcluster_high .* (cellmasks.LcFull == cellnum);

                inside = high_cluster_mask(round(y),round(x));
                if isempty(inside)
                    inside = 0;
                end
                %YES: distance to edge of cluster, cluster# & dist
                if inside
                    perim = bwperim(high_cluster_mask == inside);
                    perim_ind = find(perim ~= 0);
                    perim_ind = [ floor(perim_ind / 1608), mod(perim_ind,1608) ];

                    dist = sqrt((x - perim_ind(:,1)).^2 + (y - perim_ind(:,2)).^2);
                    [mini, ~] = min(dist);
                %NO: distance to edge of closest cluster, cluster# * -1 & dist
                else
                    perim = bwperim(high_cluster_mask);
                    perim_ind = find(perim ~= 0);
                    perim_ind = [ floor(perim_ind / 1608), mod(perim_ind,1608) ];

                    dist = sqrt((x - perim_ind(:,1)).^2 + (y - perim_ind(:,2)).^2);
                    [mini, ind] = min(dist);
                    inside = -1 * high_cluster_mask( ...
                            perim_ind(ind,1)*1608 + perim_ind(ind,2) );
                end
                if isempty(inside)
                    inside = 0;
                    mini = mini_default;
                end
                distancelist_new{1 + fake_spot}{sample_num}(counter,3) = inside;
                distancelist_new{1 + fake_spot}{sample_num}(counter,4) = mini;
                %end if in high cluster

                %closest rRNA spot, spot# & dist & int
                %   also
                %count rRNA within X rad, # of spots & total int of those spots
                rolling_int = [0 0 0];
                rolling_spotnum = [0 0 0];
                temp_cy5_data = peakdata_cy5(peakdata_cy5(:,12) == cellnum,:);
                spotcount = 1:length(peakdata_cy5);
                spotcount = spotcount(peakdata_cy5(:,12) == cellnum);

                dist = sqrt((x - temp_cy5_data(:, 2)).^2 + (y - temp_cy5_data(:, 3)).^2);
                [mini, ind] = min(dist);
                if isempty(mini)
                    mini = mini_default;
                    int = 0;
                    spot = 0;
                else
                    int = temp_cy5_data(ind, 14);
                    spot = spotcount(ind);
                end

                ind = dist < rad(1);
                rolling_int(1) = sum(temp_cy5_data(ind, 14));
                rolling_spotnum(1) = sum(ind);

                ind = dist < rad(2);
                rolling_int(2) = sum(temp_cy5_data(ind, 14));
                rolling_spotnum(2) = sum(ind);

                ind = dist < rad(3);
                rolling_int(3) = sum(temp_cy5_data(ind, 14));
                rolling_spotnum(3) = sum(ind);

                %%%
                distancelist_new{1 + fake_spot}{sample_num}(counter,5) = cellnum;
                distancelist_new{1 + fake_spot}{sample_num}(counter,6) = n_frame;
                distancelist_new{1 + fake_spot}{sample_num}(counter,7) = cy3_spotnum(spot_num);
                distancelist_new{1 + fake_spot}{sample_num}(counter,8) = peakdata_cy3(spot_num,14);
                
                %calculate GFP fluor inside the spot
                for iii = 1:length(rad2)
                    spot_mask = DrawEllipse(x,y,rad2(iii),rad2(iii),0);
                    coords = unique([uint16(round(spot_mask.YData)), uint16(round(spot_mask.XData))], 'rows');
                    
                    spot_mask = zeros(ip.image.size, ip.image.size);
                    spot_mask(sub2ind([ip.image.size, ip.image.size], coords(:,1), coords(:,2))) = 1;
                    spot_mask = uint16(imfill(spot_mask)) .* uint16(cellmasks.LcFull == cellnum);
                    
                    gfp_fluor = max_FluorImage(logical(spot_mask));
                    
                    distancelist_new{1 + fake_spot}{sample_num}(counter, 9 + (iii-1)*2) = ...
                        sum(gfp_fluor);
                    distancelist_new{1 + fake_spot}{sample_num}(counter, 10 + (iii-1)*2) = ...
                        mean(gfp_fluor);
                end
                %
                distancelist_new{1 + fake_spot}{sample_num}(counter,13) = spot;
                distancelist_new{1 + fake_spot}{sample_num}(counter,14) = mini;
                distancelist_new{1 + fake_spot}{sample_num}(counter,15) = int;

                distancelist_new{1 + fake_spot}{sample_num}(counter,16) = rolling_spotnum(1);
                distancelist_new{1 + fake_spot}{sample_num}(counter,17) = rolling_int(1);
                distancelist_new{1 + fake_spot}{sample_num}(counter,18) = rolling_spotnum(2);
                distancelist_new{1 + fake_spot}{sample_num}(counter,19) = rolling_int(2);            
                distancelist_new{1 + fake_spot}{sample_num}(counter,20) = rolling_spotnum(3);
                distancelist_new{1 + fake_spot}{sample_num}(counter,21) = rolling_int(3);
            end
        end
        %end looking at Cy3
        
        %%%
        
        %Look at Cy5 spots
        disp(' Done!');
        progress_ = '';
        for spot_num = 1:height(peakdata_cy5)
            counter2 = counter2 + 1;
            
            for d_=1:1:size(progress_,2) ; fprintf(1,'\b') ; end
            progress_ = ['    Looking at Cy5 spot #' num2str(spot_num) ' of ' num2str(height(peakdata_cy5)) '.'] ;
            fprintf(1,progress_);
            
            for fake_spot = [0 1]
                if ~fake_spot
                    x = peakdata_cy5(spot_num, 2);
                    y = peakdata_cy5(spot_num, 3);
                else
                    cell_coor = cellmasks.LcFull == cellnum;
                    
                    which_coor = randi([1 sum(cell_coor(:))], 1,1);
                    coor = find(cell_coor', which_coor);
                    coor = coor(end);
                    
                    x = mod(coor, ip.image.size);
                    y = ceil(coor / ip.image.size) + 1;
                end

                cellnum = peakdata_cy5(spot_num, 12);

                %is inside a low cluster?
                low_cluster_mask = all_low_cluster_mask.Lcluster_low .* (cellmasks.LcFull == cellnum);

                inside = low_cluster_mask(round(y),round(x));
                if isempty(inside)
                    inside = 0;
                end

                %YES: distance to edge of cluster, cluster# & dist
                if inside
                    perim = bwperim(low_cluster_mask == inside);
                    perim_ind = find(perim ~= 0);
                    perim_ind = [ floor(perim_ind / 1608), mod(perim_ind,1608) ];

                    dist = sqrt((x - perim_ind(:,1)).^2 + (y - perim_ind(:,2)).^2);
                    [mini, ~] = min(dist);
                %NO: distance to edge of closest cluster, cluster# * -1 & dist
                else
                    perim = bwperim(low_cluster_mask);
                    perim_ind = find(perim ~= 0);
                    perim_ind = [ floor(perim_ind / 1608), mod(perim_ind,1608) ];

                    dist = sqrt((x - perim_ind(:,1)).^2 + (y - perim_ind(:,2)).^2);
                    [mini, ind] = min(dist);
                    inside = -1 * low_cluster_mask( ...
                            perim_ind(ind,1)*1608 + perim_ind(ind,2) );
                end
                if isempty(inside)
                    inside = 0;
                    mini = mini_default;
                end
                distancelist_new{3 + fake_spot}{sample_num}(counter2,1) = inside;
                distancelist_new{3 + fake_spot}{sample_num}(counter2,2) = mini;
                %end if in low cluster

                %is inside high cluster?, 0 or 1
                high_cluster_mask = all_high_cluster_mask.Lcluster_high .* (cellmasks.LcFull == cellnum);

                inside = high_cluster_mask(round(y),round(x));
                if isempty(inside)
                    inside = 0;
                end

                %YES: distance to edge of cluster, cluster# & dist
                if inside
                    perim = bwperim(high_cluster_mask == inside);
                    perim_ind = find(perim ~= 0);
                    perim_ind = [ floor(perim_ind / 1608), mod(perim_ind,1608) ];

                    dist = sqrt((x - perim_ind(:,1)).^2 + (y - perim_ind(:,2)).^2);
                    [mini, ~] = min(dist);
                %NO: distance to edge of closest cluster, cluster# * -1 & dist
                else
                    perim = bwperim(high_cluster_mask);
                    perim_ind = find(perim ~= 0);
                    perim_ind = [ floor(perim_ind / 1608), mod(perim_ind,1608) ];

                    dist = sqrt((x - perim_ind(:,1)).^2 + (y - perim_ind(:,2)).^2);
                    [mini, ind] = min(dist);
                    inside = -1 * high_cluster_mask( ...
                            perim_ind(ind,1)*1608 + perim_ind(ind,2) );
                end
                if isempty(inside)
                    inside = 0;
                    mini = mini_default;
                end
                distancelist_new{3 + fake_spot}{sample_num}(counter2,3) = inside;
                distancelist_new{3 + fake_spot}{sample_num}(counter2,4) = mini;
                %end if in high cluster

                %closest rRNA spot, spot# & dist & int
                %   also
                %count rRNA within X rad, # of spots & total int of those spots
                rolling_int = [0 0 0];
                rolling_spotnum = [0 0 0];
                temp_cy5_data = peakdata_cy5(peakdata_cy5(:,12) == cellnum,:);
                spotcount = 1:length(peakdata_cy5);
                spotcount = spotcount(peakdata_cy5(:,12) == cellnum);

                dist = sqrt((x - temp_cy5_data(:, 2)).^2 + (y - temp_cy5_data(:, 3)).^2);
                %Don't include yourself & adjust other dataset to match ind
                temp_cy5_data = temp_cy5_data(dist ~= 0,:);
                dist = dist(dist ~= 0);

                [mini, ind] = min(dist);
                if isempty(mini)
                    mini = mini_default;
                    int = 0;
                    spot = 0;
                else
                    int = temp_cy5_data(ind, 14);
                    spot = spotcount(ind);
                end

                ind = dist < rad(1);
                rolling_int(1) = sum(temp_cy5_data(ind, 14));
                rolling_spotnum(1) = sum(ind);

                ind = dist < rad(2);
                rolling_int(2) = sum(temp_cy5_data(ind, 14));
                rolling_spotnum(2) = sum(ind);

                ind = dist < rad(3);
                rolling_int(3) = sum(temp_cy5_data(ind, 14));
                rolling_spotnum(3) = sum(ind);

                %%%
                distancelist_new{3 + fake_spot}{sample_num}(counter2,5) = cellnum;
                distancelist_new{3 + fake_spot}{sample_num}(counter2,6) = n_frame;
                distancelist_new{3 + fake_spot}{sample_num}(counter2,7) = cy5_spotnum(spot_num);
                distancelist_new{3 + fake_spot}{sample_num}(counter2,8) = peakdata_cy5(spot_num,14);            
                %calculate GFP fluor inside the spot
                for iii = 1:length(rad2)
                    spot_mask = DrawEllipse(x,y,rad2(iii),rad2(iii),0);
                    coords = unique([uint16(round(spot_mask.YData)), uint16(round(spot_mask.XData))], 'rows');
                    
                    spot_mask = zeros(ip.image.size, ip.image.size);
                    spot_mask(sub2ind([ip.image.size, ip.image.size], coords(:,1), coords(:,2))) = 1;
                    spot_mask = uint16(imfill(spot_mask)) .* uint16(cellmasks.LcFull == cellnum);
                    
                    gfp_fluor = max_FluorImage(logical(spot_mask));
                    
                    distancelist_new{3 + fake_spot}{sample_num}(counter2, 9 + (iii-1)*2) = ...
                        sum(gfp_fluor);
                    distancelist_new{3 + fake_spot}{sample_num}(counter2, 10 + (iii-1)*2) = ...
                        mean(gfp_fluor);
                end
                %
                distancelist_new{3 + fake_spot}{sample_num}(counter2,13) = spot;
                distancelist_new{3 + fake_spot}{sample_num}(counter2,14) = mini;
                distancelist_new{3 + fake_spot}{sample_num}(counter2,15) = int;

                distancelist_new{3 + fake_spot}{sample_num}(counter2,16) = rolling_spotnum(1);
                distancelist_new{3 + fake_spot}{sample_num}(counter2,17) = rolling_int(1);
                distancelist_new{3 + fake_spot}{sample_num}(counter2,18) = rolling_spotnum(2);
                distancelist_new{3 + fake_spot}{sample_num}(counter2,19) = rolling_int(2);            
                distancelist_new{3 + fake_spot}{sample_num}(counter2,20) = rolling_spotnum(3);
                distancelist_new{3 + fake_spot}{sample_num}(counter2,21) = rolling_int(3);
            end
        end
        disp(' Done!');
        %end looking at Cy5
    end
    %end looking at frames
end
%end looking at samples

if ~exist([pwd '\spots&cluster_quantify'],'dir')
    mkdir([pwd '\spots&cluster_quantify'])
end

save([pwd '\spots&cluster_quantify\distancelist_new.mat'], 'distancelist_new');


%% Function
function ell = DrawEllipse(XPos, YPos, minor_ax, major_ax, angle)
    r = 0:0.01:2*pi;                                        %Full circle
    p = [(minor_ax * cos(r))' (major_ax * sin(r))'];        %Ellipse before rotate and move
    alpha = [ cos(angle) sin(angle);
             -sin(angle) cos(angle)];
    p1 = p * alpha;                                         %Ellipse after rotate
    p2 = p1 + [XPos YPos];                                  %Ellipse after move
%     figure(); ell = patch('xdata', p2(:, 1), 'ydata', p2(:, 2), 'facecolor', [0 0 1], 'edgecolor', 'none');
    ell.XData = p2(:,1);
    ell.YData = p2(:,2);
end

