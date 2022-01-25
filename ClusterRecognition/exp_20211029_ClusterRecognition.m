<<<<<<< HEAD
clc
clear

ip = exp_20211029_InitializeExp ; 
% samples channels
%   1   2
%   ph3	TXRED 

spl2cnum = struct() ; 
spl2cnum.Cy5 = 2*ones(1,max(ip.exp.sampleList)) ;
spl2cnum.Cy3 = 3*ones(1,max(ip.exp.sampleList)) ;
spl2cnum.GFP = 4*ones(1,max(ip.exp.sampleList)) ;
spl2cnum.DAPI = 5*ones(1,max(ip.exp.sampleList)) ;


%% test

clc
close all;
clearvars -except spl2cnum ip;

% spot_folder = [ip.exp.path '\cluster_quantify\Run-GFP-.4low.6high-cellnormed-thresh-' date '\'];
norm_over_all = 0;

tic ;
for combo = 3
    if combo == 1
        high_threshold = .43;%prctile(norm_nobg_fluor,90);
        low_threshold = .3;%prctile(norm_nobg_fluor,50);
    elseif combo == 2
        high_threshold = .52;%prctile(norm_nobg_fluor,90);
        low_threshold = .39;%prctile(norm_nobg_fluor,50);
    elseif combo == 3
        high_threshold = .53;%prctile(norm_nobg_fluor,90);
        low_threshold = .4;%prctile(norm_nobg_fluor,50);
    elseif combo == 4
        high_threshold = .54;%prctile(norm_nobg_fluor,90);
        low_threshold = .41;%prctile(norm_nobg_fluor,50);
    elseif combo == 5
        high_threshold = .63;%prctile(norm_nobg_fluor,90);
        low_threshold = .5;%prctile(norm_nobg_fluor,50);
    end
    
    spot_folder = [ip.exp.path '\cluster_quantify\Run-GFP-' num2str(low_threshold) ...
        'low' num2str(high_threshold) 'high-cellnormed-thresh-' date '\'];
    
for n_spl = [1 2 3]
%     if n_spl == 1
%         imgs = [1 12];
%     elseif n_spl == 2
%         imgs = [7];
%     else
%         imgs = [8];
%     end
    
    for n_img =   ip.sample{n_spl}.idx
            
        n_frame = ip.exp.splimg2frm(n_spl , n_img) ;
        time = toc;
        hours = floor(time/3600);
        minutes = floor(time/60) - 60*hours;
        seconds = floor(time - 60*minutes - 3600*hours);

        progress_ = ['SPOT RECOGNITION. FRAME ' num2str(n_frame,'%4d')...
            ' of ' num2str(ip.exp.totalframes,'%4d') '. Elapsed time = ' ...
            num2str(hours,'%02d') 'hr:' ...
            num2str(minutes,'%02d') 'min:' ...
            num2str(seconds,'%02d') 'sec'] ;
        fprintf(1,[progress_ sprintf('\n')]) ;
        n_channel = spl2cnum.GFP(n_spl) ;
        sr = InitializeClusterRecognitionParameters_c234(ip,n_frame,n_channel,spot_folder);
        
        %%%
        
        load([sr.seg.dir sr.seg.name]);
        z_range = sr.image.zrange;
        LcFull = uint16(LcFull);
        ALL_cellPerim = zeros(size(LcFull)) ;
        for iCell = 1:1:max(LcFull(:))
            ALL_cellPerim = ALL_cellPerim + bwperim(LcFull==iCell) ;
        end
        N = max(double(LcFull(:)));
        
        all_cv = zeros(numel(z_range),N);
        fprintf(1,['Finding best z for each cell.' sprintf('\n')]);
        for i = 1:numel(z_range)
            %fprintf(1,['z-Slice ' num2str(i) ' out of ' num2str(numel(z_range)) '.' sprintf('\n')]);
            
            slice_name = [sr.image.fullname 'z' num2str(z_range(i), '%01d') 'c' num2str(sr.image.channel, '%01d') '.tif'];
            FluorImage = imread(slice_name);
            
            for ii = 1:N
                cell = immultiply(LcFull == ii, FluorImage);
                all_cv(i, ii) = sqrt(var(double(cell(LcFull(:) == ii)))) / mean(double(cell(LcFull(:) == ii)));
            end
        end
        
        [~, max_cv_z] = max(all_cv);
        max_FluorImage = zeros(size(LcFull));
        max_normed = zeros(size(LcFull));
        fprintf(1,['Applying best z for each cell.' sprintf('\n')]);
        for i = 1:numel(z_range)
            slice_name = [sr.image.fullname 'z' num2str(z_range(i), '%01d') 'c' num2str(sr.image.channel, '%01d') '.tif'];
            FluorImage = imread(slice_name);          
            
            for ii = 1:N
                if max_cv_z(ii) == i
                max_FluorImage = max_FluorImage + double(immultiply(LcFull == ii, FluorImage));
                
                normed_cell = double(immultiply(LcFull == ii, FluorImage));
                normed_cell = normed_cell - double(min(FluorImage(LcFull == ii)));
                normed_cell = normed_cell / max(normed_cell(:));
                max_normed = max_normed + double(immultiply(LcFull == ii, normed_cell));
                end
            end
        end
        %only cells
        cell_fluor = max_FluorImage(logical(LcFull));
        if norm_over_all
        %normalize only cells and on pixels within the percentiles
            pix_to_analyze = cell_fluor((cell_fluor < prctile(cell_fluor(:),98)) ...
                & (cell_fluor > prctile(cell_fluor(:),2)));
            min_to_sub = min(double(pix_to_analyze(:)));
            pix_to_analyze = double(pix_to_analyze) - min_to_sub;
            max_to_divide = max(pix_to_analyze(:));

            %apply normalization to whole frame
            norm_fluor = double(max_FluorImage) - min_to_sub;
            norm_fluor = norm_fluor/max_to_divide;
            
            %set threshold within normalized range
            high_threshold = .8;%prctile(norm_nobg_fluor,90);
            low_threshold = .4;%prctile(norm_nobg_fluor,50);
        else
            norm_fluor = max_normed;
            %set threshold within normalized range
%             high_threshold = .53;%prctile(norm_nobg_fluor,90);
%             low_threshold = .4;%prctile(norm_nobg_fluor,50);
        end
        
        %apply threshold
        T_high_image = logical(norm_fluor >= high_threshold);
        T_high_image = T_high_image & LcFull;
        T_low_image = logical(norm_fluor >= low_threshold);
        T_low_image = T_low_image & LcFull;

        smooth_T_high_image = imfill(T_high_image, 'holes');
        se = strel('disk',3);            
        smooth_T_high_image = imerode(smooth_T_high_image,se);
        smooth_T_high_image = imdilate(smooth_T_high_image,se);

        smooth_T_low_image = imfill(T_low_image, 'holes');
        se = strel('disk',2);            
        smooth_T_low_image = imerode(smooth_T_low_image,se);
        smooth_T_low_image = imdilate(smooth_T_low_image,se);

        f1 = figure('Units', 'Normalized', 'OuterPosition', [0.05 0.05 0.95/2 0.90]);
%         histogram(norm_fluor(norm_fluor > 0 & norm_fluor < 1));
        subplot(2,1,1);
        xlim([0 1]);
        nbin = linspace(0, 1, 100) ;
        [ya,xa] = hist(norm_fluor(norm_fluor > 0 & norm_fluor < 1),nbin) ;
%             plotDoubleGaussFit_gen(xa(1:end-1) , ya(1:end-1) , nbin) ;
%             plotThreeGaussFit_gen(xa(1:end-1) , ya(1:end-1) , nbin) ;
        plotMultiGaussFit_real(xa(1:end-1) , ya(1:end-1) , nbin) ;
        ylabel('Probability') ;
        xlabel('Fluor. intensity within cells') ;
        legend('show')
        title([ip.sample{n_spl}.name ': Norm. int. in a cell']);
        subplot(2,1,2);
        yyaxis left;
        ylabel('Probability') ; hold on;
        h = bar(xa,ya/sum(ya),'EdgeColor','none','FaceColor', 0.9 * [1 1 1] ,...
            'BarWidth',.9, 'displayname', 'Data histogram'); hold all ; 
        yyaxis right
        ylabel('Fraction of cell');
        colors = {'#F00','#F80','#FF0','#0B0','#00F'};
        nbin = linspace(0, 1 ,51) ;
        for i = 1:5
            cell_norm_fluor = norm_fluor(LcFull == i);
            [N, edges] = histcounts(cell_norm_fluor, 'BinEdges', nbin, 'Normalization', 'probability');
            edges = (edges(1:end-1) + edges(2:end) ) / 2;
            plot(edges, N, 'DisplayName', ['Sample single cell histogram #', num2str(i)], ...
                'LineWidth', 2, 'LineStyle', '-', 'Color', colors{i}, 'Marker', '.');
            hold on;
        end
        legend('Location', 'best');
        xlabel('Norm fluor. intensity within cells') ;
        %set(gca, 'YScale', 'log')
        %set(gca, 'XScale', 'log')

        f2 = figure('numbertitle', 'off', 'name' , ['Reconstruction: frame No ' num2str(n_frame)], ...
            'Units', 'Normalized', 'OuterPosition', [0.05 0.05 0.95 0.75]);
        ha = tight_subplot(1, 3, [.01 .01], [.1 .01], [.01 .01]);
        axes(ha(1)); hold on;
        imshow(norm_fluor, [0 1]); %imcontrast
        % cell outline
        spy(ALL_cellPerim,'w',0.5); hold on ; xlabel('');
        title('Raw colored image');
        axis square; daspect([1 1 1]);

        axes(ha(2)); hold on;
        rgb_image = double(cat(3,T_low_image & ~T_high_image,T_low_image,T_low_image));
        imshow(rgb_image);
        % cell outline
        spy(ALL_cellPerim,'w',0.5); hold on ; xlabel('');
        title('Reconstructed image from clusters');
        axis square; daspect([1 1 1]);

        axes(ha(3)); hold on;
        rgb_image = double(cat(3,smooth_T_low_image & ~smooth_T_high_image,smooth_T_low_image,smooth_T_low_image));
        imshow(rgb_image);
        % cell outline
        spy(ALL_cellPerim,'w',0.5); hold on ; xlabel('');
        title('Cleaned Reconstructed image from clusters');
        axis square; daspect([1 1 1]);

        linkaxes(ha, 'xy');

        savefig(f1, [sr.spotrec.output 'clusterdata_histogram' num2str(n_frame,'%03d') '.fig']);
        savefig(f2, [sr.spotrec.output 'clusterdata_reconstruction' num2str(n_frame,'%03d') '.fig']);
        saveas(f1, [sr.spotrec.output 'clusterdata_histogram' num2str(n_frame,'%03d') '.png']);
        saveas(f2, [sr.spotrec.output 'clusterdata_reconstruction' num2str(n_frame,'%03d') '.png']);

        close all
        
        %%%
        
        Lcluster_low = bwlabel(smooth_T_low_image);
        Lcluster_low = Lcluster_low .* (Lcluster_low & ~smooth_T_high_image);
        if exist([sr.output 'clusterdata_low' num2str(n_frame,'%03d') '.mat'],'file')
            save([sr.output 'clusterdata_low' num2str(n_frame,'%03d') '.mat'],'Lcluster_low', '-append');
        else
            save([sr.output 'clusterdata_low' num2str(n_frame,'%03d') '.mat'],'Lcluster_low');
        end
        
        Lcluster_high = bwlabel(smooth_T_high_image);
        if exist([sr.output 'clusterdata_high' num2str(n_frame,'%03d') '.mat'],'file')
            save([sr.output 'clusterdata_high' num2str(n_frame,'%03d') '.mat'],'Lcluster_high', '-append');
        else
            save([sr.output 'clusterdata_high' num2str(n_frame,'%03d') '.mat'],'Lcluster_high');
        end
        
        %%%
        
        clusterdata_low = zeros(max(Lcluster_low(:)), 10);
        clusterdata_high = zeros(max(Lcluster_high(:)), 10);
        for iii = 1:max(Lcluster_low(:))
            %cell number
            clusterdata_low(iii, 1) = median(LcFull(Lcluster_low == iii));
            %z slice of cell
            clusterdata_low(iii, 2) = max_cv_z(clusterdata_low(iii, 1));
            %intensity total
            clusterdata_low(iii, 3) = sum(max_FluorImage(Lcluster_low == iii));
            %area
            clusterdata_low(iii, 4) = sum(sum(Lcluster_low == iii));
            %intensity per pix
            clusterdata_low(iii, 5) = clusterdata_low(iii, 3) / clusterdata_low(iii, 4);
            %frame number
            clusterdata_low(iii, 6) = n_frame;
            %X pixel position
            x_coords = repmat(1:ip.image.size, [ip.image.size, 1]);
            clusterdata_low(iii, 7) = mean( x_coords(Lcluster_low == iii) );
            %Y pixel position
            y_coords = repmat([1:ip.image.size]', [1, ip.image.size]);
            clusterdata_low(iii, 8) = mean( y_coords(Lcluster_low == iii) );
            %cluster number
            clusterdata_low(iii, 9) = iii;
        end
        
        for iii = 1:max(Lcluster_high(:))
            %cell number
            clusterdata_high(iii, 1) = median(LcFull(Lcluster_high == iii));
            %z slice of cell
            clusterdata_high(iii, 2) = max_cv_z(clusterdata_high(iii, 1));
            %intensity total
            clusterdata_high(iii, 3) = sum(max_FluorImage(Lcluster_high == iii));
            %area
            clusterdata_high(iii, 4) = sum(sum(Lcluster_high == iii));
            %intensity per pix
            clusterdata_high(iii, 5) = clusterdata_high(iii, 3) / clusterdata_high(iii, 4);
            %frame number
            clusterdata_high(iii, 6) = n_frame;
            %X pixel position
            x_coords = repmat(1:ip.image.size, [ip.image.size, 1]);
            clusterdata_high(iii, 7) = mean( x_coords(Lcluster_high == iii) );
            %Y pixel position
            y_coords = repmat([1:ip.image.size]', [1, ip.image.size]);
            clusterdata_high(iii, 8) = mean( y_coords(Lcluster_high == iii) );
            %low cluster number
            clusterdata_high(iii, 9) = median(Lcluster_low(Lcluster_high == iii));
            %high cluster number
            clusterdata_high(iii, 10) = iii;
        end
        
        if exist([sr.output 'clusterdata_low' num2str(n_frame,'%03d') '.mat'],'file')
            save([sr.output 'clusterdata_low' num2str(n_frame,'%03d') '.mat'],'clusterdata_low','-append');
        else
            save([sr.output 'clusterdata_low' num2str(n_frame,'%03d') '.mat'],'clusterdata_low');
        end
        
        if exist([sr.output 'clusterdata_high' num2str(n_frame,'%03d') '.mat'],'file')
            save([sr.output 'clusterdata_high' num2str(n_frame,'%03d') '.mat'],'clusterdata_high','-append');
        else
            save([sr.output 'clusterdata_high' num2str(n_frame,'%03d') '.mat'],'clusterdata_high');
        end
    end
end
time = toc;
hours = floor(time/3600);
minutes = floor(time/60) - 60*hours;
seconds = floor(time - 60*minutes - 3600*hours);

progress_ = ['Done. Elapsed time = ' ...
    num2str(hours,'%02d') 'hr:' ...
    num2str(minutes,'%02d') 'min:' ...
    num2str(seconds,'%02d') 'sec'] ;
fprintf(1,[progress_ sprintf('\n')]);      
end


%% run

clc
close all;
clearvars -except spl2cnum ip;

spot_folder = [ip.exp.path '\cluster_quantify\Run-TXRED-.4low.7high-cellnormed-thresh-' date '\'];
norm_over_all = 0;

tic ;

for n_spl = [10]%ip.exp.sampleList
    for n_img =   ip.sample{n_spl}.idx
        
        close all;
        
        n_frame = ip.exp.splimg2frm(n_spl , n_img) ;
        time = toc;
        hours = floor(time/3600);
        minutes = floor(time/60) - 60*hours;
        seconds = floor(time - 60*minutes - 3600*hours);

        progress_ = ['SPOT RECOGNITION. FRAME ' num2str(n_frame,'%4d')...
            ' of ' num2str(ip.exp.totalframes,'%4d') '. Elapsed time = ' ...
            num2str(hours,'%02d') 'hr:' ...
            num2str(minutes,'%02d') 'min:' ...
            num2str(seconds,'%02d') 'sec'] ;
        fprintf(1,[progress_ sprintf('\n')]) ;
        n_channel = spl2cnum.TXRED(n_spl) ;
        sr = InitializeClusterRecognitionParameters_c234(ip,n_frame,n_channel,spot_folder);
        
        %%%
        
        load([sr.seg.dir sr.seg.name]);
        z_range = sr.image.zrange;
        LcFull = uint16(LcFull);
        ALL_cellPerim = zeros(size(LcFull)) ;
        for iCell = 1:1:max(LcFull(:))
            ALL_cellPerim = ALL_cellPerim + bwperim(LcFull==iCell) ;
        end
        N = max(double(LcFull(:)));
        
        all_cv = zeros(numel(z_range),N);
        fprintf(1,['Finding best z for each cell.' sprintf('\n')]);
        for i = 1:numel(z_range)
            %fprintf(1,['z-Slice ' num2str(i) ' out of ' num2str(numel(z_range)) '.' sprintf('\n')]);
            
            slice_name = [sr.image.fullname 'z' num2str(z_range(i), '%01d') 'c' num2str(sr.image.channel, '%01d') '.tif'];
            FluorImage = imread(slice_name);
            
            for ii = 1:N
                cell = immultiply(LcFull == ii, FluorImage);
                all_cv(i, ii) = sqrt(var(double(cell(LcFull(:) == ii)))) / mean(double(cell(LcFull(:) == ii)));
            end
        end
        
        [~, max_cv_z] = max(all_cv);
        max_FluorImage = zeros(size(LcFull));
        max_normed = zeros(size(LcFull));
        fprintf(1,['Applying best z for each cell.' sprintf('\n')]);
        for i = 1:numel(z_range)
            slice_name = [sr.image.fullname 'z' num2str(z_range(i), '%01d') 'c' num2str(sr.image.channel, '%01d') '.tif'];
            FluorImage = imread(slice_name);          
            
            for ii = 1:N
                if max_cv_z(ii) == i
                max_FluorImage = max_FluorImage + double(immultiply(LcFull == ii, FluorImage));
                
                normed_cell = double(immultiply(LcFull == ii, FluorImage));
                normed_cell = normed_cell - double(min(FluorImage(LcFull == ii)));
                normed_cell = normed_cell / max(normed_cell(:));
                max_normed = max_normed + double(immultiply(LcFull == ii, normed_cell));
                end
            end
        end
        %only cells
        cell_fluor = max_FluorImage(logical(LcFull));
        if norm_over_all
        %normalize only cells and on pixels within the percentiles
            pix_to_analyze = cell_fluor((cell_fluor < prctile(cell_fluor(:),98)) ...
                & (cell_fluor > prctile(cell_fluor(:),2)));
            min_to_sub = min(double(pix_to_analyze(:)));
            pix_to_analyze = double(pix_to_analyze) - min_to_sub;
            max_to_divide = max(pix_to_analyze(:));

            %apply normalization to whole frame
            norm_fluor = double(max_FluorImage) - min_to_sub;
            norm_fluor = norm_fluor/max_to_divide;
            
            %set threshold within normalized range
            high_threshold = .8;%prctile(norm_nobg_fluor,90);
            low_threshold = .4;%prctile(norm_nobg_fluor,50);
        else
            norm_fluor = max_normed;
            %set threshold within normalized range
            high_threshold = .70;%prctile(norm_nobg_fluor,90);
            low_threshold = .4;%prctile(norm_nobg_fluor,50);
        end
        
        %apply threshold
        T_high_image = logical(norm_fluor >= high_threshold);
        T_high_image = T_high_image & LcFull;
        T_low_image = logical(norm_fluor >= low_threshold);
        T_low_image = T_low_image & LcFull;

        smooth_T_high_image = imfill(T_high_image, 'holes');
        se = strel('disk',2);            
        smooth_T_high_image = imerode(smooth_T_high_image,se);
        smooth_T_high_image = imdilate(smooth_T_high_image,se);

        smooth_T_low_image = imfill(T_low_image, 'holes');
        se = strel('disk',2);            
        smooth_T_low_image = imerode(smooth_T_low_image,se);
        smooth_T_low_image = imdilate(smooth_T_low_image,se);

        f1 = figure();
        histogram(norm_fluor(norm_fluor > 0 & norm_fluor < 1));
        xlim([0 1]);

%             nbin = linspace(6, 10 ,100) ;
%             [ya,xa] = hist(log2(double(NoBg_FluorImage(:))),nbin) ;
%             plotDoubleGaussFit_gen(xa(1:end-1) , ya(1:end-1) , nbin) ;
%             ylabel('Probability') ;
%             xlabel('Fluor. intensity within cells') ;
%             legend('show')
        %set(gca, 'YScale', 'log')
        %set(gca, 'XScale', 'log')

        f2 = figure('numbertitle', 'off', 'name' , ['Reconstruction: frame No ' num2str(n_frame)], ...
            'Units', 'Normalized', 'OuterPosition', [0.05 0.05 0.95 0.75]);
        ha = tight_subplot(1, 3, [.01 .01], [.1 .01], [.01 .01]);
        axes(ha(1)); hold on;
        imshow(norm_fluor, [0 1]); %imcontrast
        % cell outline
        spy(ALL_cellPerim,'w',0.5); hold on ;
        axis square; daspect([1 1 1]);

        axes(ha(2)); hold on;
        rgb_image = double(cat(3,T_low_image & ~T_high_image,T_low_image,T_low_image));
        imshow(rgb_image);
        % cell outline
        spy(ALL_cellPerim,'w',0.5); hold on ;
        axis square; daspect([1 1 1]);

        axes(ha(3)); hold on;
        rgb_image = double(cat(3,smooth_T_low_image & ~smooth_T_high_image,smooth_T_low_image,smooth_T_low_image));
        imshow(rgb_image);
        % cell outline
        spy(ALL_cellPerim,'w',0.5); hold on ;
        axis square; daspect([1 1 1]);

        linkaxes(ha, 'xy');

        savefig(f1, [sr.spotrec.output 'clusterdata_histogram' num2str(n_frame,'%03d') '.fig']);
        savefig(f2, [sr.spotrec.output 'clusterdata_reconstruction' num2str(n_frame,'%03d') '.fig']);
        saveas(f1, [sr.spotrec.output 'clusterdata_histogram' num2str(n_frame,'%03d') '.png']);
        saveas(f2, [sr.spotrec.output 'clusterdata_reconstruction' num2str(n_frame,'%03d') '.png']);

        %close all
        
        %%%
        
        Lcluster_low = bwlabel(smooth_T_low_image);
        Lcluster_high = bwlabel(smooth_T_high_image);
        
        %%%
        
        clusterdata_low = zeros(max(Lcluster_low(:)), 10);
        clusterdata_high = zeros(max(Lcluster_high(:)), 10);
        for iii = 1:max(Lcluster_low(:))
            %cell number
            clusterdata_low(iii, 1) = median(LcFull(Lcluster_low == iii));
            %z slice of cell
            clusterdata_low(iii, 2) = max_cv_z(clusterdata_low(iii, 1));
            %intensity total
            clusterdata_low(iii, 3) = sum(max_FluorImage(Lcluster_low == iii));
            %area
            clusterdata_low(iii, 4) = sum(sum(Lcluster_low == iii));
            %intensity per pix
            clusterdata_low(iii, 5) = clusterdata_low(iii, 3) / clusterdata_low(iii, 4);
            %frame number
            clusterdata_low(iii, 6) = n_frame;
            %median background outside cells
            clusterdata_low(iii, 7) = -1;
        end
        
        for iii = 1:max(Lcluster_high(:))
            %cell number
            clusterdata_high(iii, 1) = median(LcFull(Lcluster_high == iii));
            %z slice of cell
            clusterdata_high(iii, 2) = max_cv_z(clusterdata_high(iii, 1));
            %intensity total
            clusterdata_high(iii, 3) = sum(max_FluorImage(Lcluster_high == iii));
            %area
            clusterdata_high(iii, 4) = sum(sum(Lcluster_high == iii));
            %intensity per pix
            clusterdata_high(iii, 5) = clusterdata_high(iii, 3) / clusterdata_high(iii, 4);
            %frame number
            clusterdata_high(iii, 6) = n_frame;
            %median background outside cells
            clusterdata_high(iii, 7) = -1;
            %cluster number
            clusterdata_high(iii, 8) = median(Lcluster_low(Lcluster_high == iii));
        end
        
        save([sr.output 'clusterdata_low' num2str(n_frame,'%03d') '.mat'],'clusterdata_low');
        save([sr.output 'clusterdata_high' num2str(n_frame,'%03d') '.mat'],'clusterdata_high');
    end
end
time = toc;
hours = floor(time/3600);
minutes = floor(time/60) - 60*hours;
seconds = floor(time - 60*minutes - 3600*hours);

progress_ = ['Done. Elapsed time = ' ...
    num2str(hours,'%02d') 'hr:' ...
    num2str(minutes,'%02d') 'min:' ...
    num2str(seconds,'%02d') 'sec'] ;
fprintf(1,[progress_ sprintf('\n')]); 
        
=======
clc
clear

ip = exp_20211029_InitializeExp ; 
% samples channels
%   1   2
%   ph3	TXRED 

spl2cnum = struct() ; 
spl2cnum.Cy5 = 2*ones(1,max(ip.exp.sampleList)) ;
spl2cnum.Cy3 = 3*ones(1,max(ip.exp.sampleList)) ;
spl2cnum.GFP = 4*ones(1,max(ip.exp.sampleList)) ;
spl2cnum.DAPI = 5*ones(1,max(ip.exp.sampleList)) ;


%% test

clc
close all;
clearvars -except spl2cnum ip;

% spot_folder = [ip.exp.path '\cluster_quantify\Run-GFP-.4low.6high-cellnormed-thresh-' date '\'];
norm_over_all = 0;

tic ;
for combo = 3
    if combo == 1
        high_threshold = .43;%prctile(norm_nobg_fluor,90);
        low_threshold = .3;%prctile(norm_nobg_fluor,50);
    elseif combo == 2
        high_threshold = .52;%prctile(norm_nobg_fluor,90);
        low_threshold = .39;%prctile(norm_nobg_fluor,50);
    elseif combo == 3
        high_threshold = .53;%prctile(norm_nobg_fluor,90);
        low_threshold = .4;%prctile(norm_nobg_fluor,50);
    elseif combo == 4
        high_threshold = .54;%prctile(norm_nobg_fluor,90);
        low_threshold = .41;%prctile(norm_nobg_fluor,50);
    elseif combo == 5
        high_threshold = .63;%prctile(norm_nobg_fluor,90);
        low_threshold = .5;%prctile(norm_nobg_fluor,50);
    end
    
    spot_folder = [ip.exp.path '\cluster_quantify\Run-GFP-' num2str(low_threshold) ...
        'low' num2str(high_threshold) 'high-cellnormed-thresh-' date '\'];
    
for n_spl = [1 2 3]
%     if n_spl == 1
%         imgs = [1 12];
%     elseif n_spl == 2
%         imgs = [7];
%     else
%         imgs = [8];
%     end
    
    for n_img =   ip.sample{n_spl}.idx
            
        n_frame = ip.exp.splimg2frm(n_spl , n_img) ;
        time = toc;
        hours = floor(time/3600);
        minutes = floor(time/60) - 60*hours;
        seconds = floor(time - 60*minutes - 3600*hours);

        progress_ = ['SPOT RECOGNITION. FRAME ' num2str(n_frame,'%4d')...
            ' of ' num2str(ip.exp.totalframes,'%4d') '. Elapsed time = ' ...
            num2str(hours,'%02d') 'hr:' ...
            num2str(minutes,'%02d') 'min:' ...
            num2str(seconds,'%02d') 'sec'] ;
        fprintf(1,[progress_ sprintf('\n')]) ;
        n_channel = spl2cnum.GFP(n_spl) ;
        sr = InitializeClusterRecognitionParameters_c234(ip,n_frame,n_channel,spot_folder);
        
        %%%
        
        load([sr.seg.dir sr.seg.name]);
        z_range = sr.image.zrange;
        LcFull = uint16(LcFull);
        ALL_cellPerim = zeros(size(LcFull)) ;
        for iCell = 1:1:max(LcFull(:))
            ALL_cellPerim = ALL_cellPerim + bwperim(LcFull==iCell) ;
        end
        N = max(double(LcFull(:)));
        
        all_cv = zeros(numel(z_range),N);
        fprintf(1,['Finding best z for each cell.' sprintf('\n')]);
        for i = 1:numel(z_range)
            %fprintf(1,['z-Slice ' num2str(i) ' out of ' num2str(numel(z_range)) '.' sprintf('\n')]);
            
            slice_name = [sr.image.fullname 'z' num2str(z_range(i), '%01d') 'c' num2str(sr.image.channel, '%01d') '.tif'];
            FluorImage = imread(slice_name);
            
            for ii = 1:N
                cell = immultiply(LcFull == ii, FluorImage);
                all_cv(i, ii) = sqrt(var(double(cell(LcFull(:) == ii)))) / mean(double(cell(LcFull(:) == ii)));
            end
        end
        
        [~, max_cv_z] = max(all_cv);
        max_FluorImage = zeros(size(LcFull));
        max_normed = zeros(size(LcFull));
        fprintf(1,['Applying best z for each cell.' sprintf('\n')]);
        for i = 1:numel(z_range)
            slice_name = [sr.image.fullname 'z' num2str(z_range(i), '%01d') 'c' num2str(sr.image.channel, '%01d') '.tif'];
            FluorImage = imread(slice_name);          
            
            for ii = 1:N
                if max_cv_z(ii) == i
                max_FluorImage = max_FluorImage + double(immultiply(LcFull == ii, FluorImage));
                
                normed_cell = double(immultiply(LcFull == ii, FluorImage));
                normed_cell = normed_cell - double(min(FluorImage(LcFull == ii)));
                normed_cell = normed_cell / max(normed_cell(:));
                max_normed = max_normed + double(immultiply(LcFull == ii, normed_cell));
                end
            end
        end
        %only cells
        cell_fluor = max_FluorImage(logical(LcFull));
        if norm_over_all
        %normalize only cells and on pixels within the percentiles
            pix_to_analyze = cell_fluor((cell_fluor < prctile(cell_fluor(:),98)) ...
                & (cell_fluor > prctile(cell_fluor(:),2)));
            min_to_sub = min(double(pix_to_analyze(:)));
            pix_to_analyze = double(pix_to_analyze) - min_to_sub;
            max_to_divide = max(pix_to_analyze(:));

            %apply normalization to whole frame
            norm_fluor = double(max_FluorImage) - min_to_sub;
            norm_fluor = norm_fluor/max_to_divide;
            
            %set threshold within normalized range
            high_threshold = .8;%prctile(norm_nobg_fluor,90);
            low_threshold = .4;%prctile(norm_nobg_fluor,50);
        else
            norm_fluor = max_normed;
            %set threshold within normalized range
%             high_threshold = .53;%prctile(norm_nobg_fluor,90);
%             low_threshold = .4;%prctile(norm_nobg_fluor,50);
        end
        
        %apply threshold
        T_high_image = logical(norm_fluor >= high_threshold);
        T_high_image = T_high_image & LcFull;
        T_low_image = logical(norm_fluor >= low_threshold);
        T_low_image = T_low_image & LcFull;

        smooth_T_high_image = imfill(T_high_image, 'holes');
        se = strel('disk',3);            
        smooth_T_high_image = imerode(smooth_T_high_image,se);
        smooth_T_high_image = imdilate(smooth_T_high_image,se);

        smooth_T_low_image = imfill(T_low_image, 'holes');
        se = strel('disk',2);            
        smooth_T_low_image = imerode(smooth_T_low_image,se);
        smooth_T_low_image = imdilate(smooth_T_low_image,se);

        f1 = figure('Units', 'Normalized', 'OuterPosition', [0.05 0.05 0.95/2 0.90]);
%         histogram(norm_fluor(norm_fluor > 0 & norm_fluor < 1));
        subplot(2,1,1);
        xlim([0 1]);
        nbin = linspace(0, 1, 100) ;
        [ya,xa] = hist(norm_fluor(norm_fluor > 0 & norm_fluor < 1),nbin) ;
%             plotDoubleGaussFit_gen(xa(1:end-1) , ya(1:end-1) , nbin) ;
%             plotThreeGaussFit_gen(xa(1:end-1) , ya(1:end-1) , nbin) ;
        plotMultiGaussFit_real(xa(1:end-1) , ya(1:end-1) , nbin) ;
        ylabel('Probability') ;
        xlabel('Fluor. intensity within cells') ;
        legend('show')
        title([ip.sample{n_spl}.name ': Norm. int. in a cell']);
        subplot(2,1,2);
        yyaxis left;
        ylabel('Probability') ; hold on;
        h = bar(xa,ya/sum(ya),'EdgeColor','none','FaceColor', 0.9 * [1 1 1] ,...
            'BarWidth',.9, 'displayname', 'Data histogram'); hold all ; 
        yyaxis right
        ylabel('Fraction of cell');
        colors = {'#F00','#F80','#FF0','#0B0','#00F'};
        nbin = linspace(0, 1 ,51) ;
        for i = 1:5
            cell_norm_fluor = norm_fluor(LcFull == i);
            [N, edges] = histcounts(cell_norm_fluor, 'BinEdges', nbin, 'Normalization', 'probability');
            edges = (edges(1:end-1) + edges(2:end) ) / 2;
            plot(edges, N, 'DisplayName', ['Sample single cell histogram #', num2str(i)], ...
                'LineWidth', 2, 'LineStyle', '-', 'Color', colors{i}, 'Marker', '.');
            hold on;
        end
        legend('Location', 'best');
        xlabel('Norm fluor. intensity within cells') ;
        %set(gca, 'YScale', 'log')
        %set(gca, 'XScale', 'log')

        f2 = figure('numbertitle', 'off', 'name' , ['Reconstruction: frame No ' num2str(n_frame)], ...
            'Units', 'Normalized', 'OuterPosition', [0.05 0.05 0.95 0.75]);
        ha = tight_subplot(1, 3, [.01 .01], [.1 .01], [.01 .01]);
        axes(ha(1)); hold on;
        imshow(norm_fluor, [0 1]); %imcontrast
        % cell outline
        spy(ALL_cellPerim,'w',0.5); hold on ; xlabel('');
        title('Raw colored image');
        axis square; daspect([1 1 1]);

        axes(ha(2)); hold on;
        rgb_image = double(cat(3,T_low_image & ~T_high_image,T_low_image,T_low_image));
        imshow(rgb_image);
        % cell outline
        spy(ALL_cellPerim,'w',0.5); hold on ; xlabel('');
        title('Reconstructed image from clusters');
        axis square; daspect([1 1 1]);

        axes(ha(3)); hold on;
        rgb_image = double(cat(3,smooth_T_low_image & ~smooth_T_high_image,smooth_T_low_image,smooth_T_low_image));
        imshow(rgb_image);
        % cell outline
        spy(ALL_cellPerim,'w',0.5); hold on ; xlabel('');
        title('Cleaned Reconstructed image from clusters');
        axis square; daspect([1 1 1]);

        linkaxes(ha, 'xy');

        savefig(f1, [sr.spotrec.output 'clusterdata_histogram' num2str(n_frame,'%03d') '.fig']);
        savefig(f2, [sr.spotrec.output 'clusterdata_reconstruction' num2str(n_frame,'%03d') '.fig']);
        saveas(f1, [sr.spotrec.output 'clusterdata_histogram' num2str(n_frame,'%03d') '.png']);
        saveas(f2, [sr.spotrec.output 'clusterdata_reconstruction' num2str(n_frame,'%03d') '.png']);

        close all
        
        %%%
        
        Lcluster_low = bwlabel(smooth_T_low_image);
        Lcluster_low = Lcluster_low .* (Lcluster_low & ~smooth_T_high_image);
        if exist([sr.output 'clusterdata_low' num2str(n_frame,'%03d') '.mat'],'file')
            save([sr.output 'clusterdata_low' num2str(n_frame,'%03d') '.mat'],'Lcluster_low', '-append');
        else
            save([sr.output 'clusterdata_low' num2str(n_frame,'%03d') '.mat'],'Lcluster_low');
        end
        
        Lcluster_high = bwlabel(smooth_T_high_image);
        if exist([sr.output 'clusterdata_high' num2str(n_frame,'%03d') '.mat'],'file')
            save([sr.output 'clusterdata_high' num2str(n_frame,'%03d') '.mat'],'Lcluster_high', '-append');
        else
            save([sr.output 'clusterdata_high' num2str(n_frame,'%03d') '.mat'],'Lcluster_high');
        end
        
        %%%
        
        clusterdata_low = zeros(max(Lcluster_low(:)), 10);
        clusterdata_high = zeros(max(Lcluster_high(:)), 10);
        for iii = 1:max(Lcluster_low(:))
            %cell number
            clusterdata_low(iii, 1) = median(LcFull(Lcluster_low == iii));
            %z slice of cell
            clusterdata_low(iii, 2) = max_cv_z(clusterdata_low(iii, 1));
            %intensity total
            clusterdata_low(iii, 3) = sum(max_FluorImage(Lcluster_low == iii));
            %area
            clusterdata_low(iii, 4) = sum(sum(Lcluster_low == iii));
            %intensity per pix
            clusterdata_low(iii, 5) = clusterdata_low(iii, 3) / clusterdata_low(iii, 4);
            %frame number
            clusterdata_low(iii, 6) = n_frame;
            %X pixel position
            x_coords = repmat(1:ip.image.size, [ip.image.size, 1]);
            clusterdata_low(iii, 7) = mean( x_coords(Lcluster_low == iii) );
            %Y pixel position
            y_coords = repmat([1:ip.image.size]', [1, ip.image.size]);
            clusterdata_low(iii, 8) = mean( y_coords(Lcluster_low == iii) );
            %cluster number
            clusterdata_low(iii, 9) = iii;
        end
        
        for iii = 1:max(Lcluster_high(:))
            %cell number
            clusterdata_high(iii, 1) = median(LcFull(Lcluster_high == iii));
            %z slice of cell
            clusterdata_high(iii, 2) = max_cv_z(clusterdata_high(iii, 1));
            %intensity total
            clusterdata_high(iii, 3) = sum(max_FluorImage(Lcluster_high == iii));
            %area
            clusterdata_high(iii, 4) = sum(sum(Lcluster_high == iii));
            %intensity per pix
            clusterdata_high(iii, 5) = clusterdata_high(iii, 3) / clusterdata_high(iii, 4);
            %frame number
            clusterdata_high(iii, 6) = n_frame;
            %X pixel position
            x_coords = repmat(1:ip.image.size, [ip.image.size, 1]);
            clusterdata_high(iii, 7) = mean( x_coords(Lcluster_high == iii) );
            %Y pixel position
            y_coords = repmat([1:ip.image.size]', [1, ip.image.size]);
            clusterdata_high(iii, 8) = mean( y_coords(Lcluster_high == iii) );
            %low cluster number
            clusterdata_high(iii, 9) = median(Lcluster_low(Lcluster_high == iii));
            %high cluster number
            clusterdata_high(iii, 10) = iii;
        end
        
        if exist([sr.output 'clusterdata_low' num2str(n_frame,'%03d') '.mat'],'file')
            save([sr.output 'clusterdata_low' num2str(n_frame,'%03d') '.mat'],'clusterdata_low','-append');
        else
            save([sr.output 'clusterdata_low' num2str(n_frame,'%03d') '.mat'],'clusterdata_low');
        end
        
        if exist([sr.output 'clusterdata_high' num2str(n_frame,'%03d') '.mat'],'file')
            save([sr.output 'clusterdata_high' num2str(n_frame,'%03d') '.mat'],'clusterdata_high','-append');
        else
            save([sr.output 'clusterdata_high' num2str(n_frame,'%03d') '.mat'],'clusterdata_high');
        end
    end
end
time = toc;
hours = floor(time/3600);
minutes = floor(time/60) - 60*hours;
seconds = floor(time - 60*minutes - 3600*hours);

progress_ = ['Done. Elapsed time = ' ...
    num2str(hours,'%02d') 'hr:' ...
    num2str(minutes,'%02d') 'min:' ...
    num2str(seconds,'%02d') 'sec'] ;
fprintf(1,[progress_ sprintf('\n')]);      
end


%% run

clc
close all;
clearvars -except spl2cnum ip;

spot_folder = [ip.exp.path '\cluster_quantify\Run-TXRED-.4low.7high-cellnormed-thresh-' date '\'];
norm_over_all = 0;

tic ;

for n_spl = [10]%ip.exp.sampleList
    for n_img =   ip.sample{n_spl}.idx
        
        close all;
        
        n_frame = ip.exp.splimg2frm(n_spl , n_img) ;
        time = toc;
        hours = floor(time/3600);
        minutes = floor(time/60) - 60*hours;
        seconds = floor(time - 60*minutes - 3600*hours);

        progress_ = ['SPOT RECOGNITION. FRAME ' num2str(n_frame,'%4d')...
            ' of ' num2str(ip.exp.totalframes,'%4d') '. Elapsed time = ' ...
            num2str(hours,'%02d') 'hr:' ...
            num2str(minutes,'%02d') 'min:' ...
            num2str(seconds,'%02d') 'sec'] ;
        fprintf(1,[progress_ sprintf('\n')]) ;
        n_channel = spl2cnum.TXRED(n_spl) ;
        sr = InitializeClusterRecognitionParameters_c234(ip,n_frame,n_channel,spot_folder);
        
        %%%
        
        load([sr.seg.dir sr.seg.name]);
        z_range = sr.image.zrange;
        LcFull = uint16(LcFull);
        ALL_cellPerim = zeros(size(LcFull)) ;
        for iCell = 1:1:max(LcFull(:))
            ALL_cellPerim = ALL_cellPerim + bwperim(LcFull==iCell) ;
        end
        N = max(double(LcFull(:)));
        
        all_cv = zeros(numel(z_range),N);
        fprintf(1,['Finding best z for each cell.' sprintf('\n')]);
        for i = 1:numel(z_range)
            %fprintf(1,['z-Slice ' num2str(i) ' out of ' num2str(numel(z_range)) '.' sprintf('\n')]);
            
            slice_name = [sr.image.fullname 'z' num2str(z_range(i), '%01d') 'c' num2str(sr.image.channel, '%01d') '.tif'];
            FluorImage = imread(slice_name);
            
            for ii = 1:N
                cell = immultiply(LcFull == ii, FluorImage);
                all_cv(i, ii) = sqrt(var(double(cell(LcFull(:) == ii)))) / mean(double(cell(LcFull(:) == ii)));
            end
        end
        
        [~, max_cv_z] = max(all_cv);
        max_FluorImage = zeros(size(LcFull));
        max_normed = zeros(size(LcFull));
        fprintf(1,['Applying best z for each cell.' sprintf('\n')]);
        for i = 1:numel(z_range)
            slice_name = [sr.image.fullname 'z' num2str(z_range(i), '%01d') 'c' num2str(sr.image.channel, '%01d') '.tif'];
            FluorImage = imread(slice_name);          
            
            for ii = 1:N
                if max_cv_z(ii) == i
                max_FluorImage = max_FluorImage + double(immultiply(LcFull == ii, FluorImage));
                
                normed_cell = double(immultiply(LcFull == ii, FluorImage));
                normed_cell = normed_cell - double(min(FluorImage(LcFull == ii)));
                normed_cell = normed_cell / max(normed_cell(:));
                max_normed = max_normed + double(immultiply(LcFull == ii, normed_cell));
                end
            end
        end
        %only cells
        cell_fluor = max_FluorImage(logical(LcFull));
        if norm_over_all
        %normalize only cells and on pixels within the percentiles
            pix_to_analyze = cell_fluor((cell_fluor < prctile(cell_fluor(:),98)) ...
                & (cell_fluor > prctile(cell_fluor(:),2)));
            min_to_sub = min(double(pix_to_analyze(:)));
            pix_to_analyze = double(pix_to_analyze) - min_to_sub;
            max_to_divide = max(pix_to_analyze(:));

            %apply normalization to whole frame
            norm_fluor = double(max_FluorImage) - min_to_sub;
            norm_fluor = norm_fluor/max_to_divide;
            
            %set threshold within normalized range
            high_threshold = .8;%prctile(norm_nobg_fluor,90);
            low_threshold = .4;%prctile(norm_nobg_fluor,50);
        else
            norm_fluor = max_normed;
            %set threshold within normalized range
            high_threshold = .70;%prctile(norm_nobg_fluor,90);
            low_threshold = .4;%prctile(norm_nobg_fluor,50);
        end
        
        %apply threshold
        T_high_image = logical(norm_fluor >= high_threshold);
        T_high_image = T_high_image & LcFull;
        T_low_image = logical(norm_fluor >= low_threshold);
        T_low_image = T_low_image & LcFull;

        smooth_T_high_image = imfill(T_high_image, 'holes');
        se = strel('disk',2);            
        smooth_T_high_image = imerode(smooth_T_high_image,se);
        smooth_T_high_image = imdilate(smooth_T_high_image,se);

        smooth_T_low_image = imfill(T_low_image, 'holes');
        se = strel('disk',2);            
        smooth_T_low_image = imerode(smooth_T_low_image,se);
        smooth_T_low_image = imdilate(smooth_T_low_image,se);

        f1 = figure();
        histogram(norm_fluor(norm_fluor > 0 & norm_fluor < 1));
        xlim([0 1]);

%             nbin = linspace(6, 10 ,100) ;
%             [ya,xa] = hist(log2(double(NoBg_FluorImage(:))),nbin) ;
%             plotDoubleGaussFit_gen(xa(1:end-1) , ya(1:end-1) , nbin) ;
%             ylabel('Probability') ;
%             xlabel('Fluor. intensity within cells') ;
%             legend('show')
        %set(gca, 'YScale', 'log')
        %set(gca, 'XScale', 'log')

        f2 = figure('numbertitle', 'off', 'name' , ['Reconstruction: frame No ' num2str(n_frame)], ...
            'Units', 'Normalized', 'OuterPosition', [0.05 0.05 0.95 0.75]);
        ha = tight_subplot(1, 3, [.01 .01], [.1 .01], [.01 .01]);
        axes(ha(1)); hold on;
        imshow(norm_fluor, [0 1]); %imcontrast
        % cell outline
        spy(ALL_cellPerim,'w',0.5); hold on ;
        axis square; daspect([1 1 1]);

        axes(ha(2)); hold on;
        rgb_image = double(cat(3,T_low_image & ~T_high_image,T_low_image,T_low_image));
        imshow(rgb_image);
        % cell outline
        spy(ALL_cellPerim,'w',0.5); hold on ;
        axis square; daspect([1 1 1]);

        axes(ha(3)); hold on;
        rgb_image = double(cat(3,smooth_T_low_image & ~smooth_T_high_image,smooth_T_low_image,smooth_T_low_image));
        imshow(rgb_image);
        % cell outline
        spy(ALL_cellPerim,'w',0.5); hold on ;
        axis square; daspect([1 1 1]);

        linkaxes(ha, 'xy');

        savefig(f1, [sr.spotrec.output 'clusterdata_histogram' num2str(n_frame,'%03d') '.fig']);
        savefig(f2, [sr.spotrec.output 'clusterdata_reconstruction' num2str(n_frame,'%03d') '.fig']);
        saveas(f1, [sr.spotrec.output 'clusterdata_histogram' num2str(n_frame,'%03d') '.png']);
        saveas(f2, [sr.spotrec.output 'clusterdata_reconstruction' num2str(n_frame,'%03d') '.png']);

        %close all
        
        %%%
        
        Lcluster_low = bwlabel(smooth_T_low_image);
        Lcluster_high = bwlabel(smooth_T_high_image);
        
        %%%
        
        clusterdata_low = zeros(max(Lcluster_low(:)), 10);
        clusterdata_high = zeros(max(Lcluster_high(:)), 10);
        for iii = 1:max(Lcluster_low(:))
            %cell number
            clusterdata_low(iii, 1) = median(LcFull(Lcluster_low == iii));
            %z slice of cell
            clusterdata_low(iii, 2) = max_cv_z(clusterdata_low(iii, 1));
            %intensity total
            clusterdata_low(iii, 3) = sum(max_FluorImage(Lcluster_low == iii));
            %area
            clusterdata_low(iii, 4) = sum(sum(Lcluster_low == iii));
            %intensity per pix
            clusterdata_low(iii, 5) = clusterdata_low(iii, 3) / clusterdata_low(iii, 4);
            %frame number
            clusterdata_low(iii, 6) = n_frame;
            %median background outside cells
            clusterdata_low(iii, 7) = -1;
        end
        
        for iii = 1:max(Lcluster_high(:))
            %cell number
            clusterdata_high(iii, 1) = median(LcFull(Lcluster_high == iii));
            %z slice of cell
            clusterdata_high(iii, 2) = max_cv_z(clusterdata_high(iii, 1));
            %intensity total
            clusterdata_high(iii, 3) = sum(max_FluorImage(Lcluster_high == iii));
            %area
            clusterdata_high(iii, 4) = sum(sum(Lcluster_high == iii));
            %intensity per pix
            clusterdata_high(iii, 5) = clusterdata_high(iii, 3) / clusterdata_high(iii, 4);
            %frame number
            clusterdata_high(iii, 6) = n_frame;
            %median background outside cells
            clusterdata_high(iii, 7) = -1;
            %cluster number
            clusterdata_high(iii, 8) = median(Lcluster_low(Lcluster_high == iii));
        end
        
        save([sr.output 'clusterdata_low' num2str(n_frame,'%03d') '.mat'],'clusterdata_low');
        save([sr.output 'clusterdata_high' num2str(n_frame,'%03d') '.mat'],'clusterdata_high');
    end
end
time = toc;
hours = floor(time/3600);
minutes = floor(time/60) - 60*hours;
seconds = floor(time - 60*minutes - 3600*hours);

progress_ = ['Done. Elapsed time = ' ...
    num2str(hours,'%02d') 'hr:' ...
    num2str(minutes,'%02d') 'min:' ...
    num2str(seconds,'%02d') 'sec'] ;
fprintf(1,[progress_ sprintf('\n')]); 
        
>>>>>>> 88dd4f4ac7eb6f5ae758552b04c47df241263be0
