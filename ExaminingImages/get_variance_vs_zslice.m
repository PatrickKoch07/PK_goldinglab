<<<<<<< HEAD
%% Initialize
clc
clear all;

ip = exp_20211029_InitializeExp();

channels = ip.image.channels.name;
date = ['(' ip.exp.date(6:7) '/' ip.exp.date(9:10) ')'];

if ~exist([ip.exp.path '\presegmentation_fluor_analysis'], 'dir')
    mkdir([ip.exp.path '\presegmentation_fluor_analysis']);
end
if ~exist([ip.exp.path '\postsegmentation_fluor_analysis'], 'dir')
    mkdir([ip.exp.path '\postsegmentation_fluor_analysis']);
end

for c = channels
    if ~exist([ip.exp.path '\presegmentation_fluor_analysis\' c{1} '_CV_var_int_vs_zslice'], 'dir')
        mkdir([ip.exp.path '\presegmentation_fluor_analysis\' c{1} '_CV_var_int_vs_zslice']);
    end
    
    if ~exist([ip.exp.path '\postsegmentation_fluor_analysis\' c{1} '_CV_var_int_vs_zslice'], 'dir')
        mkdir([ip.exp.path '\postsegmentation_fluor_analysis\' c{1} '_CV_var_int_vs_zslice']);
    end  
end

clear c


%% Variance & Intensity & CV vs. z-slice PRE
fluor_chan = 2;

phase3_bestz = zeros(ip.exp.totalframes, 1);
fluor_bestvarz = zeros(ip.exp.totalframes,1);
fluor_bestCVz = zeros(ip.exp.totalframes,1);

color_used = {};

for c = 1:length(channels)
    %%% CV vs. z-slice PRE
    
    spl2znum = ip.image.zrange .* ...
        ones(1,max(ip.exp.sampleList)) ;

    y_min(1) = 1;
    y_min(2) = 1;
    y_min(3) = 1;

    for n_frame = 1:ip.exp.totalframes
        % default value
        imdigitnum_ = 2 ;
        current_spl = ip.exp.frm2spl(n_frame) ;
        if max(ip.sample{current_spl}.idx) <= 9
            imdigitnum_ = 1 ;
        end
        var_lst_cv = zeros(1,ip.image.zrange) ;
        var_lst_var = zeros(1,ip.image.zrange) ;
        var_lst_int = zeros(1,ip.image.zrange) ;

        for iz = 1: spl2znum(current_spl)
            PH3 = imread([ip.exp.path 'images\' ...
                ip.image.base_name '_' ...
                num2str(ip.exp.frm2spl(n_frame),'%03d') ...
                'xy' num2str(ip.exp.frm2img(n_frame),['%0' num2str(imdigitnum_ ,'%01d') 'd']) ...
                ['z' num2str(iz ,'%01d') 'c' num2str(c) '.tif'] ]);
            var_lst_cv(iz) = sqrt(var(double(PH3(:)))) / mean(double(PH3(:))) ;
            var_lst_var(iz) = var(double(PH3(:))) ;
            var_lst_int(iz) = sum(double(PH3(:))) ;
        end

        if y_min(1) > min(var_lst_cv/max(var_lst_cv))
            y_min(1) = min(var_lst_cv/max(var_lst_cv));
        end
        
        if y_min(2) > min(var_lst_var/max(var_lst_var))
            y_min(2) = min(var_lst_var/max(var_lst_var));
        end
        
        if y_min(3) > min(var_lst_int/max(var_lst_int))
            y_min(3) = min(var_lst_int/max(var_lst_int));
        end
        
        if c == fluor_chan
            [~,fluor_bestCVz(n_frame)] = max(var_lst_cv);
            [~,fluor_bestvarz(n_frame)] = max(var_lst_var);
        end
        
        if strcmp(channels{c}, 'phase3')
            [~,phase3_bestz(n_frame)] = max(var_lst_var);
        end

        figure(ip.exp.frm2spl(n_frame) + (c-1)*max(ip.exp.frm2spl));
        subplot(1,3,1);
        hold on;
        if c == 1
            color_used{n_frame} = which_color(gca);
        end
        plot(1: spl2znum(current_spl), var_lst_cv/max(var_lst_cv), 'DisplayName', ...
            ['XY position ', num2str(ip.exp.frm2img(n_frame))], 'Color', color_used{n_frame});
        xlabel('z slice')
        ylabel([channels{c} ' pixel CV / max ' channels{c} ' CV']);
        title(['Normalized ' channels{c} ' pixel CV vs. z-slice'])
        
        subplot(1,3,2);
        hold on;
        %color_used{n_frame} = which_color(gca);
        plot(1: spl2znum(current_spl), var_lst_var/max(var_lst_var), 'DisplayName', ...
            ['XY position ', num2str(ip.exp.frm2img(n_frame))], 'Color', color_used{n_frame});
        xlabel('z slice')
        ylabel([channels{c} ' pixel variance / max ' channels{c} ' variance']);
        title(['Normalized ' channels{c} ' pixel var. vs. z-slice'])
        
        subplot(1,3,3);
        hold on;
        %color_used{n_frame} = which_color(gca);
        plot(1: spl2znum(current_spl), var_lst_int/max(var_lst_int), 'DisplayName', ...
            ['XY position ', num2str(ip.exp.frm2img(n_frame))], 'Color', color_used{n_frame});
        xlabel('z slice')
        ylabel([channels{c} ' pixel intensity / max ' channels{c} ' intensity']);
        title(['Normalized ' channels{c} ' pixel int. vs. z-slice'])        
    end
    %{
    %%% Variance vs. z-slice PRE
    
    spl2znum = ip.image.zrange .* ...
        ones(1,max(ip.exp.sampleList)) ;

    color_used = {};
    y_min(2) = 1;

    for n_frame = 1:ip.exp.totalframes
        % default value
        imdigitnum_ = 2 ;
        current_spl = ip.exp.frm2spl(n_frame) ;
        if max(ip.sample{current_spl}.idx) <= 9
            imdigitnum_ = 1 ;
        end
        var_lst = zeros(1,ip.image.zrange) ;

        for iz = 1: spl2znum(current_spl)
            PH3 = imread([ip.exp.path 'images\' ...
                ip.image.base_name '_' ...
                num2str(ip.exp.frm2spl(n_frame),'%03d') ...
                'xy' num2str(ip.exp.frm2img(n_frame),['%0' num2str(imdigitnum_ ,'%01d') 'd']) ...
                ['z' num2str(iz ,'%01d') 'c' num2str(c) '.tif'] ]);
            var_lst(iz) = var(double(PH3(:))) ;
        end

        if y_min(2) > min(var_lst/max(var_lst))
            y_min(2) = min(var_lst/max(var_lst));
        end
        
        if c == fluor_chan
            [~,fluor_bestvarz(n_frame)] = max(var_lst);
        end
        
        if strcmp(channels{c}, 'phase3')
            [~,phase3_bestz(n_frame)] = max(var_lst);
        end

        figure(ip.exp.frm2spl(n_frame) + (c-1)*max(ip.exp.frm2spl));
        subplot(1,3,2);
        hold on;
        color_used{n_frame} = which_color(gca);
        plot(1: spl2znum(current_spl), var_lst/max(var_lst), 'DisplayName', ['XY position ', num2str(ip.exp.frm2img(n_frame))], 'Color', color_used{n_frame});
        xlabel('z slice')
        ylabel([channels{c} ' pixel variance / max ' channels{c} ' variance']);
        title(['Normalized ' channels{c} ' pixel var. vs. z-slice'])
    end

    %%% Intensity vs. z-slice PRE

    spl2znum = ip.image.zrange .* ...
        ones(1,max(ip.exp.sampleList)) ;

    color_used = {};
    y_min(3) = 1;

    for n_frame = 1:ip.exp.totalframes

        % default value
        imdigitnum_ = 2 ;
        current_spl = ip.exp.frm2spl(n_frame) ;
        if max(ip.sample{current_spl}.idx) <= 9
            imdigitnum_ = 1 ;
        end
        int_lst = zeros(1,ip.image.zrange) ;

        for iz = 1: spl2znum(current_spl)
            PH3 = imread([ip.exp.path 'images\' ...
                ip.image.base_name '_' ...
                num2str(ip.exp.frm2spl(n_frame),'%03d') ...
                'xy' num2str(ip.exp.frm2img(n_frame),['%0' num2str(imdigitnum_ ,'%01d') 'd']) ...
                ['z' num2str(iz ,'%01d') 'c' num2str(c) '.tif'] ]);
            int_lst(iz) = sum(double(PH3(:))) ;
        end

        if y_min(3) > min(int_lst/max(int_lst))
            y_min(3) = min(int_lst/max(int_lst));
        end

        figure(ip.exp.frm2spl(n_frame) + (c-1)*max(ip.exp.frm2spl));
        subplot(1,3,3);
        hold on;
        color_used{n_frame} = which_color(gca);
        plot(1: spl2znum(current_spl), int_lst/max(int_lst), 'DisplayName', ['XY position ', num2str(ip.exp.frm2img(n_frame))], 'Color', color_used{n_frame});
        xlabel('z slice')
        ylabel([channels{c} ' pixel intensity / max ' channels{c} ' intensity']);
        title(['Normalized ' channels{c} ' pixel int. vs. z-slice'])
    end
    %}
    %%% Saving figures of current channel
    
    range_to_process = unique(ip.exp.frm2spl(1:ip.exp.totalframes));
    for i = transpose(range_to_process)
        figure(i + (c-1)*max(ip.exp.frm2spl));
        set(gcf,'units','points','position',[-1000,250,720,326]);
        sgtitle([date ' Sample ', num2str(i), ': ',ip.sample{1,i}.name]);
        for ii = 1:3
            subplot(1,3,ii);
            legend('location', 'best');
            ylim([y_min(ii) 1]);
            
            grid on;
        end

        saveas(figure(i + (c-1)*max(ip.exp.frm2spl)), ['presegmentation_fluor_analysis\' channels{c} '_CV_var_int_vs_zslice\Sample_', ...
            num2str(i, '%02d'), '__Norm' channels{c} 'Pixel__CV_var_int_vs_z.png']);
        saveas(figure(i + (c-1)*max(ip.exp.frm2spl)), ['presegmentation_fluor_analysis\' channels{c} '_CV_var_int_vs_zslice\Sample_', ...
            num2str(i, '%02d'), '__Norm' channels{c} 'Pixel__CV_var_int_vs_z.fig']);
        close(figure(i + (c-1)*max(ip.exp.frm2spl)));
    end
end

%%% Best phase3 variance vs. best FLUOR variance
cmap = hsv(max(ip.exp.frm2spl));

figure ();
subplot(2,2,1);
for i = 1:max(ip.exp.frm2spl)
    xy = [phase3_bestz(nonzeros(ip.exp.splimg2frm(i,:))) - max(ip.exp.frm2spl)/2 * .05 + i * .05,...
        fluor_bestvarz(nonzeros(ip.exp.splimg2frm(i,:))) - max(ip.exp.frm2spl)/2 * .05 + i * .05];
    [C,ia,ic] = unique(xy, 'rows');
    a_counts = accumarray(ic,1);

    scatter(C(:,1), C(:,2), a_counts*5, 'markeredgecolor', 'none',...
        'markerfacecolor', cmap(i,:), 'DisplayName', ['Sample ', num2str(i)]);

%     scatter(unique(phase3_bestz(nonzeros(ip.exp.splimg2frm(i,:))) - max(ip.exp.frm2spl)/2 * .05 + i * .05), ...
%         unique(fluor_bestvarz(nonzeros(ip.exp.splimg2frm(i,:))) - max(ip.exp.frm2spl)/2 * .05 + i * .05), ...
%         a_counts/max(a_counts) * 15, 'markerfacecolor', cmap(i,:), 'DisplayName', ['Sample ', num2str(i)]);
    hold on;
end
plot(0:8, 0:8, 'Color', 'black', 'DisplayName', 'Slope = 1');
xlim([0 8]);
ylim([0 8]);
title(['Z w/ max Ph3 var. vs. Z w/ max ' channels{fluor_chan} ' var.']);
xlabel('Phase 3 highest variance z slice');
ylabel([channels(fluor_chan) ' highest variance z slice']);
grid on;
legend('location', 'best');

subplot(2,2,3);
xy = [phase3_bestz(nonzeros(ip.exp.splimg2frm(:))), fluor_bestvarz(nonzeros(ip.exp.splimg2frm(:)))];
[C,ia,ic] = unique(xy, 'rows');
a_counts = accumarray(ic,1);
scatter(C(:,1), C(:,2), a_counts*5, 'markeredgecolor', 'none',...
    'markerfacecolor', 'black', 'DisplayName', ['All samples together']);
hold on;
plot(0:8, 0:8, 'Color', 'black', 'DisplayName', 'Slope = 1');
xlim([0 8]);
ylim([0 8]);
title(['Z w/ max Ph3 var. vs. Z w/ max ' channels{fluor_chan} ' var.']);
xlabel('Phase 3 highest variance z slice');
ylabel([channels(fluor_chan) ' highest variance z slice']);
grid on;
legend('location', 'best');


subplot(2,2,2);
for i = 1:max(ip.exp.frm2spl)
    xy = [phase3_bestz(nonzeros(ip.exp.splimg2frm(i,:))) - max(ip.exp.frm2spl)/2 * .05 + i * .05,...
        fluor_bestCVz(nonzeros(ip.exp.splimg2frm(i,:))) - max(ip.exp.frm2spl)/2 * .05 + i * .05];
    [C,ia,ic] = unique(xy, 'rows');
    a_counts = accumarray(ic,1);

    scatter(C(:,1), C(:,2), a_counts/max(a_counts) * 50, 'markeredgecolor', 'none',...
        'markerfacecolor', cmap(i,:), 'DisplayName', ['Sample ', num2str(i)]);
    
%     scatter(phase3_bestz(nonzeros(ip.exp.splimg2frm(i,:))) - max(ip.exp.frm2spl)/2 * .05 + i * .05, ...
%         fluor_bestCVz(nonzeros(ip.exp.splimg2frm(i,:))) - max(ip.exp.frm2spl)/2 * .05 + i * .05, ...
%         a_counts/max(a_counts) * 15, 'markerfacecolor', cmap(i,:), 'DisplayName', ['Sample ', num2str(i)]);
    hold on;
end
plot(0:8, 0:8, 'Color', 'black', 'DisplayName', 'Slope = 1');
xlim([0 8]);
ylim([0 8]);
title(['Z w/ max Ph3 CV vs. Z w/ max ' channels{fluor_chan} ' CV']);
xlabel('Phase 3 highest CV z slice');
ylabel([channels{fluor_chan} ' highest CV z slice']);
grid on;
legend('location', 'best');

subplot(2,2,4);
xy = [phase3_bestz(nonzeros(ip.exp.splimg2frm(:))), fluor_bestCVz(nonzeros(ip.exp.splimg2frm(:)))];
[C,ia,ic] = unique(xy, 'rows');
a_counts = accumarray(ic,1);
scatter(C(:,1), C(:,2), a_counts*5, 'markeredgecolor', 'none',...
    'markerfacecolor', 'black', 'DisplayName', ['All samples together']);
hold on;
plot(0:8, 0:8, 'Color', 'black', 'DisplayName', 'Slope = 1');
xlim([0 8]);
ylim([0 8]);
title(['Z w/ max Ph3 CV vs. Z w/ max ' channels{fluor_chan} ' CV']);
xlabel('Phase 3 highest CV z slice');
ylabel([channels{fluor_chan} ' highest CV z slice']);
grid on;
legend('location', 'best');

set(gcf,'units','points','position',[-1200,50,720*1.5,275*1.25*2]);

saveas(figure(i + (c-1)*max(ip.exp.frm2spl)), ['presegmentation_fluor_analysis\phase3_maxvar_z_vs_' ...
    channels{fluor_chan} '_maxvar_z_maxcv_z.png']);
saveas(figure(i + (c-1)*max(ip.exp.frm2spl)), ['presegmentation_fluor_analysis\phase3_maxvar_z_vs_' ...
    channels{fluor_chan} '_maxvar_z_maxcv_z.fig']);


%% Variance & Intensity & CV vs. z-slice POST (whole frame focus comparison)
fluor_chan = 2;
neg_samp = 8;
wholexy = 0;

phase3_bestz = zeros(ip.exp.totalframes, 1);
phase3_bestzPRE = zeros(ip.exp.totalframes, 1);
fluor_bestvarz = zeros(ip.exp.totalframes,1);
fluor_bestCVz = zeros(ip.exp.totalframes,1);

[num_sam, num_img] = size(ip.exp.splimg2frm);

fluor_varz = zeros([num_sam - 1, num_img, ip.image.zrange]);
fluor_CVz = zeros([num_sam - 1, num_img, ip.image.zrange]);

color_used = {};

for c = 1:length(channels)
    %%% CV vs. z-slice PRE
    
    spl2znum = ip.image.zrange .* ...
        ones(1,max(ip.exp.sampleList)) ;

    y_min(1) = 1;
    y_min(2) = 1;
    y_min(3) = 1;
    
    for n_frame = 1:ip.exp.totalframes
        n_image = ip.exp.frm2img(n_frame);
        %figure(ip.exp.frm2spl(n_frame));    
        %set(gcf,'units','points','position',[-1000,31,720,715]);
        
        if n_image == 1
            fprintf(1, newline);
            fprintf(1,['Sample ' num2str(ip.exp.frm2spl(n_frame)) ' of ' ...
                num2str(max(ip.exp.sampleList)) '.' newline]);
            progress_2 = [];
        end
        for d_=1:1:size(progress_2,2) ; fprintf(1,'\b') ; end
        progress_2 = [sprintf('\t') 'Image ' num2str(n_image) ' in ' ...
            num2str(min(ip.sample{ip.exp.frm2spl(n_frame)}.idx)) '-' ...
            num2str(max(ip.sample{ip.exp.frm2spl(n_frame)}.idx)) '.'] ;
        fprintf(1,progress_2) ;

        load([ip.seg.dir ip.seg.base_name 'seg' num2str(n_frame,'%03d') '.mat'],'LcFull') ;

        % default value
        imdigitnum_ = 2 ;
        current_spl = ip.exp.frm2spl(n_frame) ;
        if max(ip.sample{current_spl}.idx) <= 9
            imdigitnum_ = 1 ;
        end
        var_lst = zeros(1,ip.image.zrange) ;

        CellMap = LcFull;
        nnn=double(LcFull);
        N=max(nnn(:));
        
        cv_temp = zeros(spl2znum(current_spl),N);
        var_temp = zeros(spl2znum(current_spl),N);
        int_temp = zeros(spl2znum(current_spl),N);
        
        for index=1:N
            for iz = 1: spl2znum(current_spl)
                PH3 = imread([ip.exp.path 'images\' ...
                    ip.image.base_name '_' ...
                    num2str(ip.exp.frm2spl(n_frame),'%03d') ...
                    'xy' num2str(ip.exp.frm2img(n_frame),['%0' num2str(imdigitnum_ ,'%01d') 'd']) ...
                    ['z' num2str(iz ,'%01d') 'c' num2str(c) '.tif'] ]);

                Cell = immultiply( CellMap == index, PH3);
                int_per_pix = Cell(CellMap == index);

                cv_temp(iz,index) = sqrt(var(double(int_per_pix(:)))) / mean(double(int_per_pix(:))) ;
                var_temp(iz,index) = var(double(int_per_pix(:))) ;
                int_temp(iz,index) = sum(double(int_per_pix(:))) ;
            end
        end
        
        % sum whole image, keep z seperate
        var_lst_cv = sum(transpose(cv_temp)) ;
        var_lst_var = sum(transpose(var_temp)) ;
        var_lst_int = sum(transpose(int_temp)) ;
    
        %for plot limits
        if y_min(1) > min(var_lst_cv/max(var_lst_cv))
            y_min(1) = min(var_lst_cv/max(var_lst_cv));
        end
        
        if y_min(2) > min(var_lst_var/max(var_lst_var))
            y_min(2) = min(var_lst_var/max(var_lst_var));
        end
        
        if y_min(3) > min(var_lst_int/max(var_lst_int))
            y_min(3) = min(var_lst_int/max(var_lst_int));
        end
        
        %for specific channel of interest only, get var and cv normalized
        %for bar plots later
        if c == fluor_chan
            [~,fluor_bestCVz(n_frame)] = max(var_lst_cv);
            [~,fluor_bestvarz(n_frame)] = max(var_lst_var);
            if ip.exp.frm2spl(n_frame) < neg_samp
                fluor_varz(ip.exp.frm2spl(n_frame),ip.exp.frm2img(n_frame), :) = var_lst_var/max(var_lst_var);
                fluor_CVz(ip.exp.frm2spl(n_frame),ip.exp.frm2img(n_frame), :) = var_lst_cv/max(var_lst_cv);
            elseif ip.exp.frm2spl(n_frame) > neg_samp
                fluor_varz(ip.exp.frm2spl(n_frame)-1,ip.exp.frm2img(n_frame), :) = var_lst_var/max(var_lst_var);
                fluor_CVz(ip.exp.frm2spl(n_frame)-1,ip.exp.frm2img(n_frame), :) = var_lst_cv/max(var_lst_cv);
            end
        end
        %go over the whole image only (rather than cell by cell)
        if strcmp(channels{c}, 'phase3')
            [~,phase3_bestz(n_frame)] = max(var_lst_var);
            
            var_temp = zeros(spl2znum(current_spl),1);
            for iz = 1: spl2znum(current_spl)
                PH3 = imread([ip.exp.path 'images\' ...
                    ip.image.base_name '_' ...
                    num2str(ip.exp.frm2spl(n_frame),'%03d') ...
                    'xy' num2str(ip.exp.frm2img(n_frame),['%0' num2str(imdigitnum_ ,'%01d') 'd']) ...
                    ['z' num2str(iz ,'%01d') 'c' num2str(c) '.tif'] ]);
                
                cv_temp(iz) = sqrt(var(double(PH3(:)))) / mean(double(PH3(:))) ;
                var_temp(iz) = var(double(PH3(:))) ;
                int_temp(iz) = sum(double(PH3(:))) ;
            end
            if wholexy
                var_lst_cv = cv_temp ;
                var_lst_var = var_temp ;
                var_lst_int = int_temp ;
            end
            [~,phase3_bestzPRE(n_frame)] = max(var_temp);
        end

        figure(ip.exp.frm2spl(n_frame) + (c-1)*max(ip.exp.frm2spl));
        subplot(1,3,1);
        hold on;
        %phase3 defines what colors to use
        if c == 1
            color_used{n_frame} = which_color(gca);
        end

            plot(1: spl2znum(current_spl), var_lst_cv/max(var_lst_cv), 'DisplayName', ...
                ['XY position ', num2str(ip.exp.frm2img(n_frame))], 'Color', color_used{n_frame});
            xlabel('z slice')
            ylabel([channels{c} ' pixel CV / max ' channels{c} ' CV']);
            title(['Normalized ' channels{c} ' pixel CV vs. z-slice'])

            subplot(1,3,2);
            hold on;
            %color_used{n_frame} = which_color(gca);
            plot(1: spl2znum(current_spl), var_lst_var/max(var_lst_var), 'DisplayName', ...
                ['XY position ', num2str(ip.exp.frm2img(n_frame))], 'Color', color_used{n_frame});
            xlabel('z slice')
            ylabel([channels{c} ' pixel variance / max ' channels{c} ' variance']);
            title(['Normalized ' channels{c} ' pixel var. vs. z-slice'])

            subplot(1,3,3);
            hold on;
            %color_used{n_frame} = which_color(gca);
            plot(1: spl2znum(current_spl), var_lst_int/max(var_lst_int), 'DisplayName', ...
                ['XY position ', num2str(ip.exp.frm2img(n_frame))], 'Color', color_used{n_frame});
            xlabel('z slice')
            ylabel([channels{c} ' pixel intensity / max ' channels{c} ' intensity']);
            title(['Normalized ' channels{c} ' pixel int. vs. z-slice'])
    end
    
    %%% Saving figures of current channel
    
    range_to_process = unique(ip.exp.frm2spl(1:ip.exp.totalframes));
    for i = transpose(range_to_process)
        figure(i + (c-1)*max(ip.exp.frm2spl));
        set(gcf,'units','points','position',[-1000,250,720,326]);
        if strcmp(channels{c}, 'phase3') && wholexy
            sgtitle([date ' Sample ', num2str(i), ': ',ip.sample{1,i}.name, ' [Analysis over whole XY]']);
        else
            sgtitle([date ' Sample ', num2str(i), ': ',ip.sample{1,i}.name, ' [Analysis over each cell]']);
        end
        for ii = 1:3
            subplot(1,3,ii);
            legend('location', 'best');
            ylim([y_min(ii) 1]);
            
            grid on;
        end

        saveas(figure(i + (c-1)*max(ip.exp.frm2spl)), ...
            ['postsegmentation_fluor_analysis\' channels{c} '_CV_var_int_vs_zslice\Sample_', ...
            num2str(i, '%02d'), '__Norm' channels{c} 'Pixel__CV_var_int_vs_z.png']);
        saveas(figure(i + (c-1)*max(ip.exp.frm2spl)), ...
            ['postsegmentation_fluor_analysis\' channels{c} '_CV_var_int_vs_zslice\Sample_', ...
            num2str(i, '%02d'), '__Norm' channels{c} 'Pixel__CV_var_int_vs_z.fig']);
        
        close(figure(i + (c-1)*max(ip.exp.frm2spl)));
    end
end

%%% Best phase3 variance vs. best FLUOR variance
cmap = hsv(max(ip.exp.frm2spl));
phase3_bestzPRE = phase3_bestz;

f = figure ();
subplot(2,2,1);
for i = 1:max(ip.exp.frm2spl)
    xy = [phase3_bestzPRE(nonzeros(ip.exp.splimg2frm(i,:))) - max(ip.exp.frm2spl)/2 * .05 + i * .05,...
        fluor_bestvarz(nonzeros(ip.exp.splimg2frm(i,:))) - max(ip.exp.frm2spl)/2 * .05 + i * .05];
    [C,ia,ic] = unique(xy, 'rows');
    a_counts = accumarray(ic,1);

    scatter(C(:,1), C(:,2), a_counts*5, 'markeredgecolor', 'none',...
        'markerfacecolor', cmap(i,:), 'DisplayName', ['Sample ', num2str(i)]);

%     scatter(unique(phase3_bestz(nonzeros(ip.exp.splimg2frm(i,:))) - max(ip.exp.frm2spl)/2 * .05 + i * .05), ...
%         unique(fluor_bestvarz(nonzeros(ip.exp.splimg2frm(i,:))) - max(ip.exp.frm2spl)/2 * .05 + i * .05), ...
%         a_counts/max(a_counts) * 15, 'markerfacecolor', cmap(i,:), 'DisplayName', ['Sample ', num2str(i)]);
    hold on;
end
plot(0:8, 0:8, 'Color', 'black', 'DisplayName', 'Slope = 1');
xlim([0 8]);
ylim([0 8]);
title(['Z w/ max Ph3 var. pre-seg. vs. Z w/ max ' channels{fluor_chan} ' var. post-seg.']);
xlabel('Phase 3 highest variance z slice');
ylabel([channels(fluor_chan) ' highest variance z slice']);
grid on;
legend('location', 'best');

subplot(2,2,3);
xy = [phase3_bestzPRE(nonzeros(ip.exp.splimg2frm(:))), fluor_bestvarz(nonzeros(ip.exp.splimg2frm(:)))];
[C,ia,ic] = unique(xy, 'rows');
a_counts = accumarray(ic,1);
scatter(C(:,1), C(:,2), a_counts*5, 'markeredgecolor', 'none',...
    'markerfacecolor', 'black', 'DisplayName', ['All samples together']);
hold on;
plot(0:8, 0:8, 'Color', 'black', 'DisplayName', 'Slope = 1');
xlim([0 8]);
ylim([0 8]);
title(['Z w/ max Ph3 var. pre-seg. vs. Z w/ max ' channels{fluor_chan} ' var. post-seg']);
xlabel('Phase 3 highest variance z slice');
ylabel([channels(fluor_chan) ' highest variance z slice']);
grid on;
legend('location', 'best');


subplot(2,2,2);
for i = 1:max(ip.exp.frm2spl)
    xy = [phase3_bestzPRE(nonzeros(ip.exp.splimg2frm(i,:))) - max(ip.exp.frm2spl)/2 * .05 + i * .05,...
        fluor_bestCVz(nonzeros(ip.exp.splimg2frm(i,:))) - max(ip.exp.frm2spl)/2 * .05 + i * .05];
    [C,ia,ic] = unique(xy, 'rows');
    a_counts = accumarray(ic,1);

    scatter(C(:,1), C(:,2), a_counts/max(a_counts) * 50, 'markeredgecolor', 'none',...
        'markerfacecolor', cmap(i,:), 'DisplayName', ['Sample ', num2str(i)]);
    
%     scatter(phase3_bestz(nonzeros(ip.exp.splimg2frm(i,:))) - max(ip.exp.frm2spl)/2 * .05 + i * .05, ...
%         fluor_bestCVz(nonzeros(ip.exp.splimg2frm(i,:))) - max(ip.exp.frm2spl)/2 * .05 + i * .05, ...
%         a_counts/max(a_counts) * 15, 'markerfacecolor', cmap(i,:), 'DisplayName', ['Sample ', num2str(i)]);
    hold on;
end
plot(0:8, 0:8, 'Color', 'black', 'DisplayName', 'Slope = 1');
xlim([0 8]);
ylim([0 8]);
title(['Z w/ max Ph3 var pre-seg. vs. Z w/ max ' channels{fluor_chan} ' CV post-seg.']);
xlabel('Phase 3 highest CV z slice');
ylabel([channels{fluor_chan} ' highest CV z slice']);
grid on;
legend('location', 'best');

subplot(2,2,4);
xy = [phase3_bestzPRE(nonzeros(ip.exp.splimg2frm(:))), fluor_bestCVz(nonzeros(ip.exp.splimg2frm(:)))];
[C,ia,ic] = unique(xy, 'rows');
a_counts = accumarray(ic,1);
scatter(C(:,1), C(:,2), a_counts*5, 'markeredgecolor', 'none',...
    'markerfacecolor', 'black', 'DisplayName', ['All samples together']);
hold on;
plot(0:8, 0:8, 'Color', 'black', 'DisplayName', 'Slope = 1');
xlim([0 8]);
ylim([0 8]);
title(['Z w/ max Ph3 var pre-seg vs. Z w/ max ' channels{fluor_chan} ' CV post-seg.']);
xlabel('Phase 3 highest CV z slice');
ylabel([channels{fluor_chan} ' highest CV z slice']);
grid on;
legend('location', 'best');

set(gcf,'units','points','position',[-1200,50,720*1.5,275*1.25*2]);

saveas(f, ['postsegmentation_fluor_analysis\phase3_maxvar_z_vs_' ...
    channels{fluor_chan} '_maxvar_z_maxcv_z.png']);
saveas(f, ['postsegmentation_fluor_analysis\phase3_maxvar_z_vs_' ...
    channels{fluor_chan} '_maxvar_z_maxcv_z.fig']);


%%% Best phase3 variance vs. best FLUOR variance bar of difference
%{
pick_z_fluorvar = transpose(fluor_bestvarz(ip.exp.frm2spl ~= 8));
[all_x, all_y] = ind2sub(size(fluor_varz(:,:,1)), 1:numel(fluor_varz(:,:,1)));
all_x = all_x(fluor_varz(:,:,1) ~= 0);
all_y = all_y(fluor_varz(:,:,1) ~= 0);
best_var = fluor_varz(sub2ind(size(fluor_varz), all_x, all_y, pick_z_fluorvar));

best_var_t = best_var(pick_z_fluorvar > 2 & pick_z_fluorvar < ip.image.zrange-1);
pick_z_fluorvar_t = pick_z_fluorvar(pick_z_fluorvar > 2 & pick_z_fluorvar < ip.image.zrange-1);
all_d1_var = [abs(best_var - fluor_varz(sub2ind(size(fluor_varz), all_x, all_y, pick_z_fluorvar+1))), ...
    abs(best_var - fluor_varz(sub2ind(size(fluor_varz), all_x, all_y, pick_z_fluorvar-1)))];
all_d2_var = [abs(best_var_t - fluor_varz(sub2ind(size(fluor_varz), ...
    all_x(pick_z_fluorvar > 2 & pick_z_fluorvar < ip.image.zrange-1), ...
    all_y(pick_z_fluorvar > 2 & pick_z_fluorvar < ip.image.zrange-1), pick_z_fluorvar_t+2))), ...
    abs(best_var_t - fluor_varz(sub2ind(size(fluor_varz), ...
    all_x(pick_z_fluorvar > 2 & pick_z_fluorvar < ip.image.zrange-1), ...
    all_y(pick_z_fluorvar > 2 & pick_z_fluorvar < ip.image.zrange-1), pick_z_fluorvar_t-2)))];

pick_z_ph3 = transpose(phase3_bestzPRE(ip.exp.frm2spl ~= 8));
ph3_var = fluor_varz(sub2ind(size(fluor_varz), all_x, all_y, pick_z_ph3));
d_var = abs(best_var - ph3_var);
idx_var = abs(pick_z_fluorvar - pick_z_ph3);
%
pick_z_fluorcv = transpose(fluor_bestCVz(ip.exp.frm2spl ~= 8));
[all_x, all_y] = ind2sub(size(fluor_CVz(:,:,1)), 1:numel(fluor_CVz(:,:,1)));
all_x = all_x(fluor_CVz(:,:,1) ~= 0);
all_y = all_y(fluor_CVz(:,:,1) ~=0);
best_CV = fluor_CVz(sub2ind(size(fluor_CVz), all_x, all_y, pick_z_fluorcv));

best_CV_t = best_CV(pick_z_fluorcv > 2 & pick_z_fluorcv < ip.image.zrange-1);
pick_z_fluorcv_t = pick_z_fluorcv(pick_z_fluorcv > 2 & pick_z_fluorcv < ip.image.zrange-1);
all_d1_CV = [abs(best_CV(pick_z_fluorcv > 1 & pick_z_fluorcv < ip.image.zrange) - fluor_CVz(sub2ind(size(fluor_CVz), ...
    all_x(pick_z_fluorcv > 1 & pick_z_fluorcv < ip.image.zrange), ...
    all_y(pick_z_fluorcv > 1 & pick_z_fluorcv < ip.image.zrange), ...
    pick_z_fluorcv(pick_z_fluorcv > 1 & pick_z_fluorcv < ip.image.zrange)+1))), ...
    abs(best_CV(pick_z_fluorcv > 1 & pick_z_fluorcv < ip.image.zrange) - fluor_CVz(sub2ind(size(fluor_CVz), ...
    all_x(pick_z_fluorcv > 1 & pick_z_fluorcv < ip.image.zrange), ...
    all_y(pick_z_fluorcv > 1 & pick_z_fluorcv < ip.image.zrange), ...
    pick_z_fluorcv(pick_z_fluorcv > 1 & pick_z_fluorcv < ip.image.zrange)-1)))];
all_d2_CV = [abs(best_CV_t - fluor_CVz(sub2ind(size(fluor_CVz), ...
    all_x(pick_z_fluorcv > 2 & pick_z_fluorcv < ip.image.zrange-1), ...
    all_y(pick_z_fluorcv > 2 & pick_z_fluorcv < ip.image.zrange-1), pick_z_fluorcv_t+2))), ...
    abs(best_CV_t - fluor_CVz(sub2ind(size(fluor_CVz), ...
    all_x(pick_z_fluorcv > 2 & pick_z_fluorcv < ip.image.zrange-1), ...
    all_y(pick_z_fluorcv > 2 & pick_z_fluorcv < ip.image.zrange-1), pick_z_fluorcv_t-2)))];

ph3_CV = fluor_CVz(sub2ind(size(fluor_CVz), all_x, all_y, pick_z_ph3));
d_CV = abs(best_CV - ph3_CV);
idx_cv = abs(pick_z_fluorcv - pick_z_ph3);
%
f = figure ();
subplot(1,2,1);
bar(1, mean(all_d1_var), 'DisplayName', 'Mean \Deltavar. 1 zstep from max', 'facecolor', cmap(1,:));
hold on;
errorbar(1, mean(all_d1_var), std(all_d1_var), 'DisplayName', 'Std of all (+/- 1 from best Z)', 'color', cmap(1,:));
hold on;
swarmchart(1*ones(size(all_d1_var)), all_d1_var, 'HandleVisibility', 'off', 'MarkerFaceColor', cmap(1,:),'MarkerEdgeAlpha',0.4);
hold on;
%
bar(2, mean(all_d2_var), 'DisplayName', 'Mean \Deltavar. 2 zstep from max', 'facecolor', cmap(2,:));
hold on;
errorbar(2, mean(all_d2_var), std(all_d2_var), 'DisplayName', 'Std of all (+/- 2 from best Z)', 'color', cmap(2,:));
hold on;
swarmchart(2*ones(size(all_d2_var)), all_d2_var, 'HandleVisibility', 'off', 'MarkerFaceColor', cmap(2,:),'MarkerEdgeAlpha',0.4);
hold on;
%
bar(3, mean(d_var(idx_var == 1)), 'DisplayName', 'Mean \Deltavar. in samp. w/ 1 zstep diff.', 'facecolor', cmap(3,:));
hold on;
errorbar(3, mean(d_var(idx_var == 1)), std(d_var(idx_var == 1)), 'DisplayName', 'Std of all +/- 1 z diff', 'color', cmap(3,:));
hold on;
swarmchart(3*ones(size(d_var(idx_var == 1))), d_var(idx_var == 1), 'HandleVisibility', 'off', 'MarkerFaceColor', cmap(3,:),'MarkerEdgeAlpha',0.4);
hold on;
%
bar(4, mean(d_var(idx_var == 2)), 'DisplayName', 'Mean \Deltavar. in samp. w/ 2 zstep diff', 'facecolor', cmap(4,:));
hold on;
errorbar(4, mean(d_var(idx_var == 2)), std(d_var(idx_var == 2)), 'DisplayName', 'Std of all +/- 2 z diff', 'color', cmap(4,:));
hold on;
swarmchart(4*ones(size(d_var(idx_var == 2))), d_var(idx_var == 2), 'HandleVisibility', 'off', 'MarkerFaceColor', cmap(4,:),'MarkerEdgeAlpha',0.4);
hold on;

xlim([.5 4.5]);
ylim([0 .6]);
title(['Z w/ max Ph3 var. pre-seg. vs. Z w/ max ' channels{fluor_chan} ' var. post-seg.']);
xticks([1,2,3,4]);
xticklabels({'Fluor focus \newline 1 slice from max', 'Fluor focus \newline 2slices from max', ...
    'Fluor focus \newline 1 slice from Ph3', 'Fluor focus \newline 2 slices from Ph3'})
ylabel([channels(fluor_chan) ' difference in normalized variance']);
grid on;
legend('location', 'best');


subplot(1,2,2);
bar(1, mean(all_d1_CV), 'DisplayName', 'Mean \Deltacv 1 zstep from max', 'facecolor', cmap(1,:));
hold on;
errorbar(1, mean(all_d1_CV), std(all_d1_CV), 'DisplayName', 'Std of all (+/- 1 from best Z)', 'color', cmap(1,:));
hold on;
swarmchart(1*ones(size(all_d1_CV)), all_d1_CV, 'HandleVisibility', 'off', 'MarkerFaceColor', cmap(1,:),'MarkerEdgeAlpha',0.4);
hold on;
%
bar(2, mean(all_d2_CV), 'DisplayName', 'Mean \Deltacv 2 zstep from max', 'facecolor', cmap(2,:));
hold on;
errorbar(2, mean(all_d2_CV), std(all_d2_CV), 'DisplayName', 'Std of all (+/- 2 from best Z)', 'color', cmap(2,:));
hold on;
swarmchart(2*ones(size(all_d2_CV)), all_d2_CV, 'HandleVisibility', 'off', 'MarkerFaceColor', cmap(2,:),'MarkerEdgeAlpha',0.4);
hold on;
%
bar(3, mean(d_CV(idx_cv == 1)), 'DisplayName', 'Mean \Deltacv in samp. w/ 1 zstep diff.', 'facecolor', cmap(3,:));
hold on;
errorbar(3, mean(d_CV(idx_cv == 1)), std(d_CV(idx_cv == 1)), 'DisplayName', 'Std of all +/- 1 z diff', 'color', cmap(3,:));
hold on;
swarmchart(3*ones(size(d_CV(idx_cv == 1))), d_CV(idx_cv == 1), 'HandleVisibility', 'off', 'MarkerFaceColor', cmap(3,:),'MarkerEdgeAlpha',0.4);
hold on;
%
bar(4, mean(d_CV(idx_cv == 2)), 'DisplayName', 'Mean \Deltacv in samp. w/ 2 zstep diff', 'facecolor', cmap(4,:));
hold on;
errorbar(4, mean(d_CV(idx_cv == 2)), std(d_CV(idx_cv == 2)), 'DisplayName', 'Std of all +/- 2 z diff', 'color', cmap(4,:));
hold on;
swarmchart(4*ones(size(d_CV(idx_cv == 2))), d_CV(idx_cv == 2), 'HandleVisibility', 'off', 'MarkerFaceColor', cmap(4,:),'MarkerEdgeAlpha',0.4);
hold on;

xlim([.5 4.5]);
ylim([0 .6]);
title(['Z w/ max Ph3 CV pre-seg. vs. Z w/ max ' channels{fluor_chan} ' CV post-seg.']);
xticks([1,2,3,4]);
xticklabels({'Fluor focus \newline 1 slice from max', 'Fluor focus \newline 2slices from max', ...
    'Fluor focus \newline 1 slice from Ph3', 'Fluor focus \newline 2 slices from Ph3'})
ylabel([channels(fluor_chan) 'Difference in normalized CV ']);
grid on;
legend('location', 'best');


set(gcf,'units','points','position',[-1200,50,720*1.5,275*1.25*2]);

saveas(f, ['postsegmentation_fluor_analysis\phase3_maxvar_z_vs_' ...
    channels{fluor_chan} '_maxvar_z_maxcv_z_BAR.png']);
saveas(f, ['postsegmentation_fluor_analysis\phase3_maxvar_z_vs_' ...
    channels{fluor_chan} '_maxvar_z_maxcv_z_BAR.fig']);


% function padded_vector = pad_vector(input_vector, model_matrix)
%     padding_per_row = sum(model_matrix == 0, 2);
%     [rows,cols] = size(model_matrix);
%     padded_vector = [];
%     for r = 1:rows
%         padded_vector = [padded_vector, input_vector((r-1)*cols + 1:(cols-padding_per_row)), zeros(1,padding_per_row)];
%     end
% end
%}


%% Variance & Intensity & CV vs. z-slice POST (cell by cell comparison)
[num_sam, num_img] = size(ip.exp.splimg2frm);

var_lst_cv = cell(length(channels),1);
var_lst_var = cell(length(channels),1);
phase3_bestzPRE = [];

keeptrackofthesample = [];
for c = 1:length(channels)
    %z slice #
    spl2znum = ip.image.zrange .* ...
        ones(1,max(ip.exp.sampleList)) ;
    
    %work our way down for y axis limits
    y_min(1) = 1;
    y_min(2) = 1;
    y_min(3) = 1;
    
    var_lst_cv{c} = [];
    var_lst_var{c} = [];
    %loop frames
    for n_frame = 1:ip.exp.totalframes
        n_image = ip.exp.frm2img(n_frame);
        %progress updates
        if n_image == 1
            fprintf(1, newline);
            fprintf(1,['Sample ' num2str(ip.exp.frm2spl(n_frame)) ' of ' ...
                num2str(max(ip.exp.sampleList)) '.' newline]);
            progress_2 = [];
        end
        for d_=1:1:size(progress_2,2) ; fprintf(1,'\b') ; end
        progress_2 = [sprintf('\t') 'Image ' num2str(n_image) ' in ' ...
            num2str(min(ip.sample{ip.exp.frm2spl(n_frame)}.idx)) '-' ...
            num2str(max(ip.sample{ip.exp.frm2spl(n_frame)}.idx)) '.'] ;
        fprintf(1,progress_2) ;
        %load image
        load([ip.seg.dir ip.seg.base_name 'seg' num2str(n_frame,'%03d') '.mat'],'LcFull') ;
        imdigitnum_ = 2 ;
        current_spl = ip.exp.frm2spl(n_frame) ;
        if max(ip.sample{current_spl}.idx) <= 9
            imdigitnum_ = 1 ;
        end
        var_lst = zeros(1,ip.image.zrange) ;

        CellMap = LcFull;
        nnn=double(LcFull);
        N=max(nnn(:));
        
        cv_temp = zeros(spl2znum(current_spl),N);
        var_temp = zeros(spl2znum(current_spl),N);
        %loop z
        for iz = 1: spl2znum(current_spl)
        PH3 = imread([ip.exp.path 'images\' ...
                ip.image.base_name '_' ...
                num2str(ip.exp.frm2spl(n_frame),'%03d') ...
                'xy' num2str(ip.exp.frm2img(n_frame),['%0' num2str(imdigitnum_ ,'%01d') 'd']) ...
                ['z' num2str(iz ,'%01d') 'c' num2str(c) '.tif'] ]);

            %loop cells
            for index=1:N
                Cell = immultiply( CellMap == index, PH3);
                int_per_pix = Cell(CellMap == index);

                cv_temp(iz,index) = sqrt(var(double(int_per_pix(:)))) / mean(double(int_per_pix(:))) ;
                var_temp(iz,index) = var(double(int_per_pix(:))) ;
            end
        end
        
        % max for each cell
        [~, cv_temp] = max(cv_temp) ;
        var_lst_cv{c} = [var_lst_cv{c}, cv_temp];
        [~, var_temp] = max(var_temp) ;
        var_lst_var{c} = [var_lst_var{c}, var_temp];
        
%         cv_temp = zeros(spl2znum(current_spl),N);
%         var_temp = zeros(spl2znum(current_spl),N);
%         %go over the whole image only (rather than cell by cell)
%         if strcmp(channels{c}, 'phase3')
%             var_temp = zeros(spl2znum(current_spl),1);
%             for iz = 1: spl2znum(current_spl)
%                 PH3 = imread([ip.exp.path 'images\' ...
%                     ip.image.base_name '_' ...
%                     num2str(ip.exp.frm2spl(n_frame),'%03d') ...
%                     'xy' num2str(ip.exp.frm2img(n_frame),['%0' num2str(imdigitnum_ ,'%01d') 'd']) ...
%                     ['z' num2str(iz ,'%01d') 'c' num2str(c) '.tif'] ]);
%                 
% %                 cv_temp(iz) = sqrt(var(double(PH3(:)))) / mean(double(PH3(:))) ;
%                 var_temp(iz) = var(double(PH3(:))) ;
%             end
%             [~,var_temp] = max(var_temp);
%             phase3_bestzPRE = [phase3_bestzPRE, var_temp];
%         end
%         
        keeptrackofthesample = [keeptrackofthesample; ip.exp.frm2spl(n_frame)*ones(length(var_temp),1)];
    end
end
%
%plotting
dot_scale = 1;
for fluor_chan = 2:length(channels)
    %%% Best phase3 variance vs. best FLUOR variance
    cmap = hsv(max(ip.exp.frm2spl));
    f = figure ();
    subplot(2,2,1);
    for i = 1:max(ip.exp.frm2spl)
        x = var_lst_var{1}(keeptrackofthesample(1:length(keeptrackofthesample)/length(channels)) == i);
        y = var_lst_var{fluor_chan}(keeptrackofthesample(1:length(keeptrackofthesample)/length(channels)) == i);
        xy = [x' - max(ip.exp.frm2spl)/2 * .1 + i * .1,...
            y' - max(ip.exp.frm2spl)/2 * .1 + i * .1];
        [C,~,ic] = unique(xy, 'rows');
        a_counts = accumarray(ic,1);

        scatter(C(:,1), C(:,2), a_counts*dot_scale, 'markeredgecolor', 'none',...
            'markerfacecolor', cmap(i,:), 'DisplayName', ['Sample ', num2str(i)]);
        hold on;
    end
    plot(0:8, 0:8, 'Color', 'black', 'DisplayName', 'Slope = 1');
    xlim([0 8]);
    ylim([0 8]);
    title(['Z w/ max Ph3 var. cell vs. Z w/ max ' channels{fluor_chan} ' var. cell']);
    xlabel('Phase 3 highest variance z slice');
    ylabel([channels(fluor_chan) ' highest variance z slice']);
    grid on;
    legend('location', 'best');

    subplot(2,2,3);
    xy = [var_lst_var{1}(:), var_lst_var{fluor_chan}(:)];
    [C,~,ic] = unique(xy, 'rows');
    a_counts = accumarray(ic,1);
    scatter(C(:,1), C(:,2), a_counts*dot_scale, 'markeredgecolor', 'none',...
        'markerfacecolor', 'black', 'DisplayName', ['All samples together']);
    hold on;
    plot(0:8, 0:8, 'Color', 'black', 'DisplayName', 'Slope = 1');
    hold on;
    x = 1:ip.image.zrange;
    y = zeros(1,ip.image.zrange);
    for i = x
        y(i) = mean(var_lst_var{fluor_chan}(var_lst_var{1} == i));
        y_er(i) = std(var_lst_var{fluor_chan}(var_lst_var{1} == i));
    end
    errorbar(x,y,y_er, 'displayname', 'Mean fluor z slice +/- std'); hold on;
    xlim([0 8]);
    ylim([0 8]);
    title(['Z w/ max Ph3 var. cell vs. Z w/ max ' channels{fluor_chan} ' var. cell']);
    xlabel('Phase 3 highest variance z slice');
    ylabel([channels(fluor_chan) ' highest variance z slice']);
    grid on;
    legend('location', 'best');


    subplot(2,2,2);
    for i = 1:max(ip.exp.frm2spl)
        x = var_lst_cv{1}(keeptrackofthesample(1:length(keeptrackofthesample)/length(channels)) == i);
        y = var_lst_cv{fluor_chan}(keeptrackofthesample(1:length(keeptrackofthesample)/length(channels)) == i);
        xy = [x' - max(ip.exp.frm2spl)/2 * .1 + i * .1,...
            y' - max(ip.exp.frm2spl)/2 * .1 + i * .1];
        [C,~,ic] = unique(xy, 'rows');
        a_counts = accumarray(ic,1);

        scatter(C(:,1), C(:,2), a_counts*dot_scale, 'markeredgecolor', 'none',...
            'markerfacecolor', cmap(i,:), 'DisplayName', ['Sample ', num2str(i)]);
        hold on;
    end
    plot(0:8, 0:8, 'Color', 'black', 'DisplayName', 'Slope = 1');
    xlim([0 8]);
    ylim([0 8]);
    title(['Z w/ max Ph3 var cell vs. Z w/ max ' channels{fluor_chan} ' CV cell']);
    xlabel('Phase 3 highest CV z slice');
    ylabel([channels{fluor_chan} ' highest CV z slice']);
    grid on;
    legend('location', 'best');

    subplot(2,2,4);
    xy = [var_lst_var{1}(:), var_lst_cv{fluor_chan}(:)];
    [C,~,ic] = unique(xy, 'rows');
    a_counts = accumarray(ic,1);
    scatter(C(:,1), C(:,2), a_counts*dot_scale, 'markeredgecolor', 'none',...
        'markerfacecolor', 'black', 'DisplayName', ['All samples together']);
    hold on;
    plot(0:8, 0:8, 'Color', 'black', 'DisplayName', 'Slope = 1'); 
    hold on;
    x = 1:ip.image.zrange;
    y = zeros(1,ip.image.zrange);
    for i = x
        y(i) = mean(var_lst_cv{fluor_chan}(var_lst_cv{1} == i));
        y_er(i) = std(var_lst_cv{fluor_chan}(var_lst_cv{1} == i));
    end
    errorbar(x,y,y_er, 'displayname', 'Mean fluor z slice +/- std'); hold on;
    xlim([0 8]);
    ylim([0 8]);
    title(['Z w/ max Ph3 var cell vs. Z w/ max ' channels{fluor_chan} ' CV cell']);
    xlabel('Phase 3 highest CV z slice');
    ylabel([channels{fluor_chan} ' highest CV z slice']);
    grid on;
    legend('location', 'best');

    set(gcf,'units','points','position',[-1200,50,720*1.5,275*1.25*2]);

    saveas(f, ['postsegmentation_fluor_analysis\phase3_maxvar_z_vs_' ...
        channels{fluor_chan} '_maxvar_z_maxcv_z__cell_by_cell.png']);
    saveas(f, ['postsegmentation_fluor_analysis\phase3_maxvar_z_vs_' ...
        channels{fluor_chan} '_maxvar_z_maxcv_z__cell_by_cell.fig']);


%     %%% Best phase3 variance vs. best FLUOR variance bar of difference
%     pick_z_fluorvar = transpose(fluor_bestvarz(ip.exp.frm2spl ~= 8));
%     [all_x, all_y] = ind2sub(size(fluor_varz(:,:,1)), 1:numel(fluor_varz(:,:,1)));
%     all_x = all_x(fluor_varz(:,:,1) ~= 0);
%     all_y = all_y(fluor_varz(:,:,1) ~= 0);
%     best_var = fluor_varz(sub2ind(size(fluor_varz), all_x, all_y, pick_z_fluorvar));
% 
%     best_var_t = best_var(pick_z_fluorvar > 2 & pick_z_fluorvar < ip.image.zrange-1);
%     pick_z_fluorvar_t = pick_z_fluorvar(pick_z_fluorvar > 2 & pick_z_fluorvar < ip.image.zrange-1);
%     all_d1_var = [abs(best_var - fluor_varz(sub2ind(size(fluor_varz), all_x, all_y, pick_z_fluorvar+1))), ...
%         abs(best_var - fluor_varz(sub2ind(size(fluor_varz), all_x, all_y, pick_z_fluorvar-1)))];
%     all_d2_var = [abs(best_var_t - fluor_varz(sub2ind(size(fluor_varz), ...
%         all_x(pick_z_fluorvar > 2 & pick_z_fluorvar < ip.image.zrange-1), ...
%         all_y(pick_z_fluorvar > 2 & pick_z_fluorvar < ip.image.zrange-1), pick_z_fluorvar_t+2))), ...
%         abs(best_var_t - fluor_varz(sub2ind(size(fluor_varz), ...
%         all_x(pick_z_fluorvar > 2 & pick_z_fluorvar < ip.image.zrange-1), ...
%         all_y(pick_z_fluorvar > 2 & pick_z_fluorvar < ip.image.zrange-1), pick_z_fluorvar_t-2)))];
% 
%     pick_z_ph3 = transpose(phase3_bestzPRE(ip.exp.frm2spl ~= 8));
%     ph3_var = fluor_varz(sub2ind(size(fluor_varz), all_x, all_y, pick_z_ph3));
%     d_var = abs(best_var - ph3_var);
%     idx_var = abs(pick_z_fluorvar - pick_z_ph3);
%     %
%     pick_z_fluorcv = transpose(fluor_bestCVz(ip.exp.frm2spl ~= 8));
%     [all_x, all_y] = ind2sub(size(fluor_CVz(:,:,1)), 1:numel(fluor_CVz(:,:,1)));
%     all_x = all_x(fluor_CVz(:,:,1) ~= 0);
%     all_y = all_y(fluor_CVz(:,:,1) ~=0);
%     best_CV = fluor_CVz(sub2ind(size(fluor_CVz), all_x, all_y, pick_z_fluorcv));
% 
%     best_CV_t = best_CV(pick_z_fluorcv > 2 & pick_z_fluorcv < ip.image.zrange-1);
%     pick_z_fluorcv_t = pick_z_fluorcv(pick_z_fluorcv > 2 & pick_z_fluorcv < ip.image.zrange-1);
%     all_d1_CV = [abs(best_CV(pick_z_fluorcv > 1 & pick_z_fluorcv < ip.image.zrange) - fluor_CVz(sub2ind(size(fluor_CVz), ...
%         all_x(pick_z_fluorcv > 1 & pick_z_fluorcv < ip.image.zrange), ...
%         all_y(pick_z_fluorcv > 1 & pick_z_fluorcv < ip.image.zrange), ...
%         pick_z_fluorcv(pick_z_fluorcv > 1 & pick_z_fluorcv < ip.image.zrange)+1))), ...
%         abs(best_CV(pick_z_fluorcv > 1 & pick_z_fluorcv < ip.image.zrange) - fluor_CVz(sub2ind(size(fluor_CVz), ...
%         all_x(pick_z_fluorcv > 1 & pick_z_fluorcv < ip.image.zrange), ...
%         all_y(pick_z_fluorcv > 1 & pick_z_fluorcv < ip.image.zrange), ...
%         pick_z_fluorcv(pick_z_fluorcv > 1 & pick_z_fluorcv < ip.image.zrange)-1)))];
%     all_d2_CV = [abs(best_CV_t - fluor_CVz(sub2ind(size(fluor_CVz), ...
%         all_x(pick_z_fluorcv > 2 & pick_z_fluorcv < ip.image.zrange-1), ...
%         all_y(pick_z_fluorcv > 2 & pick_z_fluorcv < ip.image.zrange-1), pick_z_fluorcv_t+2))), ...
%         abs(best_CV_t - fluor_CVz(sub2ind(size(fluor_CVz), ...
%         all_x(pick_z_fluorcv > 2 & pick_z_fluorcv < ip.image.zrange-1), ...
%         all_y(pick_z_fluorcv > 2 & pick_z_fluorcv < ip.image.zrange-1), pick_z_fluorcv_t-2)))];
% 
%     ph3_CV = fluor_CVz(sub2ind(size(fluor_CVz), all_x, all_y, pick_z_ph3));
%     d_CV = abs(best_CV - ph3_CV);
%     idx_cv = abs(pick_z_fluorcv - pick_z_ph3);
%     %
%     f = figure ();
%     subplot(1,2,1);
%     bar(1, mean(all_d1_var), 'DisplayName', 'Mean \Deltavar. 1 zstep from max', 'facecolor', cmap(1,:));
%     hold on;
%     errorbar(1, mean(all_d1_var), std(all_d1_var), 'DisplayName', 'Std of all (+/- 1 from best Z)', 'color', cmap(1,:));
%     hold on;
%     swarmchart(1*ones(size(all_d1_var)), all_d1_var, 'HandleVisibility', 'off', 'MarkerFaceColor', cmap(1,:),'MarkerEdgeAlpha',0.4);
%     hold on;
%     %
%     bar(2, mean(all_d2_var), 'DisplayName', 'Mean \Deltavar. 2 zstep from max', 'facecolor', cmap(2,:));
%     hold on;
%     errorbar(2, mean(all_d2_var), std(all_d2_var), 'DisplayName', 'Std of all (+/- 2 from best Z)', 'color', cmap(2,:));
%     hold on;
%     swarmchart(2*ones(size(all_d2_var)), all_d2_var, 'HandleVisibility', 'off', 'MarkerFaceColor', cmap(2,:),'MarkerEdgeAlpha',0.4);
%     hold on;
%     %
%     bar(3, mean(d_var(idx_var == 1)), 'DisplayName', 'Mean \Deltavar. in samp. w/ 1 zstep diff.', 'facecolor', cmap(3,:));
%     hold on;
%     errorbar(3, mean(d_var(idx_var == 1)), std(d_var(idx_var == 1)), 'DisplayName', 'Std of all +/- 1 z diff', 'color', cmap(3,:));
%     hold on;
%     swarmchart(3*ones(size(d_var(idx_var == 1))), d_var(idx_var == 1), 'HandleVisibility', 'off', 'MarkerFaceColor', cmap(3,:),'MarkerEdgeAlpha',0.4);
%     hold on;
%     %
%     bar(4, mean(d_var(idx_var == 2)), 'DisplayName', 'Mean \Deltavar. in samp. w/ 2 zstep diff', 'facecolor', cmap(4,:));
%     hold on;
%     errorbar(4, mean(d_var(idx_var == 2)), std(d_var(idx_var == 2)), 'DisplayName', 'Std of all +/- 2 z diff', 'color', cmap(4,:));
%     hold on;
%     swarmchart(4*ones(size(d_var(idx_var == 2))), d_var(idx_var == 2), 'HandleVisibility', 'off', 'MarkerFaceColor', cmap(4,:),'MarkerEdgeAlpha',0.4);
%     hold on;
% 
%     xlim([.5 4.5]);
%     ylim([0 .6]);
%     title(['Z w/ max Ph3 var. pre-seg. vs. Z w/ max ' channels{fluor_chan} ' var. post-seg.']);
%     xticks([1,2,3,4]);
%     xticklabels({'Fluor focus \newline 1 slice from max', 'Fluor focus \newline 2slices from max', ...
%         'Fluor focus \newline 1 slice from Ph3', 'Fluor focus \newline 2 slices from Ph3'})
%     ylabel([channels(fluor_chan) ' difference in normalized variance']);
%     grid on;
%     legend('location', 'best');
% 
% 
%     subplot(1,2,2);
%     bar(1, mean(all_d1_CV), 'DisplayName', 'Mean \Deltacv 1 zstep from max', 'facecolor', cmap(1,:));
%     hold on;
%     errorbar(1, mean(all_d1_CV), std(all_d1_CV), 'DisplayName', 'Std of all (+/- 1 from best Z)', 'color', cmap(1,:));
%     hold on;
%     swarmchart(1*ones(size(all_d1_CV)), all_d1_CV, 'HandleVisibility', 'off', 'MarkerFaceColor', cmap(1,:),'MarkerEdgeAlpha',0.4);
%     hold on;
%     %
%     bar(2, mean(all_d2_CV), 'DisplayName', 'Mean \Deltacv 2 zstep from max', 'facecolor', cmap(2,:));
%     hold on;
%     errorbar(2, mean(all_d2_CV), std(all_d2_CV), 'DisplayName', 'Std of all (+/- 2 from best Z)', 'color', cmap(2,:));
%     hold on;
%     swarmchart(2*ones(size(all_d2_CV)), all_d2_CV, 'HandleVisibility', 'off', 'MarkerFaceColor', cmap(2,:),'MarkerEdgeAlpha',0.4);
%     hold on;
%     %
%     bar(3, mean(d_CV(idx_cv == 1)), 'DisplayName', 'Mean \Deltacv in samp. w/ 1 zstep diff.', 'facecolor', cmap(3,:));
%     hold on;
%     errorbar(3, mean(d_CV(idx_cv == 1)), std(d_CV(idx_cv == 1)), 'DisplayName', 'Std of all +/- 1 z diff', 'color', cmap(3,:));
%     hold on;
%     swarmchart(3*ones(size(d_CV(idx_cv == 1))), d_CV(idx_cv == 1), 'HandleVisibility', 'off', 'MarkerFaceColor', cmap(3,:),'MarkerEdgeAlpha',0.4);
%     hold on;
%     %
%     bar(4, mean(d_CV(idx_cv == 2)), 'DisplayName', 'Mean \Deltacv in samp. w/ 2 zstep diff', 'facecolor', cmap(4,:));
%     hold on;
%     errorbar(4, mean(d_CV(idx_cv == 2)), std(d_CV(idx_cv == 2)), 'DisplayName', 'Std of all +/- 2 z diff', 'color', cmap(4,:));
%     hold on;
%     swarmchart(4*ones(size(d_CV(idx_cv == 2))), d_CV(idx_cv == 2), 'HandleVisibility', 'off', 'MarkerFaceColor', cmap(4,:),'MarkerEdgeAlpha',0.4);
%     hold on;
% 
%     xlim([.5 4.5]);
%     ylim([0 .6]);
%     title(['Z w/ max Ph3 CV pre-seg. vs. Z w/ max ' channels{fluor_chan} ' CV post-seg.']);
%     xticks([1,2,3,4]);
%     xticklabels({'Fluor focus \newline 1 slice from max', 'Fluor focus \newline 2slices from max', ...
%         'Fluor focus \newline 1 slice from Ph3', 'Fluor focus \newline 2 slices from Ph3'})
%     ylabel([channels(fluor_chan) 'Difference in normalized CV ']);
%     grid on;
%     legend('location', 'best');
% 
% 
%     set(gcf,'units','points','position',[-1200,50,720*1.5,275*1.25*2]);
% 
%     saveas(f, ['postsegmentation_fluor_analysis\phase3_maxvar_z_vs_' ...
%         channels{fluor_chan} '_maxvar_z_maxcv_z_BAR.png']);
%     saveas(f, ['postsegmentation_fluor_analysis\phase3_maxvar_z_vs_' ...
%         channels{fluor_chan} '_maxvar_z_maxcv_z_BAR.fig']);
end

=======
%% Initialize
clc
clear all;

ip = exp_20211029_InitializeExp();

channels = ip.image.channels.name;
date = ['(' ip.exp.date(6:7) '/' ip.exp.date(9:10) ')'];

if ~exist([ip.exp.path '\presegmentation_fluor_analysis'], 'dir')
    mkdir([ip.exp.path '\presegmentation_fluor_analysis']);
end
if ~exist([ip.exp.path '\postsegmentation_fluor_analysis'], 'dir')
    mkdir([ip.exp.path '\postsegmentation_fluor_analysis']);
end

for c = channels
    if ~exist([ip.exp.path '\presegmentation_fluor_analysis\' c{1} '_CV_var_int_vs_zslice'], 'dir')
        mkdir([ip.exp.path '\presegmentation_fluor_analysis\' c{1} '_CV_var_int_vs_zslice']);
    end
    
    if ~exist([ip.exp.path '\postsegmentation_fluor_analysis\' c{1} '_CV_var_int_vs_zslice'], 'dir')
        mkdir([ip.exp.path '\postsegmentation_fluor_analysis\' c{1} '_CV_var_int_vs_zslice']);
    end  
end

clear c


%% Variance & Intensity & CV vs. z-slice PRE
fluor_chan = 2;

phase3_bestz = zeros(ip.exp.totalframes, 1);
fluor_bestvarz = zeros(ip.exp.totalframes,1);
fluor_bestCVz = zeros(ip.exp.totalframes,1);

color_used = {};

for c = 1:length(channels)
    %%% CV vs. z-slice PRE
    
    spl2znum = ip.image.zrange .* ...
        ones(1,max(ip.exp.sampleList)) ;

    y_min(1) = 1;
    y_min(2) = 1;
    y_min(3) = 1;

    for n_frame = 1:ip.exp.totalframes
        % default value
        imdigitnum_ = 2 ;
        current_spl = ip.exp.frm2spl(n_frame) ;
        if max(ip.sample{current_spl}.idx) <= 9
            imdigitnum_ = 1 ;
        end
        var_lst_cv = zeros(1,ip.image.zrange) ;
        var_lst_var = zeros(1,ip.image.zrange) ;
        var_lst_int = zeros(1,ip.image.zrange) ;

        for iz = 1: spl2znum(current_spl)
            PH3 = imread([ip.exp.path 'images\' ...
                ip.image.base_name '_' ...
                num2str(ip.exp.frm2spl(n_frame),'%03d') ...
                'xy' num2str(ip.exp.frm2img(n_frame),['%0' num2str(imdigitnum_ ,'%01d') 'd']) ...
                ['z' num2str(iz ,'%01d') 'c' num2str(c) '.tif'] ]);
            var_lst_cv(iz) = sqrt(var(double(PH3(:)))) / mean(double(PH3(:))) ;
            var_lst_var(iz) = var(double(PH3(:))) ;
            var_lst_int(iz) = sum(double(PH3(:))) ;
        end

        if y_min(1) > min(var_lst_cv/max(var_lst_cv))
            y_min(1) = min(var_lst_cv/max(var_lst_cv));
        end
        
        if y_min(2) > min(var_lst_var/max(var_lst_var))
            y_min(2) = min(var_lst_var/max(var_lst_var));
        end
        
        if y_min(3) > min(var_lst_int/max(var_lst_int))
            y_min(3) = min(var_lst_int/max(var_lst_int));
        end
        
        if c == fluor_chan
            [~,fluor_bestCVz(n_frame)] = max(var_lst_cv);
            [~,fluor_bestvarz(n_frame)] = max(var_lst_var);
        end
        
        if strcmp(channels{c}, 'phase3')
            [~,phase3_bestz(n_frame)] = max(var_lst_var);
        end

        figure(ip.exp.frm2spl(n_frame) + (c-1)*max(ip.exp.frm2spl));
        subplot(1,3,1);
        hold on;
        if c == 1
            color_used{n_frame} = which_color(gca);
        end
        plot(1: spl2znum(current_spl), var_lst_cv/max(var_lst_cv), 'DisplayName', ...
            ['XY position ', num2str(ip.exp.frm2img(n_frame))], 'Color', color_used{n_frame});
        xlabel('z slice')
        ylabel([channels{c} ' pixel CV / max ' channels{c} ' CV']);
        title(['Normalized ' channels{c} ' pixel CV vs. z-slice'])
        
        subplot(1,3,2);
        hold on;
        %color_used{n_frame} = which_color(gca);
        plot(1: spl2znum(current_spl), var_lst_var/max(var_lst_var), 'DisplayName', ...
            ['XY position ', num2str(ip.exp.frm2img(n_frame))], 'Color', color_used{n_frame});
        xlabel('z slice')
        ylabel([channels{c} ' pixel variance / max ' channels{c} ' variance']);
        title(['Normalized ' channels{c} ' pixel var. vs. z-slice'])
        
        subplot(1,3,3);
        hold on;
        %color_used{n_frame} = which_color(gca);
        plot(1: spl2znum(current_spl), var_lst_int/max(var_lst_int), 'DisplayName', ...
            ['XY position ', num2str(ip.exp.frm2img(n_frame))], 'Color', color_used{n_frame});
        xlabel('z slice')
        ylabel([channels{c} ' pixel intensity / max ' channels{c} ' intensity']);
        title(['Normalized ' channels{c} ' pixel int. vs. z-slice'])        
    end
    %{
    %%% Variance vs. z-slice PRE
    
    spl2znum = ip.image.zrange .* ...
        ones(1,max(ip.exp.sampleList)) ;

    color_used = {};
    y_min(2) = 1;

    for n_frame = 1:ip.exp.totalframes
        % default value
        imdigitnum_ = 2 ;
        current_spl = ip.exp.frm2spl(n_frame) ;
        if max(ip.sample{current_spl}.idx) <= 9
            imdigitnum_ = 1 ;
        end
        var_lst = zeros(1,ip.image.zrange) ;

        for iz = 1: spl2znum(current_spl)
            PH3 = imread([ip.exp.path 'images\' ...
                ip.image.base_name '_' ...
                num2str(ip.exp.frm2spl(n_frame),'%03d') ...
                'xy' num2str(ip.exp.frm2img(n_frame),['%0' num2str(imdigitnum_ ,'%01d') 'd']) ...
                ['z' num2str(iz ,'%01d') 'c' num2str(c) '.tif'] ]);
            var_lst(iz) = var(double(PH3(:))) ;
        end

        if y_min(2) > min(var_lst/max(var_lst))
            y_min(2) = min(var_lst/max(var_lst));
        end
        
        if c == fluor_chan
            [~,fluor_bestvarz(n_frame)] = max(var_lst);
        end
        
        if strcmp(channels{c}, 'phase3')
            [~,phase3_bestz(n_frame)] = max(var_lst);
        end

        figure(ip.exp.frm2spl(n_frame) + (c-1)*max(ip.exp.frm2spl));
        subplot(1,3,2);
        hold on;
        color_used{n_frame} = which_color(gca);
        plot(1: spl2znum(current_spl), var_lst/max(var_lst), 'DisplayName', ['XY position ', num2str(ip.exp.frm2img(n_frame))], 'Color', color_used{n_frame});
        xlabel('z slice')
        ylabel([channels{c} ' pixel variance / max ' channels{c} ' variance']);
        title(['Normalized ' channels{c} ' pixel var. vs. z-slice'])
    end

    %%% Intensity vs. z-slice PRE

    spl2znum = ip.image.zrange .* ...
        ones(1,max(ip.exp.sampleList)) ;

    color_used = {};
    y_min(3) = 1;

    for n_frame = 1:ip.exp.totalframes

        % default value
        imdigitnum_ = 2 ;
        current_spl = ip.exp.frm2spl(n_frame) ;
        if max(ip.sample{current_spl}.idx) <= 9
            imdigitnum_ = 1 ;
        end
        int_lst = zeros(1,ip.image.zrange) ;

        for iz = 1: spl2znum(current_spl)
            PH3 = imread([ip.exp.path 'images\' ...
                ip.image.base_name '_' ...
                num2str(ip.exp.frm2spl(n_frame),'%03d') ...
                'xy' num2str(ip.exp.frm2img(n_frame),['%0' num2str(imdigitnum_ ,'%01d') 'd']) ...
                ['z' num2str(iz ,'%01d') 'c' num2str(c) '.tif'] ]);
            int_lst(iz) = sum(double(PH3(:))) ;
        end

        if y_min(3) > min(int_lst/max(int_lst))
            y_min(3) = min(int_lst/max(int_lst));
        end

        figure(ip.exp.frm2spl(n_frame) + (c-1)*max(ip.exp.frm2spl));
        subplot(1,3,3);
        hold on;
        color_used{n_frame} = which_color(gca);
        plot(1: spl2znum(current_spl), int_lst/max(int_lst), 'DisplayName', ['XY position ', num2str(ip.exp.frm2img(n_frame))], 'Color', color_used{n_frame});
        xlabel('z slice')
        ylabel([channels{c} ' pixel intensity / max ' channels{c} ' intensity']);
        title(['Normalized ' channels{c} ' pixel int. vs. z-slice'])
    end
    %}
    %%% Saving figures of current channel
    
    range_to_process = unique(ip.exp.frm2spl(1:ip.exp.totalframes));
    for i = transpose(range_to_process)
        figure(i + (c-1)*max(ip.exp.frm2spl));
        set(gcf,'units','points','position',[-1000,250,720,326]);
        sgtitle([date ' Sample ', num2str(i), ': ',ip.sample{1,i}.name]);
        for ii = 1:3
            subplot(1,3,ii);
            legend('location', 'best');
            ylim([y_min(ii) 1]);
            
            grid on;
        end

        saveas(figure(i + (c-1)*max(ip.exp.frm2spl)), ['presegmentation_fluor_analysis\' channels{c} '_CV_var_int_vs_zslice\Sample_', ...
            num2str(i, '%02d'), '__Norm' channels{c} 'Pixel__CV_var_int_vs_z.png']);
        saveas(figure(i + (c-1)*max(ip.exp.frm2spl)), ['presegmentation_fluor_analysis\' channels{c} '_CV_var_int_vs_zslice\Sample_', ...
            num2str(i, '%02d'), '__Norm' channels{c} 'Pixel__CV_var_int_vs_z.fig']);
        close(figure(i + (c-1)*max(ip.exp.frm2spl)));
    end
end

%%% Best phase3 variance vs. best FLUOR variance
cmap = hsv(max(ip.exp.frm2spl));

figure ();
subplot(2,2,1);
for i = 1:max(ip.exp.frm2spl)
    xy = [phase3_bestz(nonzeros(ip.exp.splimg2frm(i,:))) - max(ip.exp.frm2spl)/2 * .05 + i * .05,...
        fluor_bestvarz(nonzeros(ip.exp.splimg2frm(i,:))) - max(ip.exp.frm2spl)/2 * .05 + i * .05];
    [C,ia,ic] = unique(xy, 'rows');
    a_counts = accumarray(ic,1);

    scatter(C(:,1), C(:,2), a_counts*5, 'markeredgecolor', 'none',...
        'markerfacecolor', cmap(i,:), 'DisplayName', ['Sample ', num2str(i)]);

%     scatter(unique(phase3_bestz(nonzeros(ip.exp.splimg2frm(i,:))) - max(ip.exp.frm2spl)/2 * .05 + i * .05), ...
%         unique(fluor_bestvarz(nonzeros(ip.exp.splimg2frm(i,:))) - max(ip.exp.frm2spl)/2 * .05 + i * .05), ...
%         a_counts/max(a_counts) * 15, 'markerfacecolor', cmap(i,:), 'DisplayName', ['Sample ', num2str(i)]);
    hold on;
end
plot(0:8, 0:8, 'Color', 'black', 'DisplayName', 'Slope = 1');
xlim([0 8]);
ylim([0 8]);
title(['Z w/ max Ph3 var. vs. Z w/ max ' channels{fluor_chan} ' var.']);
xlabel('Phase 3 highest variance z slice');
ylabel([channels(fluor_chan) ' highest variance z slice']);
grid on;
legend('location', 'best');

subplot(2,2,3);
xy = [phase3_bestz(nonzeros(ip.exp.splimg2frm(:))), fluor_bestvarz(nonzeros(ip.exp.splimg2frm(:)))];
[C,ia,ic] = unique(xy, 'rows');
a_counts = accumarray(ic,1);
scatter(C(:,1), C(:,2), a_counts*5, 'markeredgecolor', 'none',...
    'markerfacecolor', 'black', 'DisplayName', ['All samples together']);
hold on;
plot(0:8, 0:8, 'Color', 'black', 'DisplayName', 'Slope = 1');
xlim([0 8]);
ylim([0 8]);
title(['Z w/ max Ph3 var. vs. Z w/ max ' channels{fluor_chan} ' var.']);
xlabel('Phase 3 highest variance z slice');
ylabel([channels(fluor_chan) ' highest variance z slice']);
grid on;
legend('location', 'best');


subplot(2,2,2);
for i = 1:max(ip.exp.frm2spl)
    xy = [phase3_bestz(nonzeros(ip.exp.splimg2frm(i,:))) - max(ip.exp.frm2spl)/2 * .05 + i * .05,...
        fluor_bestCVz(nonzeros(ip.exp.splimg2frm(i,:))) - max(ip.exp.frm2spl)/2 * .05 + i * .05];
    [C,ia,ic] = unique(xy, 'rows');
    a_counts = accumarray(ic,1);

    scatter(C(:,1), C(:,2), a_counts/max(a_counts) * 50, 'markeredgecolor', 'none',...
        'markerfacecolor', cmap(i,:), 'DisplayName', ['Sample ', num2str(i)]);
    
%     scatter(phase3_bestz(nonzeros(ip.exp.splimg2frm(i,:))) - max(ip.exp.frm2spl)/2 * .05 + i * .05, ...
%         fluor_bestCVz(nonzeros(ip.exp.splimg2frm(i,:))) - max(ip.exp.frm2spl)/2 * .05 + i * .05, ...
%         a_counts/max(a_counts) * 15, 'markerfacecolor', cmap(i,:), 'DisplayName', ['Sample ', num2str(i)]);
    hold on;
end
plot(0:8, 0:8, 'Color', 'black', 'DisplayName', 'Slope = 1');
xlim([0 8]);
ylim([0 8]);
title(['Z w/ max Ph3 CV vs. Z w/ max ' channels{fluor_chan} ' CV']);
xlabel('Phase 3 highest CV z slice');
ylabel([channels{fluor_chan} ' highest CV z slice']);
grid on;
legend('location', 'best');

subplot(2,2,4);
xy = [phase3_bestz(nonzeros(ip.exp.splimg2frm(:))), fluor_bestCVz(nonzeros(ip.exp.splimg2frm(:)))];
[C,ia,ic] = unique(xy, 'rows');
a_counts = accumarray(ic,1);
scatter(C(:,1), C(:,2), a_counts*5, 'markeredgecolor', 'none',...
    'markerfacecolor', 'black', 'DisplayName', ['All samples together']);
hold on;
plot(0:8, 0:8, 'Color', 'black', 'DisplayName', 'Slope = 1');
xlim([0 8]);
ylim([0 8]);
title(['Z w/ max Ph3 CV vs. Z w/ max ' channels{fluor_chan} ' CV']);
xlabel('Phase 3 highest CV z slice');
ylabel([channels{fluor_chan} ' highest CV z slice']);
grid on;
legend('location', 'best');

set(gcf,'units','points','position',[-1200,50,720*1.5,275*1.25*2]);

saveas(figure(i + (c-1)*max(ip.exp.frm2spl)), ['presegmentation_fluor_analysis\phase3_maxvar_z_vs_' ...
    channels{fluor_chan} '_maxvar_z_maxcv_z.png']);
saveas(figure(i + (c-1)*max(ip.exp.frm2spl)), ['presegmentation_fluor_analysis\phase3_maxvar_z_vs_' ...
    channels{fluor_chan} '_maxvar_z_maxcv_z.fig']);


%% Variance & Intensity & CV vs. z-slice POST (whole frame focus comparison)
fluor_chan = 2;
neg_samp = 8;
wholexy = 0;

phase3_bestz = zeros(ip.exp.totalframes, 1);
phase3_bestzPRE = zeros(ip.exp.totalframes, 1);
fluor_bestvarz = zeros(ip.exp.totalframes,1);
fluor_bestCVz = zeros(ip.exp.totalframes,1);

[num_sam, num_img] = size(ip.exp.splimg2frm);

fluor_varz = zeros([num_sam - 1, num_img, ip.image.zrange]);
fluor_CVz = zeros([num_sam - 1, num_img, ip.image.zrange]);

color_used = {};

for c = 1:length(channels)
    %%% CV vs. z-slice PRE
    
    spl2znum = ip.image.zrange .* ...
        ones(1,max(ip.exp.sampleList)) ;

    y_min(1) = 1;
    y_min(2) = 1;
    y_min(3) = 1;
    
    for n_frame = 1:ip.exp.totalframes
        n_image = ip.exp.frm2img(n_frame);
        %figure(ip.exp.frm2spl(n_frame));    
        %set(gcf,'units','points','position',[-1000,31,720,715]);
        
        if n_image == 1
            fprintf(1, newline);
            fprintf(1,['Sample ' num2str(ip.exp.frm2spl(n_frame)) ' of ' ...
                num2str(max(ip.exp.sampleList)) '.' newline]);
            progress_2 = [];
        end
        for d_=1:1:size(progress_2,2) ; fprintf(1,'\b') ; end
        progress_2 = [sprintf('\t') 'Image ' num2str(n_image) ' in ' ...
            num2str(min(ip.sample{ip.exp.frm2spl(n_frame)}.idx)) '-' ...
            num2str(max(ip.sample{ip.exp.frm2spl(n_frame)}.idx)) '.'] ;
        fprintf(1,progress_2) ;

        load([ip.seg.dir ip.seg.base_name 'seg' num2str(n_frame,'%03d') '.mat'],'LcFull') ;

        % default value
        imdigitnum_ = 2 ;
        current_spl = ip.exp.frm2spl(n_frame) ;
        if max(ip.sample{current_spl}.idx) <= 9
            imdigitnum_ = 1 ;
        end
        var_lst = zeros(1,ip.image.zrange) ;

        CellMap = LcFull;
        nnn=double(LcFull);
        N=max(nnn(:));
        
        cv_temp = zeros(spl2znum(current_spl),N);
        var_temp = zeros(spl2znum(current_spl),N);
        int_temp = zeros(spl2znum(current_spl),N);
        
        for index=1:N
            for iz = 1: spl2znum(current_spl)
                PH3 = imread([ip.exp.path 'images\' ...
                    ip.image.base_name '_' ...
                    num2str(ip.exp.frm2spl(n_frame),'%03d') ...
                    'xy' num2str(ip.exp.frm2img(n_frame),['%0' num2str(imdigitnum_ ,'%01d') 'd']) ...
                    ['z' num2str(iz ,'%01d') 'c' num2str(c) '.tif'] ]);

                Cell = immultiply( CellMap == index, PH3);
                int_per_pix = Cell(CellMap == index);

                cv_temp(iz,index) = sqrt(var(double(int_per_pix(:)))) / mean(double(int_per_pix(:))) ;
                var_temp(iz,index) = var(double(int_per_pix(:))) ;
                int_temp(iz,index) = sum(double(int_per_pix(:))) ;
            end
        end
        
        % sum whole image, keep z seperate
        var_lst_cv = sum(transpose(cv_temp)) ;
        var_lst_var = sum(transpose(var_temp)) ;
        var_lst_int = sum(transpose(int_temp)) ;
    
        %for plot limits
        if y_min(1) > min(var_lst_cv/max(var_lst_cv))
            y_min(1) = min(var_lst_cv/max(var_lst_cv));
        end
        
        if y_min(2) > min(var_lst_var/max(var_lst_var))
            y_min(2) = min(var_lst_var/max(var_lst_var));
        end
        
        if y_min(3) > min(var_lst_int/max(var_lst_int))
            y_min(3) = min(var_lst_int/max(var_lst_int));
        end
        
        %for specific channel of interest only, get var and cv normalized
        %for bar plots later
        if c == fluor_chan
            [~,fluor_bestCVz(n_frame)] = max(var_lst_cv);
            [~,fluor_bestvarz(n_frame)] = max(var_lst_var);
            if ip.exp.frm2spl(n_frame) < neg_samp
                fluor_varz(ip.exp.frm2spl(n_frame),ip.exp.frm2img(n_frame), :) = var_lst_var/max(var_lst_var);
                fluor_CVz(ip.exp.frm2spl(n_frame),ip.exp.frm2img(n_frame), :) = var_lst_cv/max(var_lst_cv);
            elseif ip.exp.frm2spl(n_frame) > neg_samp
                fluor_varz(ip.exp.frm2spl(n_frame)-1,ip.exp.frm2img(n_frame), :) = var_lst_var/max(var_lst_var);
                fluor_CVz(ip.exp.frm2spl(n_frame)-1,ip.exp.frm2img(n_frame), :) = var_lst_cv/max(var_lst_cv);
            end
        end
        %go over the whole image only (rather than cell by cell)
        if strcmp(channels{c}, 'phase3')
            [~,phase3_bestz(n_frame)] = max(var_lst_var);
            
            var_temp = zeros(spl2znum(current_spl),1);
            for iz = 1: spl2znum(current_spl)
                PH3 = imread([ip.exp.path 'images\' ...
                    ip.image.base_name '_' ...
                    num2str(ip.exp.frm2spl(n_frame),'%03d') ...
                    'xy' num2str(ip.exp.frm2img(n_frame),['%0' num2str(imdigitnum_ ,'%01d') 'd']) ...
                    ['z' num2str(iz ,'%01d') 'c' num2str(c) '.tif'] ]);
                
                cv_temp(iz) = sqrt(var(double(PH3(:)))) / mean(double(PH3(:))) ;
                var_temp(iz) = var(double(PH3(:))) ;
                int_temp(iz) = sum(double(PH3(:))) ;
            end
            if wholexy
                var_lst_cv = cv_temp ;
                var_lst_var = var_temp ;
                var_lst_int = int_temp ;
            end
            [~,phase3_bestzPRE(n_frame)] = max(var_temp);
        end

        figure(ip.exp.frm2spl(n_frame) + (c-1)*max(ip.exp.frm2spl));
        subplot(1,3,1);
        hold on;
        %phase3 defines what colors to use
        if c == 1
            color_used{n_frame} = which_color(gca);
        end

            plot(1: spl2znum(current_spl), var_lst_cv/max(var_lst_cv), 'DisplayName', ...
                ['XY position ', num2str(ip.exp.frm2img(n_frame))], 'Color', color_used{n_frame});
            xlabel('z slice')
            ylabel([channels{c} ' pixel CV / max ' channels{c} ' CV']);
            title(['Normalized ' channels{c} ' pixel CV vs. z-slice'])

            subplot(1,3,2);
            hold on;
            %color_used{n_frame} = which_color(gca);
            plot(1: spl2znum(current_spl), var_lst_var/max(var_lst_var), 'DisplayName', ...
                ['XY position ', num2str(ip.exp.frm2img(n_frame))], 'Color', color_used{n_frame});
            xlabel('z slice')
            ylabel([channels{c} ' pixel variance / max ' channels{c} ' variance']);
            title(['Normalized ' channels{c} ' pixel var. vs. z-slice'])

            subplot(1,3,3);
            hold on;
            %color_used{n_frame} = which_color(gca);
            plot(1: spl2znum(current_spl), var_lst_int/max(var_lst_int), 'DisplayName', ...
                ['XY position ', num2str(ip.exp.frm2img(n_frame))], 'Color', color_used{n_frame});
            xlabel('z slice')
            ylabel([channels{c} ' pixel intensity / max ' channels{c} ' intensity']);
            title(['Normalized ' channels{c} ' pixel int. vs. z-slice'])
    end
    
    %%% Saving figures of current channel
    
    range_to_process = unique(ip.exp.frm2spl(1:ip.exp.totalframes));
    for i = transpose(range_to_process)
        figure(i + (c-1)*max(ip.exp.frm2spl));
        set(gcf,'units','points','position',[-1000,250,720,326]);
        if strcmp(channels{c}, 'phase3') && wholexy
            sgtitle([date ' Sample ', num2str(i), ': ',ip.sample{1,i}.name, ' [Analysis over whole XY]']);
        else
            sgtitle([date ' Sample ', num2str(i), ': ',ip.sample{1,i}.name, ' [Analysis over each cell]']);
        end
        for ii = 1:3
            subplot(1,3,ii);
            legend('location', 'best');
            ylim([y_min(ii) 1]);
            
            grid on;
        end

        saveas(figure(i + (c-1)*max(ip.exp.frm2spl)), ...
            ['postsegmentation_fluor_analysis\' channels{c} '_CV_var_int_vs_zslice\Sample_', ...
            num2str(i, '%02d'), '__Norm' channels{c} 'Pixel__CV_var_int_vs_z.png']);
        saveas(figure(i + (c-1)*max(ip.exp.frm2spl)), ...
            ['postsegmentation_fluor_analysis\' channels{c} '_CV_var_int_vs_zslice\Sample_', ...
            num2str(i, '%02d'), '__Norm' channels{c} 'Pixel__CV_var_int_vs_z.fig']);
        
        close(figure(i + (c-1)*max(ip.exp.frm2spl)));
    end
end

%%% Best phase3 variance vs. best FLUOR variance
cmap = hsv(max(ip.exp.frm2spl));
phase3_bestzPRE = phase3_bestz;

f = figure ();
subplot(2,2,1);
for i = 1:max(ip.exp.frm2spl)
    xy = [phase3_bestzPRE(nonzeros(ip.exp.splimg2frm(i,:))) - max(ip.exp.frm2spl)/2 * .05 + i * .05,...
        fluor_bestvarz(nonzeros(ip.exp.splimg2frm(i,:))) - max(ip.exp.frm2spl)/2 * .05 + i * .05];
    [C,ia,ic] = unique(xy, 'rows');
    a_counts = accumarray(ic,1);

    scatter(C(:,1), C(:,2), a_counts*5, 'markeredgecolor', 'none',...
        'markerfacecolor', cmap(i,:), 'DisplayName', ['Sample ', num2str(i)]);

%     scatter(unique(phase3_bestz(nonzeros(ip.exp.splimg2frm(i,:))) - max(ip.exp.frm2spl)/2 * .05 + i * .05), ...
%         unique(fluor_bestvarz(nonzeros(ip.exp.splimg2frm(i,:))) - max(ip.exp.frm2spl)/2 * .05 + i * .05), ...
%         a_counts/max(a_counts) * 15, 'markerfacecolor', cmap(i,:), 'DisplayName', ['Sample ', num2str(i)]);
    hold on;
end
plot(0:8, 0:8, 'Color', 'black', 'DisplayName', 'Slope = 1');
xlim([0 8]);
ylim([0 8]);
title(['Z w/ max Ph3 var. pre-seg. vs. Z w/ max ' channels{fluor_chan} ' var. post-seg.']);
xlabel('Phase 3 highest variance z slice');
ylabel([channels(fluor_chan) ' highest variance z slice']);
grid on;
legend('location', 'best');

subplot(2,2,3);
xy = [phase3_bestzPRE(nonzeros(ip.exp.splimg2frm(:))), fluor_bestvarz(nonzeros(ip.exp.splimg2frm(:)))];
[C,ia,ic] = unique(xy, 'rows');
a_counts = accumarray(ic,1);
scatter(C(:,1), C(:,2), a_counts*5, 'markeredgecolor', 'none',...
    'markerfacecolor', 'black', 'DisplayName', ['All samples together']);
hold on;
plot(0:8, 0:8, 'Color', 'black', 'DisplayName', 'Slope = 1');
xlim([0 8]);
ylim([0 8]);
title(['Z w/ max Ph3 var. pre-seg. vs. Z w/ max ' channels{fluor_chan} ' var. post-seg']);
xlabel('Phase 3 highest variance z slice');
ylabel([channels(fluor_chan) ' highest variance z slice']);
grid on;
legend('location', 'best');


subplot(2,2,2);
for i = 1:max(ip.exp.frm2spl)
    xy = [phase3_bestzPRE(nonzeros(ip.exp.splimg2frm(i,:))) - max(ip.exp.frm2spl)/2 * .05 + i * .05,...
        fluor_bestCVz(nonzeros(ip.exp.splimg2frm(i,:))) - max(ip.exp.frm2spl)/2 * .05 + i * .05];
    [C,ia,ic] = unique(xy, 'rows');
    a_counts = accumarray(ic,1);

    scatter(C(:,1), C(:,2), a_counts/max(a_counts) * 50, 'markeredgecolor', 'none',...
        'markerfacecolor', cmap(i,:), 'DisplayName', ['Sample ', num2str(i)]);
    
%     scatter(phase3_bestz(nonzeros(ip.exp.splimg2frm(i,:))) - max(ip.exp.frm2spl)/2 * .05 + i * .05, ...
%         fluor_bestCVz(nonzeros(ip.exp.splimg2frm(i,:))) - max(ip.exp.frm2spl)/2 * .05 + i * .05, ...
%         a_counts/max(a_counts) * 15, 'markerfacecolor', cmap(i,:), 'DisplayName', ['Sample ', num2str(i)]);
    hold on;
end
plot(0:8, 0:8, 'Color', 'black', 'DisplayName', 'Slope = 1');
xlim([0 8]);
ylim([0 8]);
title(['Z w/ max Ph3 var pre-seg. vs. Z w/ max ' channels{fluor_chan} ' CV post-seg.']);
xlabel('Phase 3 highest CV z slice');
ylabel([channels{fluor_chan} ' highest CV z slice']);
grid on;
legend('location', 'best');

subplot(2,2,4);
xy = [phase3_bestzPRE(nonzeros(ip.exp.splimg2frm(:))), fluor_bestCVz(nonzeros(ip.exp.splimg2frm(:)))];
[C,ia,ic] = unique(xy, 'rows');
a_counts = accumarray(ic,1);
scatter(C(:,1), C(:,2), a_counts*5, 'markeredgecolor', 'none',...
    'markerfacecolor', 'black', 'DisplayName', ['All samples together']);
hold on;
plot(0:8, 0:8, 'Color', 'black', 'DisplayName', 'Slope = 1');
xlim([0 8]);
ylim([0 8]);
title(['Z w/ max Ph3 var pre-seg vs. Z w/ max ' channels{fluor_chan} ' CV post-seg.']);
xlabel('Phase 3 highest CV z slice');
ylabel([channels{fluor_chan} ' highest CV z slice']);
grid on;
legend('location', 'best');

set(gcf,'units','points','position',[-1200,50,720*1.5,275*1.25*2]);

saveas(f, ['postsegmentation_fluor_analysis\phase3_maxvar_z_vs_' ...
    channels{fluor_chan} '_maxvar_z_maxcv_z.png']);
saveas(f, ['postsegmentation_fluor_analysis\phase3_maxvar_z_vs_' ...
    channels{fluor_chan} '_maxvar_z_maxcv_z.fig']);


%%% Best phase3 variance vs. best FLUOR variance bar of difference
%{
pick_z_fluorvar = transpose(fluor_bestvarz(ip.exp.frm2spl ~= 8));
[all_x, all_y] = ind2sub(size(fluor_varz(:,:,1)), 1:numel(fluor_varz(:,:,1)));
all_x = all_x(fluor_varz(:,:,1) ~= 0);
all_y = all_y(fluor_varz(:,:,1) ~= 0);
best_var = fluor_varz(sub2ind(size(fluor_varz), all_x, all_y, pick_z_fluorvar));

best_var_t = best_var(pick_z_fluorvar > 2 & pick_z_fluorvar < ip.image.zrange-1);
pick_z_fluorvar_t = pick_z_fluorvar(pick_z_fluorvar > 2 & pick_z_fluorvar < ip.image.zrange-1);
all_d1_var = [abs(best_var - fluor_varz(sub2ind(size(fluor_varz), all_x, all_y, pick_z_fluorvar+1))), ...
    abs(best_var - fluor_varz(sub2ind(size(fluor_varz), all_x, all_y, pick_z_fluorvar-1)))];
all_d2_var = [abs(best_var_t - fluor_varz(sub2ind(size(fluor_varz), ...
    all_x(pick_z_fluorvar > 2 & pick_z_fluorvar < ip.image.zrange-1), ...
    all_y(pick_z_fluorvar > 2 & pick_z_fluorvar < ip.image.zrange-1), pick_z_fluorvar_t+2))), ...
    abs(best_var_t - fluor_varz(sub2ind(size(fluor_varz), ...
    all_x(pick_z_fluorvar > 2 & pick_z_fluorvar < ip.image.zrange-1), ...
    all_y(pick_z_fluorvar > 2 & pick_z_fluorvar < ip.image.zrange-1), pick_z_fluorvar_t-2)))];

pick_z_ph3 = transpose(phase3_bestzPRE(ip.exp.frm2spl ~= 8));
ph3_var = fluor_varz(sub2ind(size(fluor_varz), all_x, all_y, pick_z_ph3));
d_var = abs(best_var - ph3_var);
idx_var = abs(pick_z_fluorvar - pick_z_ph3);
%
pick_z_fluorcv = transpose(fluor_bestCVz(ip.exp.frm2spl ~= 8));
[all_x, all_y] = ind2sub(size(fluor_CVz(:,:,1)), 1:numel(fluor_CVz(:,:,1)));
all_x = all_x(fluor_CVz(:,:,1) ~= 0);
all_y = all_y(fluor_CVz(:,:,1) ~=0);
best_CV = fluor_CVz(sub2ind(size(fluor_CVz), all_x, all_y, pick_z_fluorcv));

best_CV_t = best_CV(pick_z_fluorcv > 2 & pick_z_fluorcv < ip.image.zrange-1);
pick_z_fluorcv_t = pick_z_fluorcv(pick_z_fluorcv > 2 & pick_z_fluorcv < ip.image.zrange-1);
all_d1_CV = [abs(best_CV(pick_z_fluorcv > 1 & pick_z_fluorcv < ip.image.zrange) - fluor_CVz(sub2ind(size(fluor_CVz), ...
    all_x(pick_z_fluorcv > 1 & pick_z_fluorcv < ip.image.zrange), ...
    all_y(pick_z_fluorcv > 1 & pick_z_fluorcv < ip.image.zrange), ...
    pick_z_fluorcv(pick_z_fluorcv > 1 & pick_z_fluorcv < ip.image.zrange)+1))), ...
    abs(best_CV(pick_z_fluorcv > 1 & pick_z_fluorcv < ip.image.zrange) - fluor_CVz(sub2ind(size(fluor_CVz), ...
    all_x(pick_z_fluorcv > 1 & pick_z_fluorcv < ip.image.zrange), ...
    all_y(pick_z_fluorcv > 1 & pick_z_fluorcv < ip.image.zrange), ...
    pick_z_fluorcv(pick_z_fluorcv > 1 & pick_z_fluorcv < ip.image.zrange)-1)))];
all_d2_CV = [abs(best_CV_t - fluor_CVz(sub2ind(size(fluor_CVz), ...
    all_x(pick_z_fluorcv > 2 & pick_z_fluorcv < ip.image.zrange-1), ...
    all_y(pick_z_fluorcv > 2 & pick_z_fluorcv < ip.image.zrange-1), pick_z_fluorcv_t+2))), ...
    abs(best_CV_t - fluor_CVz(sub2ind(size(fluor_CVz), ...
    all_x(pick_z_fluorcv > 2 & pick_z_fluorcv < ip.image.zrange-1), ...
    all_y(pick_z_fluorcv > 2 & pick_z_fluorcv < ip.image.zrange-1), pick_z_fluorcv_t-2)))];

ph3_CV = fluor_CVz(sub2ind(size(fluor_CVz), all_x, all_y, pick_z_ph3));
d_CV = abs(best_CV - ph3_CV);
idx_cv = abs(pick_z_fluorcv - pick_z_ph3);
%
f = figure ();
subplot(1,2,1);
bar(1, mean(all_d1_var), 'DisplayName', 'Mean \Deltavar. 1 zstep from max', 'facecolor', cmap(1,:));
hold on;
errorbar(1, mean(all_d1_var), std(all_d1_var), 'DisplayName', 'Std of all (+/- 1 from best Z)', 'color', cmap(1,:));
hold on;
swarmchart(1*ones(size(all_d1_var)), all_d1_var, 'HandleVisibility', 'off', 'MarkerFaceColor', cmap(1,:),'MarkerEdgeAlpha',0.4);
hold on;
%
bar(2, mean(all_d2_var), 'DisplayName', 'Mean \Deltavar. 2 zstep from max', 'facecolor', cmap(2,:));
hold on;
errorbar(2, mean(all_d2_var), std(all_d2_var), 'DisplayName', 'Std of all (+/- 2 from best Z)', 'color', cmap(2,:));
hold on;
swarmchart(2*ones(size(all_d2_var)), all_d2_var, 'HandleVisibility', 'off', 'MarkerFaceColor', cmap(2,:),'MarkerEdgeAlpha',0.4);
hold on;
%
bar(3, mean(d_var(idx_var == 1)), 'DisplayName', 'Mean \Deltavar. in samp. w/ 1 zstep diff.', 'facecolor', cmap(3,:));
hold on;
errorbar(3, mean(d_var(idx_var == 1)), std(d_var(idx_var == 1)), 'DisplayName', 'Std of all +/- 1 z diff', 'color', cmap(3,:));
hold on;
swarmchart(3*ones(size(d_var(idx_var == 1))), d_var(idx_var == 1), 'HandleVisibility', 'off', 'MarkerFaceColor', cmap(3,:),'MarkerEdgeAlpha',0.4);
hold on;
%
bar(4, mean(d_var(idx_var == 2)), 'DisplayName', 'Mean \Deltavar. in samp. w/ 2 zstep diff', 'facecolor', cmap(4,:));
hold on;
errorbar(4, mean(d_var(idx_var == 2)), std(d_var(idx_var == 2)), 'DisplayName', 'Std of all +/- 2 z diff', 'color', cmap(4,:));
hold on;
swarmchart(4*ones(size(d_var(idx_var == 2))), d_var(idx_var == 2), 'HandleVisibility', 'off', 'MarkerFaceColor', cmap(4,:),'MarkerEdgeAlpha',0.4);
hold on;

xlim([.5 4.5]);
ylim([0 .6]);
title(['Z w/ max Ph3 var. pre-seg. vs. Z w/ max ' channels{fluor_chan} ' var. post-seg.']);
xticks([1,2,3,4]);
xticklabels({'Fluor focus \newline 1 slice from max', 'Fluor focus \newline 2slices from max', ...
    'Fluor focus \newline 1 slice from Ph3', 'Fluor focus \newline 2 slices from Ph3'})
ylabel([channels(fluor_chan) ' difference in normalized variance']);
grid on;
legend('location', 'best');


subplot(1,2,2);
bar(1, mean(all_d1_CV), 'DisplayName', 'Mean \Deltacv 1 zstep from max', 'facecolor', cmap(1,:));
hold on;
errorbar(1, mean(all_d1_CV), std(all_d1_CV), 'DisplayName', 'Std of all (+/- 1 from best Z)', 'color', cmap(1,:));
hold on;
swarmchart(1*ones(size(all_d1_CV)), all_d1_CV, 'HandleVisibility', 'off', 'MarkerFaceColor', cmap(1,:),'MarkerEdgeAlpha',0.4);
hold on;
%
bar(2, mean(all_d2_CV), 'DisplayName', 'Mean \Deltacv 2 zstep from max', 'facecolor', cmap(2,:));
hold on;
errorbar(2, mean(all_d2_CV), std(all_d2_CV), 'DisplayName', 'Std of all (+/- 2 from best Z)', 'color', cmap(2,:));
hold on;
swarmchart(2*ones(size(all_d2_CV)), all_d2_CV, 'HandleVisibility', 'off', 'MarkerFaceColor', cmap(2,:),'MarkerEdgeAlpha',0.4);
hold on;
%
bar(3, mean(d_CV(idx_cv == 1)), 'DisplayName', 'Mean \Deltacv in samp. w/ 1 zstep diff.', 'facecolor', cmap(3,:));
hold on;
errorbar(3, mean(d_CV(idx_cv == 1)), std(d_CV(idx_cv == 1)), 'DisplayName', 'Std of all +/- 1 z diff', 'color', cmap(3,:));
hold on;
swarmchart(3*ones(size(d_CV(idx_cv == 1))), d_CV(idx_cv == 1), 'HandleVisibility', 'off', 'MarkerFaceColor', cmap(3,:),'MarkerEdgeAlpha',0.4);
hold on;
%
bar(4, mean(d_CV(idx_cv == 2)), 'DisplayName', 'Mean \Deltacv in samp. w/ 2 zstep diff', 'facecolor', cmap(4,:));
hold on;
errorbar(4, mean(d_CV(idx_cv == 2)), std(d_CV(idx_cv == 2)), 'DisplayName', 'Std of all +/- 2 z diff', 'color', cmap(4,:));
hold on;
swarmchart(4*ones(size(d_CV(idx_cv == 2))), d_CV(idx_cv == 2), 'HandleVisibility', 'off', 'MarkerFaceColor', cmap(4,:),'MarkerEdgeAlpha',0.4);
hold on;

xlim([.5 4.5]);
ylim([0 .6]);
title(['Z w/ max Ph3 CV pre-seg. vs. Z w/ max ' channels{fluor_chan} ' CV post-seg.']);
xticks([1,2,3,4]);
xticklabels({'Fluor focus \newline 1 slice from max', 'Fluor focus \newline 2slices from max', ...
    'Fluor focus \newline 1 slice from Ph3', 'Fluor focus \newline 2 slices from Ph3'})
ylabel([channels(fluor_chan) 'Difference in normalized CV ']);
grid on;
legend('location', 'best');


set(gcf,'units','points','position',[-1200,50,720*1.5,275*1.25*2]);

saveas(f, ['postsegmentation_fluor_analysis\phase3_maxvar_z_vs_' ...
    channels{fluor_chan} '_maxvar_z_maxcv_z_BAR.png']);
saveas(f, ['postsegmentation_fluor_analysis\phase3_maxvar_z_vs_' ...
    channels{fluor_chan} '_maxvar_z_maxcv_z_BAR.fig']);


% function padded_vector = pad_vector(input_vector, model_matrix)
%     padding_per_row = sum(model_matrix == 0, 2);
%     [rows,cols] = size(model_matrix);
%     padded_vector = [];
%     for r = 1:rows
%         padded_vector = [padded_vector, input_vector((r-1)*cols + 1:(cols-padding_per_row)), zeros(1,padding_per_row)];
%     end
% end
%}


%% Variance & Intensity & CV vs. z-slice POST (cell by cell comparison)
[num_sam, num_img] = size(ip.exp.splimg2frm);

var_lst_cv = cell(length(channels),1);
var_lst_var = cell(length(channels),1);
phase3_bestzPRE = [];

keeptrackofthesample = [];
for c = 1:length(channels)
    %z slice #
    spl2znum = ip.image.zrange .* ...
        ones(1,max(ip.exp.sampleList)) ;
    
    %work our way down for y axis limits
    y_min(1) = 1;
    y_min(2) = 1;
    y_min(3) = 1;
    
    var_lst_cv{c} = [];
    var_lst_var{c} = [];
    %loop frames
    for n_frame = 1:ip.exp.totalframes
        n_image = ip.exp.frm2img(n_frame);
        %progress updates
        if n_image == 1
            fprintf(1, newline);
            fprintf(1,['Sample ' num2str(ip.exp.frm2spl(n_frame)) ' of ' ...
                num2str(max(ip.exp.sampleList)) '.' newline]);
            progress_2 = [];
        end
        for d_=1:1:size(progress_2,2) ; fprintf(1,'\b') ; end
        progress_2 = [sprintf('\t') 'Image ' num2str(n_image) ' in ' ...
            num2str(min(ip.sample{ip.exp.frm2spl(n_frame)}.idx)) '-' ...
            num2str(max(ip.sample{ip.exp.frm2spl(n_frame)}.idx)) '.'] ;
        fprintf(1,progress_2) ;
        %load image
        load([ip.seg.dir ip.seg.base_name 'seg' num2str(n_frame,'%03d') '.mat'],'LcFull') ;
        imdigitnum_ = 2 ;
        current_spl = ip.exp.frm2spl(n_frame) ;
        if max(ip.sample{current_spl}.idx) <= 9
            imdigitnum_ = 1 ;
        end
        var_lst = zeros(1,ip.image.zrange) ;

        CellMap = LcFull;
        nnn=double(LcFull);
        N=max(nnn(:));
        
        cv_temp = zeros(spl2znum(current_spl),N);
        var_temp = zeros(spl2znum(current_spl),N);
        %loop z
        for iz = 1: spl2znum(current_spl)
        PH3 = imread([ip.exp.path 'images\' ...
                ip.image.base_name '_' ...
                num2str(ip.exp.frm2spl(n_frame),'%03d') ...
                'xy' num2str(ip.exp.frm2img(n_frame),['%0' num2str(imdigitnum_ ,'%01d') 'd']) ...
                ['z' num2str(iz ,'%01d') 'c' num2str(c) '.tif'] ]);

            %loop cells
            for index=1:N
                Cell = immultiply( CellMap == index, PH3);
                int_per_pix = Cell(CellMap == index);

                cv_temp(iz,index) = sqrt(var(double(int_per_pix(:)))) / mean(double(int_per_pix(:))) ;
                var_temp(iz,index) = var(double(int_per_pix(:))) ;
            end
        end
        
        % max for each cell
        [~, cv_temp] = max(cv_temp) ;
        var_lst_cv{c} = [var_lst_cv{c}, cv_temp];
        [~, var_temp] = max(var_temp) ;
        var_lst_var{c} = [var_lst_var{c}, var_temp];
        
%         cv_temp = zeros(spl2znum(current_spl),N);
%         var_temp = zeros(spl2znum(current_spl),N);
%         %go over the whole image only (rather than cell by cell)
%         if strcmp(channels{c}, 'phase3')
%             var_temp = zeros(spl2znum(current_spl),1);
%             for iz = 1: spl2znum(current_spl)
%                 PH3 = imread([ip.exp.path 'images\' ...
%                     ip.image.base_name '_' ...
%                     num2str(ip.exp.frm2spl(n_frame),'%03d') ...
%                     'xy' num2str(ip.exp.frm2img(n_frame),['%0' num2str(imdigitnum_ ,'%01d') 'd']) ...
%                     ['z' num2str(iz ,'%01d') 'c' num2str(c) '.tif'] ]);
%                 
% %                 cv_temp(iz) = sqrt(var(double(PH3(:)))) / mean(double(PH3(:))) ;
%                 var_temp(iz) = var(double(PH3(:))) ;
%             end
%             [~,var_temp] = max(var_temp);
%             phase3_bestzPRE = [phase3_bestzPRE, var_temp];
%         end
%         
        keeptrackofthesample = [keeptrackofthesample; ip.exp.frm2spl(n_frame)*ones(length(var_temp),1)];
    end
end
%
%plotting
dot_scale = 1;
for fluor_chan = 2:length(channels)
    %%% Best phase3 variance vs. best FLUOR variance
    cmap = hsv(max(ip.exp.frm2spl));
    f = figure ();
    subplot(2,2,1);
    for i = 1:max(ip.exp.frm2spl)
        x = var_lst_var{1}(keeptrackofthesample(1:length(keeptrackofthesample)/length(channels)) == i);
        y = var_lst_var{fluor_chan}(keeptrackofthesample(1:length(keeptrackofthesample)/length(channels)) == i);
        xy = [x' - max(ip.exp.frm2spl)/2 * .1 + i * .1,...
            y' - max(ip.exp.frm2spl)/2 * .1 + i * .1];
        [C,~,ic] = unique(xy, 'rows');
        a_counts = accumarray(ic,1);

        scatter(C(:,1), C(:,2), a_counts*dot_scale, 'markeredgecolor', 'none',...
            'markerfacecolor', cmap(i,:), 'DisplayName', ['Sample ', num2str(i)]);
        hold on;
    end
    plot(0:8, 0:8, 'Color', 'black', 'DisplayName', 'Slope = 1');
    xlim([0 8]);
    ylim([0 8]);
    title(['Z w/ max Ph3 var. cell vs. Z w/ max ' channels{fluor_chan} ' var. cell']);
    xlabel('Phase 3 highest variance z slice');
    ylabel([channels(fluor_chan) ' highest variance z slice']);
    grid on;
    legend('location', 'best');

    subplot(2,2,3);
    xy = [var_lst_var{1}(:), var_lst_var{fluor_chan}(:)];
    [C,~,ic] = unique(xy, 'rows');
    a_counts = accumarray(ic,1);
    scatter(C(:,1), C(:,2), a_counts*dot_scale, 'markeredgecolor', 'none',...
        'markerfacecolor', 'black', 'DisplayName', ['All samples together']);
    hold on;
    plot(0:8, 0:8, 'Color', 'black', 'DisplayName', 'Slope = 1');
    hold on;
    x = 1:ip.image.zrange;
    y = zeros(1,ip.image.zrange);
    for i = x
        y(i) = mean(var_lst_var{fluor_chan}(var_lst_var{1} == i));
        y_er(i) = std(var_lst_var{fluor_chan}(var_lst_var{1} == i));
    end
    errorbar(x,y,y_er, 'displayname', 'Mean fluor z slice +/- std'); hold on;
    xlim([0 8]);
    ylim([0 8]);
    title(['Z w/ max Ph3 var. cell vs. Z w/ max ' channels{fluor_chan} ' var. cell']);
    xlabel('Phase 3 highest variance z slice');
    ylabel([channels(fluor_chan) ' highest variance z slice']);
    grid on;
    legend('location', 'best');


    subplot(2,2,2);
    for i = 1:max(ip.exp.frm2spl)
        x = var_lst_cv{1}(keeptrackofthesample(1:length(keeptrackofthesample)/length(channels)) == i);
        y = var_lst_cv{fluor_chan}(keeptrackofthesample(1:length(keeptrackofthesample)/length(channels)) == i);
        xy = [x' - max(ip.exp.frm2spl)/2 * .1 + i * .1,...
            y' - max(ip.exp.frm2spl)/2 * .1 + i * .1];
        [C,~,ic] = unique(xy, 'rows');
        a_counts = accumarray(ic,1);

        scatter(C(:,1), C(:,2), a_counts*dot_scale, 'markeredgecolor', 'none',...
            'markerfacecolor', cmap(i,:), 'DisplayName', ['Sample ', num2str(i)]);
        hold on;
    end
    plot(0:8, 0:8, 'Color', 'black', 'DisplayName', 'Slope = 1');
    xlim([0 8]);
    ylim([0 8]);
    title(['Z w/ max Ph3 var cell vs. Z w/ max ' channels{fluor_chan} ' CV cell']);
    xlabel('Phase 3 highest CV z slice');
    ylabel([channels{fluor_chan} ' highest CV z slice']);
    grid on;
    legend('location', 'best');

    subplot(2,2,4);
    xy = [var_lst_var{1}(:), var_lst_cv{fluor_chan}(:)];
    [C,~,ic] = unique(xy, 'rows');
    a_counts = accumarray(ic,1);
    scatter(C(:,1), C(:,2), a_counts*dot_scale, 'markeredgecolor', 'none',...
        'markerfacecolor', 'black', 'DisplayName', ['All samples together']);
    hold on;
    plot(0:8, 0:8, 'Color', 'black', 'DisplayName', 'Slope = 1'); 
    hold on;
    x = 1:ip.image.zrange;
    y = zeros(1,ip.image.zrange);
    for i = x
        y(i) = mean(var_lst_cv{fluor_chan}(var_lst_cv{1} == i));
        y_er(i) = std(var_lst_cv{fluor_chan}(var_lst_cv{1} == i));
    end
    errorbar(x,y,y_er, 'displayname', 'Mean fluor z slice +/- std'); hold on;
    xlim([0 8]);
    ylim([0 8]);
    title(['Z w/ max Ph3 var cell vs. Z w/ max ' channels{fluor_chan} ' CV cell']);
    xlabel('Phase 3 highest CV z slice');
    ylabel([channels{fluor_chan} ' highest CV z slice']);
    grid on;
    legend('location', 'best');

    set(gcf,'units','points','position',[-1200,50,720*1.5,275*1.25*2]);

    saveas(f, ['postsegmentation_fluor_analysis\phase3_maxvar_z_vs_' ...
        channels{fluor_chan} '_maxvar_z_maxcv_z__cell_by_cell.png']);
    saveas(f, ['postsegmentation_fluor_analysis\phase3_maxvar_z_vs_' ...
        channels{fluor_chan} '_maxvar_z_maxcv_z__cell_by_cell.fig']);


%     %%% Best phase3 variance vs. best FLUOR variance bar of difference
%     pick_z_fluorvar = transpose(fluor_bestvarz(ip.exp.frm2spl ~= 8));
%     [all_x, all_y] = ind2sub(size(fluor_varz(:,:,1)), 1:numel(fluor_varz(:,:,1)));
%     all_x = all_x(fluor_varz(:,:,1) ~= 0);
%     all_y = all_y(fluor_varz(:,:,1) ~= 0);
%     best_var = fluor_varz(sub2ind(size(fluor_varz), all_x, all_y, pick_z_fluorvar));
% 
%     best_var_t = best_var(pick_z_fluorvar > 2 & pick_z_fluorvar < ip.image.zrange-1);
%     pick_z_fluorvar_t = pick_z_fluorvar(pick_z_fluorvar > 2 & pick_z_fluorvar < ip.image.zrange-1);
%     all_d1_var = [abs(best_var - fluor_varz(sub2ind(size(fluor_varz), all_x, all_y, pick_z_fluorvar+1))), ...
%         abs(best_var - fluor_varz(sub2ind(size(fluor_varz), all_x, all_y, pick_z_fluorvar-1)))];
%     all_d2_var = [abs(best_var_t - fluor_varz(sub2ind(size(fluor_varz), ...
%         all_x(pick_z_fluorvar > 2 & pick_z_fluorvar < ip.image.zrange-1), ...
%         all_y(pick_z_fluorvar > 2 & pick_z_fluorvar < ip.image.zrange-1), pick_z_fluorvar_t+2))), ...
%         abs(best_var_t - fluor_varz(sub2ind(size(fluor_varz), ...
%         all_x(pick_z_fluorvar > 2 & pick_z_fluorvar < ip.image.zrange-1), ...
%         all_y(pick_z_fluorvar > 2 & pick_z_fluorvar < ip.image.zrange-1), pick_z_fluorvar_t-2)))];
% 
%     pick_z_ph3 = transpose(phase3_bestzPRE(ip.exp.frm2spl ~= 8));
%     ph3_var = fluor_varz(sub2ind(size(fluor_varz), all_x, all_y, pick_z_ph3));
%     d_var = abs(best_var - ph3_var);
%     idx_var = abs(pick_z_fluorvar - pick_z_ph3);
%     %
%     pick_z_fluorcv = transpose(fluor_bestCVz(ip.exp.frm2spl ~= 8));
%     [all_x, all_y] = ind2sub(size(fluor_CVz(:,:,1)), 1:numel(fluor_CVz(:,:,1)));
%     all_x = all_x(fluor_CVz(:,:,1) ~= 0);
%     all_y = all_y(fluor_CVz(:,:,1) ~=0);
%     best_CV = fluor_CVz(sub2ind(size(fluor_CVz), all_x, all_y, pick_z_fluorcv));
% 
%     best_CV_t = best_CV(pick_z_fluorcv > 2 & pick_z_fluorcv < ip.image.zrange-1);
%     pick_z_fluorcv_t = pick_z_fluorcv(pick_z_fluorcv > 2 & pick_z_fluorcv < ip.image.zrange-1);
%     all_d1_CV = [abs(best_CV(pick_z_fluorcv > 1 & pick_z_fluorcv < ip.image.zrange) - fluor_CVz(sub2ind(size(fluor_CVz), ...
%         all_x(pick_z_fluorcv > 1 & pick_z_fluorcv < ip.image.zrange), ...
%         all_y(pick_z_fluorcv > 1 & pick_z_fluorcv < ip.image.zrange), ...
%         pick_z_fluorcv(pick_z_fluorcv > 1 & pick_z_fluorcv < ip.image.zrange)+1))), ...
%         abs(best_CV(pick_z_fluorcv > 1 & pick_z_fluorcv < ip.image.zrange) - fluor_CVz(sub2ind(size(fluor_CVz), ...
%         all_x(pick_z_fluorcv > 1 & pick_z_fluorcv < ip.image.zrange), ...
%         all_y(pick_z_fluorcv > 1 & pick_z_fluorcv < ip.image.zrange), ...
%         pick_z_fluorcv(pick_z_fluorcv > 1 & pick_z_fluorcv < ip.image.zrange)-1)))];
%     all_d2_CV = [abs(best_CV_t - fluor_CVz(sub2ind(size(fluor_CVz), ...
%         all_x(pick_z_fluorcv > 2 & pick_z_fluorcv < ip.image.zrange-1), ...
%         all_y(pick_z_fluorcv > 2 & pick_z_fluorcv < ip.image.zrange-1), pick_z_fluorcv_t+2))), ...
%         abs(best_CV_t - fluor_CVz(sub2ind(size(fluor_CVz), ...
%         all_x(pick_z_fluorcv > 2 & pick_z_fluorcv < ip.image.zrange-1), ...
%         all_y(pick_z_fluorcv > 2 & pick_z_fluorcv < ip.image.zrange-1), pick_z_fluorcv_t-2)))];
% 
%     ph3_CV = fluor_CVz(sub2ind(size(fluor_CVz), all_x, all_y, pick_z_ph3));
%     d_CV = abs(best_CV - ph3_CV);
%     idx_cv = abs(pick_z_fluorcv - pick_z_ph3);
%     %
%     f = figure ();
%     subplot(1,2,1);
%     bar(1, mean(all_d1_var), 'DisplayName', 'Mean \Deltavar. 1 zstep from max', 'facecolor', cmap(1,:));
%     hold on;
%     errorbar(1, mean(all_d1_var), std(all_d1_var), 'DisplayName', 'Std of all (+/- 1 from best Z)', 'color', cmap(1,:));
%     hold on;
%     swarmchart(1*ones(size(all_d1_var)), all_d1_var, 'HandleVisibility', 'off', 'MarkerFaceColor', cmap(1,:),'MarkerEdgeAlpha',0.4);
%     hold on;
%     %
%     bar(2, mean(all_d2_var), 'DisplayName', 'Mean \Deltavar. 2 zstep from max', 'facecolor', cmap(2,:));
%     hold on;
%     errorbar(2, mean(all_d2_var), std(all_d2_var), 'DisplayName', 'Std of all (+/- 2 from best Z)', 'color', cmap(2,:));
%     hold on;
%     swarmchart(2*ones(size(all_d2_var)), all_d2_var, 'HandleVisibility', 'off', 'MarkerFaceColor', cmap(2,:),'MarkerEdgeAlpha',0.4);
%     hold on;
%     %
%     bar(3, mean(d_var(idx_var == 1)), 'DisplayName', 'Mean \Deltavar. in samp. w/ 1 zstep diff.', 'facecolor', cmap(3,:));
%     hold on;
%     errorbar(3, mean(d_var(idx_var == 1)), std(d_var(idx_var == 1)), 'DisplayName', 'Std of all +/- 1 z diff', 'color', cmap(3,:));
%     hold on;
%     swarmchart(3*ones(size(d_var(idx_var == 1))), d_var(idx_var == 1), 'HandleVisibility', 'off', 'MarkerFaceColor', cmap(3,:),'MarkerEdgeAlpha',0.4);
%     hold on;
%     %
%     bar(4, mean(d_var(idx_var == 2)), 'DisplayName', 'Mean \Deltavar. in samp. w/ 2 zstep diff', 'facecolor', cmap(4,:));
%     hold on;
%     errorbar(4, mean(d_var(idx_var == 2)), std(d_var(idx_var == 2)), 'DisplayName', 'Std of all +/- 2 z diff', 'color', cmap(4,:));
%     hold on;
%     swarmchart(4*ones(size(d_var(idx_var == 2))), d_var(idx_var == 2), 'HandleVisibility', 'off', 'MarkerFaceColor', cmap(4,:),'MarkerEdgeAlpha',0.4);
%     hold on;
% 
%     xlim([.5 4.5]);
%     ylim([0 .6]);
%     title(['Z w/ max Ph3 var. pre-seg. vs. Z w/ max ' channels{fluor_chan} ' var. post-seg.']);
%     xticks([1,2,3,4]);
%     xticklabels({'Fluor focus \newline 1 slice from max', 'Fluor focus \newline 2slices from max', ...
%         'Fluor focus \newline 1 slice from Ph3', 'Fluor focus \newline 2 slices from Ph3'})
%     ylabel([channels(fluor_chan) ' difference in normalized variance']);
%     grid on;
%     legend('location', 'best');
% 
% 
%     subplot(1,2,2);
%     bar(1, mean(all_d1_CV), 'DisplayName', 'Mean \Deltacv 1 zstep from max', 'facecolor', cmap(1,:));
%     hold on;
%     errorbar(1, mean(all_d1_CV), std(all_d1_CV), 'DisplayName', 'Std of all (+/- 1 from best Z)', 'color', cmap(1,:));
%     hold on;
%     swarmchart(1*ones(size(all_d1_CV)), all_d1_CV, 'HandleVisibility', 'off', 'MarkerFaceColor', cmap(1,:),'MarkerEdgeAlpha',0.4);
%     hold on;
%     %
%     bar(2, mean(all_d2_CV), 'DisplayName', 'Mean \Deltacv 2 zstep from max', 'facecolor', cmap(2,:));
%     hold on;
%     errorbar(2, mean(all_d2_CV), std(all_d2_CV), 'DisplayName', 'Std of all (+/- 2 from best Z)', 'color', cmap(2,:));
%     hold on;
%     swarmchart(2*ones(size(all_d2_CV)), all_d2_CV, 'HandleVisibility', 'off', 'MarkerFaceColor', cmap(2,:),'MarkerEdgeAlpha',0.4);
%     hold on;
%     %
%     bar(3, mean(d_CV(idx_cv == 1)), 'DisplayName', 'Mean \Deltacv in samp. w/ 1 zstep diff.', 'facecolor', cmap(3,:));
%     hold on;
%     errorbar(3, mean(d_CV(idx_cv == 1)), std(d_CV(idx_cv == 1)), 'DisplayName', 'Std of all +/- 1 z diff', 'color', cmap(3,:));
%     hold on;
%     swarmchart(3*ones(size(d_CV(idx_cv == 1))), d_CV(idx_cv == 1), 'HandleVisibility', 'off', 'MarkerFaceColor', cmap(3,:),'MarkerEdgeAlpha',0.4);
%     hold on;
%     %
%     bar(4, mean(d_CV(idx_cv == 2)), 'DisplayName', 'Mean \Deltacv in samp. w/ 2 zstep diff', 'facecolor', cmap(4,:));
%     hold on;
%     errorbar(4, mean(d_CV(idx_cv == 2)), std(d_CV(idx_cv == 2)), 'DisplayName', 'Std of all +/- 2 z diff', 'color', cmap(4,:));
%     hold on;
%     swarmchart(4*ones(size(d_CV(idx_cv == 2))), d_CV(idx_cv == 2), 'HandleVisibility', 'off', 'MarkerFaceColor', cmap(4,:),'MarkerEdgeAlpha',0.4);
%     hold on;
% 
%     xlim([.5 4.5]);
%     ylim([0 .6]);
%     title(['Z w/ max Ph3 CV pre-seg. vs. Z w/ max ' channels{fluor_chan} ' CV post-seg.']);
%     xticks([1,2,3,4]);
%     xticklabels({'Fluor focus \newline 1 slice from max', 'Fluor focus \newline 2slices from max', ...
%         'Fluor focus \newline 1 slice from Ph3', 'Fluor focus \newline 2 slices from Ph3'})
%     ylabel([channels(fluor_chan) 'Difference in normalized CV ']);
%     grid on;
%     legend('location', 'best');
% 
% 
%     set(gcf,'units','points','position',[-1200,50,720*1.5,275*1.25*2]);
% 
%     saveas(f, ['postsegmentation_fluor_analysis\phase3_maxvar_z_vs_' ...
%         channels{fluor_chan} '_maxvar_z_maxcv_z_BAR.png']);
%     saveas(f, ['postsegmentation_fluor_analysis\phase3_maxvar_z_vs_' ...
%         channels{fluor_chan} '_maxvar_z_maxcv_z_BAR.fig']);
end

>>>>>>> 88dd4f4ac7eb6f5ae758552b04c47df241263be0
