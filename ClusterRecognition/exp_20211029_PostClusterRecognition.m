<<<<<<< HEAD
% Initialize variables

clc
clear

ip = exp_20211029_InitializeExp();

spl2cnum = struct() ; 

spl2cnum.GFP = 4*ones(1,max(ip.exp.sampleList)) ;
spl2cnum.DAPI = 5*ones(1,max(ip.exp.sampleList)) ;

whichchanel = 'GFP' ;  

whichspl2cnum = spl2cnum.(whichchanel) ; 

% datafolder = [pwd '\cluster_quantify\Run-TXRED-.4low.7high-cellnormed-thresh-02-Nov-2021'] ;

whichsamples = ip.exp.sampleList ;

plot_col = 5 ;
plot_row = 2 ;

% if exist([datafolder '\data\' 'spots_data.mat'],'file')
%     load([datafolder '\data\' 'spots_data.mat'],'clusterlistlow_new');
% end
% if exist([datafolder '\data\' 'spots_data.mat'],'file')
%     load([datafolder '\data\' 'spots_data.mat'],'clusterlisthigh_new');
% end 
short_name = {'1mM IPTG', '10 \muM IPTG', '0mM IPTG'};

%% group clusterdata_X.mat of individual frames
for combo = 3%1:5
    
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
    
    datafolder = [ip.exp.path '\cluster_quantify\Run-GFP-' num2str(low_threshold) ...
        'low' num2str(high_threshold) 'high-cellnormed-thresh-07-Dec-2021\'];

whichsamples = ip.exp.sampleList ;

plot_col = 5 ;
plot_row = 3 ;

if exist([datafolder '\data\' 'spots_data.mat'],'file')
    load([datafolder '\data\' 'spots_data.mat'],'clusterlistlow_new');
end
if exist([datafolder '\data\' 'spots_data.mat'],'file')
    load([datafolder '\data\' 'spots_data.mat'],'clusterlisthigh_new');
end 



% only run for frames I segmented.
%ip.sample{2}.idx = 1:8 ;

clusterlistlow_new = cell(1,max(ip.exp.sampleList)) ;
clusterlisthigh_new = cell(1,max(ip.exp.sampleList)) ;

for n_sample = whichsamples
    
    fprintf(1,['Group ' num2str(n_sample) ' of ' num2str(max(ip.exp.sampleList)) '.' sprintf('\n')]);
    
    spots_in_group1 = [] ;
    spots_in_group2 = [] ;
    
    r1 = 0 ;
    r2 = 0 ;
    
    progress_2 = [];
    for n_image = ip.sample{n_sample}.idx

        n_frame = ip.exp.splimg2frm(n_sample,n_image);
        
        % get the channel number from the sample numebr
        SpotChannel = whichspl2cnum(n_sample) ; 
        sr = InitializeSpotRecognitionParameters_c234(ip,n_image,SpotChannel,datafolder); %% update the file name here. 
        
        load([sr.output 'clusterdata_low' num2str(n_frame,'%03d') '.mat']);
        
        if ~isempty(clusterdata_low)
            spots_in_group1 = [spots_in_group1 ; ...
                clusterdata_low(:,2) ...            % Z-slice  (1)
                [1:size(clusterdata_low,1)]' ...    % vector with a list of spot numbers for the frame (2)
                clusterdata_low(:,[3 4 5]) ...      % SpotInt area AvgSpotInt(3 4 5)
                clusterdata_low(:,1) ...    % cellnum (6)
                clusterdata_low(:,[6]) ... % frame (7)
                clusterdata_low(:,[9]) ... % low cluster number (8)
                clusterdata_low(:,[7 8])] ; % average x y position (9 10)
        end
        
        load([sr.output 'clusterdata_high' num2str(n_frame,'%03d') '.mat'],'clusterdata_high');
        
        if ~isempty(clusterdata_high)
            spots_in_group2 = [spots_in_group2 ; ...
                clusterdata_high(:,2) ...         % Z-slice  (1)
                [1:size(clusterdata_high,1)]' ...  % vector with a list of spot numbers for the frame (2)
                clusterdata_high(:,[3 4 5]) ...   % SpotInt area AvgSpotInt(3 4 5)
                clusterdata_high(:,1) ...  % cellnum (6)
                clusterdata_high(:, 6) ... % frame (7)
                clusterdata_high(:, 9) ... % low cluster number (8)
                clusterdata_high(:,[7 8])] ; % average x y position (9 10)
        end
    end
    
    if ~isempty(spots_in_group1)
        clusterlistlow_new{n_sample} = spots_in_group1 ;
    else
        clusterlistlow_new{n_sample} = zeros(1,7) ;
    end
    
    if ~isempty(spots_in_group2)
        clusterlisthigh_new{n_sample} = spots_in_group2 ;
    else
        clusterlisthigh_new{n_sample} = zeros(1,8) ;
    end
    
end

if exist([datafolder '\data\' 'spots_data.mat'],'file')
    save([datafolder '\data\' 'spots_data.mat'],'clusterlistlow_new','-append');
else
    save([datafolder '\data\' 'spots_data.mat'],'clusterlistlow_new');
end
if exist([datafolder '\data\' 'spots_data.mat'],'file')
    save([datafolder '\data\' 'spots_data.mat'],'clusterlisthigh_new','-append');
else
    save([datafolder '\data\' 'spots_data.mat'],'clusterlisthigh_new');
end
end

%% 1D histogram of spot int/pix
figure('Units','normalized',...
    'color','w',...
    'OuterPosition',[0.0    0.045    0.3    .6],...
    'NumberTitle','off','Name','1D spot intensity') ;
ha = tight_subplot(3, ceil(length(whichsamples)/3), [.02 .05], [.08 .01], [.05 .01]);
    
for combo = 1:5
    
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
    
    datafolder = [ip.exp.path '\cluster_quantify\Run-GFP-' num2str(low_threshold) ...
        'low' num2str(high_threshold) 'high-cellnormed-thresh-16-Nov-2021\'];
    
    if exist([datafolder '\data\' 'spots_data.mat'],'file')
        load([datafolder '\data\' 'spots_data.mat'],'clusterlistlow_new');
    end
    if exist([datafolder '\data\' 'spots_data.mat'],'file')
        load([datafolder '\data\' 'spots_data.mat'],'clusterlisthigh_new');
    end 
    
%     figure('Units','normalized',...
%         'color','w',...
%         'OuterPosition',[0.0    0.045    0.3    .6],...
%         'NumberTitle','off','Name','1D spot intensity') ;
    binN = 100 ;
%     ha = tight_subplot(3, ceil(length(whichsamples)/3), [.02 .05], [.08 .01], [.05 .01]);

    iplot = 1 ;
    for isample = whichsamples

        if ~isempty(clusterlistlow_new{isample})

            axes(ha(iplot)) ;
            set(gca, 'xscale','log','fontsize',6)

            box on ; hold all ;

            minP = min(clusterlistlow_new{isample}(:,3) ./ clusterlistlow_new{isample}(:,4)) ;
            maxP = max(clusterlistlow_new{isample}(:,3) ./ clusterlistlow_new{isample}(:,4)) ;
            nbin = logspace(log10(minP),log10(maxP),binN) ;
            %
            [ya,xa] = hist(clusterlistlow_new{isample}(:,3) ./ clusterlistlow_new{isample}(:,4),nbin) ;
            %
            plot(xa,ya/sum(ya),'-',...
                'Color',[1 0 0] * (combo + 5 - 2) / (2*5 - 2),...
                'LineWidth',1,...
                'DisplayName',[short_name{isample} '(low threshold:' num2str(low_threshold) ')']) ;
            %%%
            y = clusterlisthigh_new{isample}(:,3) ./ clusterlisthigh_new{isample}(:,4);
            minP = min(y(~isnan(y))) ;
            maxP = max(y(~isnan(y))) ;
            if isempty(minP) || isempty(maxP)
                iplot = iplot + 1 ;
                continue;
            end
            nbin = logspace(log10(minP),log10(maxP),binN) ;
            %
            [ya,xa] = hist(y(~isnan(y)),nbin) ;
            %
            plot(xa,ya/sum(ya),'-',...
                'Color',[0 0 1] * (combo + 5 - 2) / (2*5 - 2),...
                'LineWidth',1,...
                'DisplayName',[short_name{isample} '(high threshold:' num2str(high_threshold) ')']) ;

            lh=legend('show');
            set(lh,'location','northwest');
            lh.FontSize = 8;

            xlim([.7e3 2e3]) ;        
    %         ylabel('Probability') ;
            set(gca,'XTickLabel',[]) ;
            iplot = iplot + 1 ;
%             set(gca,'XTickLabelMode','auto') ;
%             xlabel('Spot intensity/pixel') ;
        end
    end
    set(gca,'XTickLabelMode','auto') ;
    xlabel('Spot intensity/pixel') ;
    
%     savefig(gcf, [datafolder '\fit\1D_spot_intensity.fig']);
%     saveas(gcf, [datafolder '\fit\1D_spot_intensity.png']);
end
    savefig(gcf, [pwd '\cluster_rec\1D_spot_intensity.fig']);
    saveas(gcf, [pwd '\cluster_rec\1D_spot_intensity.png']);

%% 1D histogram of spot area
%{
figure('Units','normalized',...
    'color','w',...
    'OuterPosition',[0.0    0.045    0.3*2    .95],...
    'NumberTitle','off','Name','1D spot intensity') ;
binN = 100 ;
ha = tight_subplot(ceil(length(whichsamples)/2), 2, [.02 .05], [.05 .01], [.05 .01]);

iplot = 1 ;
for isample = whichsamples
    
    if ~isempty(clusterlistlow_new{isample})
        
        axes(ha(iplot)) ;
        set(gca, 'xscale','log','fontsize',6)

%         subplot(length(whichsamples),1,iplot,...
%             'XScale','log','FontSize',8) ;
        
        box on ; hold all ;
        
        minP = min(clusterlistlow_new{isample}(:,4)) ;
        maxP = max(clusterlistlow_new{isample}(:,4)) ;
        nbin = logspace(log10(minP),log10(maxP),binN) ;
        
        [ya,xa] = hist(clusterlistlow_new{isample}(:,4),nbin) ;
        
        plot(xa,ya/sum(ya),'-',...
            'Color','r',...
            'LineWidth',1,...
            'DisplayName',[ip.sample{isample}.name '(low threshold)']) ;
        %%%
        minP = min(clusterlisthigh_new{isample}(:,4)) ;
        maxP = max(clusterlisthigh_new{isample}(:,4)) ;
        nbin = logspace(log10(minP),log10(maxP),binN) ;
        
        [ya,xa] = hist(clusterlisthigh_new{isample}(:,4),nbin) ;
        
        plot(xa,ya/sum(ya),'-',...
            'Color','b',...
            'LineWidth',1,...
            'DisplayName',[ip.sample{isample}.name '(high threshold)']) ;
        %
        lh=legend('show');
        set(lh,'location','northeast');
        lh.FontSize = 8;

        xlim([1e1 1e4]) ;        
%         ylabel('Probability') ;
        set(gca,'XTickLabel',[]) ;
        iplot = iplot + 1 ;
        set(gca,'XTickLabelMode','auto') ;
        xlabel('Spot area') ;
    end
end
set(gca,'XTickLabelMode','auto') ;
xlabel('Spot area') ;
%}
figure('Units','normalized',...
    'color','w',...
    'OuterPosition',[0.0    0.045    0.3    .6],...
    'NumberTitle','off','Name','1D spot intensity') ;
ha = tight_subplot(3, ceil(length(whichsamples)/3), [.02 .05], [.08 .01], [.05 .01]);
    
for combo = 1:5
    
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
    
    datafolder = [ip.exp.path '\cluster_quantify\Run-GFP-' num2str(low_threshold) ...
        'low' num2str(high_threshold) 'high-cellnormed-thresh-16-Nov-2021\'];
    
    if exist([datafolder '\data\' 'spots_data.mat'],'file')
        load([datafolder '\data\' 'spots_data.mat'],'clusterlistlow_new');
    end
    if exist([datafolder '\data\' 'spots_data.mat'],'file')
        load([datafolder '\data\' 'spots_data.mat'],'clusterlisthigh_new');
    end 
    
%     figure('Units','normalized',...
%         'color','w',...
%         'OuterPosition',[0.0    0.045    0.3    .6],...
%         'NumberTitle','off','Name','1D spot intensity') ;
    binN = 100 ;
%     ha = tight_subplot(3, ceil(length(whichsamples)/3), [.02 .05], [.08 .01], [.05 .01]);

    iplot = 1 ;
    for isample = whichsamples

        if ~isempty(clusterlistlow_new{isample})

            axes(ha(iplot)) ;
            set(gca, 'xscale','log','fontsize',6)

            box on ; hold all ;

            minP = min(clusterlistlow_new{isample}(:,4)) ;
            maxP = max(clusterlistlow_new{isample}(:,4)) ;
            nbin = logspace(log10(minP),log10(maxP),binN) ;
            %
            [ya,xa] = hist(clusterlistlow_new{isample}(:,4),nbin) ;
            %
            plot(xa,ya/sum(ya),'-',...
                'Color',[1 0 0] * (combo + 5 - 2) / (2*5 - 2),...
                'LineWidth',1,...
                'DisplayName',[short_name{isample} '(low threshold:' num2str(low_threshold) ')']) ;
            %%%
            y = clusterlisthigh_new{isample}(:,4);
            minP = min(max(y(~isnan(y)),1)) ;
            maxP = max(y(~isnan(y))) ;
            if isempty(minP) || isempty(maxP)
                iplot = iplot + 1 ;
                continue;
            end
            nbin = logspace(log10(minP),log10(maxP),binN) ;
            %
            [ya,xa] = hist(y(~isnan(y)),nbin) ;
            %
            plot(xa,ya/sum(ya),'-',...
                'Color',[0 0 1] * (combo + 5 - 2) / (2*5 - 2),...
                'LineWidth',1,...
                'DisplayName',[short_name{isample} '(high threshold:' num2str(high_threshold) ')']) ;

            lh=legend('show');
            set(lh,'location','northwest');
            lh.FontSize = 8;

            xlim([1e0 5e3]) ; 
            ylim([0 .1])
    %         ylabel('Probability') ;
            set(gca,'XTickLabel',[]) ;
            iplot = iplot + 1 ;
%             set(gca,'XTickLabelMode','auto') ;
%             xlabel('Spot intensity/pixel') ;
        end
    end
    set(gca,'XTickLabelMode','auto') ;
    xlabel('Spot pixel area') ;
    
%     savefig(gcf, [datafolder '\fit\1D_spot_intensity.fig']);
%     saveas(gcf, [datafolder '\fit\1D_spot_intensity.png']);
end
    savefig(gcf, [pwd '\cluster_rec\1D_spot_area.fig']);
    saveas(gcf, [pwd '\cluster_rec\1D_spot_area.png']);
    
%% 1D histogram of spot int
%{
figure('Units','normalized',...
    'color','w',...
    'OuterPosition',[0.0    0.045    0.3*2    .95],...
    'NumberTitle','off','Name','1D spot intensity') ;
binN = 100 ;
ha = tight_subplot(ceil(length(whichsamples)/2), 2, [.02 .05], [.05 .01], [.05 .01]);

iplot = 1 ;
for isample = whichsamples
    
    if ~isempty(clusterlistlow_new{isample})
        
        axes(ha(iplot)) ;
        set(gca, 'xscale','log','fontsize',6)

        
        box on ; hold all ;
        
        minP = min(clusterlistlow_new{isample}(:,3)) ;
        maxP = max(clusterlistlow_new{isample}(:,3)) ;
        nbin = logspace(log10(minP),log10(maxP),binN) ;
        
        [ya,xa] = hist(clusterlistlow_new{isample}(:,3),nbin) ;
        
        plot(xa,ya/sum(ya),'-',...
            'Color','r',...
            'LineWidth',1,...
            'DisplayName',[ip.sample{isample}.name '(low threshold)']) ;
        %%%
        minP = min(clusterlisthigh_new{isample}(:,3)) ;
        maxP = max(clusterlisthigh_new{isample}(:,3)) ;
        nbin = logspace(log10(minP),log10(maxP),binN) ;
        
        [ya,xa] = hist(clusterlisthigh_new{isample}(:,3),nbin) ;
        
        plot(xa,ya/sum(ya),'-',...
            'Color','b',...
            'LineWidth',1,...
            'DisplayName',[ip.sample{isample}.name '(high threshold)']) ;
        
%         [ya,xa] = hist(spotlist_new{isample}(~inside_,5),nbin) ;
% %         
%         plot(xa,ya/sum(ya),'-',...
%             'Color',0.7*[1 1 1],...
%             'LineWidth',1,...
%             'DisplayName',[ip.sample{isample}.name '(outside)']) ;
        lh=legend('show');
        set(lh,'location','northeast');
        lh.FontSize = 8;

        xlim([1e4 1e7]) ;        
%         ylabel('Probability') ;
        set(gca,'XTickLabel',[]) ;
        iplot = iplot + 1 ;
        set(gca,'XTickLabelMode','auto') ;
        xlabel('Spot intensity') ;
    end
end
set(gca,'XTickLabelMode','auto') ;
xlabel('Spot intensity') ;
%}
figure('Units','normalized',...
    'color','w',...
    'OuterPosition',[0.0    0.045    0.3    .6],...
    'NumberTitle','off','Name','1D spot intensity') ;
ha = tight_subplot(3, ceil(length(whichsamples)/3), [.02 .05], [.08 .01], [.05 .01]);
    
for combo = 1:5
    
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
    
    datafolder = [ip.exp.path '\cluster_quantify\Run-GFP-' num2str(low_threshold) ...
        'low' num2str(high_threshold) 'high-cellnormed-thresh-16-Nov-2021\'];
    
    if exist([datafolder '\data\' 'spots_data.mat'],'file')
        load([datafolder '\data\' 'spots_data.mat'],'clusterlistlow_new');
    end
    if exist([datafolder '\data\' 'spots_data.mat'],'file')
        load([datafolder '\data\' 'spots_data.mat'],'clusterlisthigh_new');
    end 
    
%     figure('Units','normalized',...
%         'color','w',...
%         'OuterPosition',[0.0    0.045    0.3    .6],...
%         'NumberTitle','off','Name','1D spot intensity') ;
    binN = 100 ;
%     ha = tight_subplot(3, ceil(length(whichsamples)/3), [.02 .05], [.08 .01], [.05 .01]);

    iplot = 1 ;
    for isample = whichsamples

        if ~isempty(clusterlistlow_new{isample})

            axes(ha(iplot)) ;
            set(gca, 'xscale','log','fontsize',6)

            box on ; hold all ;

            minP = min(clusterlistlow_new{isample}(:,3));
            maxP = max(clusterlistlow_new{isample}(:,3));
            nbin = logspace(log10(minP),log10(maxP),binN) ;
            %
            [ya,xa] = hist(clusterlistlow_new{isample}(:,3),nbin) ;
            %
            plot(xa,ya/sum(ya),'-',...
                'Color',[1 0 0] * (combo + 5 - 2) / (2*5 - 2),...
                'LineWidth',1,...
                'DisplayName',[short_name{isample} '(low threshold:' num2str(low_threshold) ')']) ;
            %%%
            y = clusterlisthigh_new{isample}(:,3);
            minP = min(max(y(~isnan(y)),1)) ;
            maxP = max(y(~isnan(y))) ;
            if isempty(minP) || isempty(maxP)
                iplot = iplot + 1 ;
                continue;
            end
            nbin = logspace(log10(minP),log10(maxP),binN) ;
            %
            [ya,xa] = hist(y(~isnan(y)),nbin) ;
            %
            plot(xa,ya/sum(ya),'-',...
                'Color',[0 0 1] * (combo + 5 - 2) / (2*5 - 2),...
                'LineWidth',1,...
                'DisplayName',[short_name{isample} '(high threshold:' num2str(high_threshold) ')']) ;

            lh=legend('show');
            set(lh,'location','northwest');
            lh.FontSize = 8;

            xlim([1e3 1e7]) ;        
    %         ylabel('Probability') ;
            set(gca,'XTickLabel',[]) ;
            iplot = iplot + 1 ;
%             set(gca,'XTickLabelMode','auto') ;
%             xlabel('Spot intensity/pixel') ;
        end
    end
    set(gca,'XTickLabelMode','auto') ;
    xlabel('Spot intensity') ;
    
%     savefig(gcf, [datafolder '\fit\1D_spot_intensity.fig']);
%     saveas(gcf, [datafolder '\fit\1D_spot_intensity.png']);
end
    savefig(gcf, [pwd '\cluster_rec\1D_spot_intensity.fig']);
    saveas(gcf, [pwd '\cluster_rec\1D_spot_intensity.png']);
    
%% 2D plot only spot area vs peak height (not changed)
%{
markersize = 2;
colors = [[255 0 0 ];[200 200 200];[99 184 255];]/255;

figure('Units','normalized',...
    'color','w',...
    'OuterPosition',[ 0 0 1 1],...
    'NumberTitle','off','Name','spot area VS peakheight scatter') ;
% figure;
iplot = 1 ;
for isample = whichsamples
    
    % spotlist fields
    % 1  Z-slice
    % 2  spot number in the frame
    % 3  BG
    % 4  peakInt
    % 5  SpotInt
    % 6  pi*major*minor (area of spot)
    % 7  X
    % 8  Y
    % 9  cellnum
    % 10 pos_long_axis (normalized)
    % 11 pos_across (normalized)
    % 12 pos_long_axis_pix
    % 13 pos_across_axis_pix
    % 14 frame
    % 15 Major_axis
    % 16 Minor_axis
    % 17 imaginary spot flag
    
    % 3. spot area vs peak height
    subplot(plot_row , plot_col , iplot); hold on;
    
    inside_ = spotlist_new{isample}(:,9) ~= 0 ;
    plot(spotlist_new{isample}(~inside_,6),spotlist_new{isample}(~inside_,4), ...
        '.','color',0.7*[1 1 1],'markersize',markersize);
    hold on ;
    plot(spotlist_new{isample}(inside_,6),spotlist_new{isample}(inside_,4), ...
        '.','color','r','markersize',markersize);
    
%     if ~strcmp(whichchanel,'CFP')
%         plot([1e0 1e4] , [1 1]*NegPeakThreshold{1},'r--')
%     end
    
    xlabel('Spot area');
    ylabel('Peak height');
    ylim([1e1 1e6]);
    xlim([1e0 1e4]);
    set(gca , 'ytick', [1e1 1e2 1e3 1e4 1e5]) ;
    set(gca , 'xtick', [1e0 1e1 1e2 1e3 1e4]) ;
    xl = get(gca, 'xlim') ;
    
    hold off
    set(gca,'XScale','log','YScale','log')
    title(ip.sample{isample}.name) ;
%     legend([num2str(length(spotlist_new{isample}(:,6))) 'spots'],...
%         'location','northeast')
%     legend('outside','inside') ;
    grid on
    grid minor
    
    iplot = iplot + 1 ;
end
%}

%% heatmap to show the histograms (not changed)

% close all

markersize = 2;
colors = [[255 0 0 ];[200 200 200];[99 184 255];]/255;

figure('Units','normalized',...
    'color','w',...
    'OuterPosition',[0 0 1 1],...
    'NumberTitle','off','Name','spot area VS peakheight heatmap') ;

iplot = 1 ;
for isample = whichsamples
    
    % spotlist fields
    % 1  Z-slice
    % 2  spot number in the frame
    % 3  BG
    % 4  peakInt
    % 5  SpotInt
    % 6  pi*major*minor (area of spot)
    % 7  X
    % 8  Y
    % 9  cellnum
    % 10 pos_long_axis (normalized)
    % 11 pos_across (normalized)
    % 12 pos_long_axis_pix
    % 13 pos_across_axis_pix
    % 14 frame
    % 15 Major_axis
    % 16 Minor_axis
    % 17 imaginary spot flag
    
    
    % 3. spot area vs peak height
    subplot(plot_row , plot_col , iplot); hold on;
    xedges = logspace(0,4,100);
    yedges = logspace(1,5,100);
    edges = {xedges,yedges};
    H = hist3(spotlist_new{isample}(:,[6 4]),'Edges',edges);
    imagesc(H'); colormap jet;
    
    xlabel('Spot area');
    ylabel('Peak height');
    set(gca,'XTick',[0 25 50 75 100],'XTickLabel',{'0','1','2','3','4'},...
        'YTick',[0 25 50 75 100],'YTickLabel',{'1','2','3','4','5'});
    xlabel('log_{10}( Spot area )');
    ylabel('log_{10}( Peak height )');
    title(ip.sample{isample}.name) ;
    xlim([0 100]);
    ylim([0 100]);
    iplot = iplot + 1 ;
end

%% % 1. Re-group cluster data into cell data
for combo = 3%1:5
    
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
    
    datafolder = [ip.exp.path '\cluster_quantify\Run-GFP-' num2str(low_threshold) ...
        'low' num2str(high_threshold) 'high-cellnormed-thresh-07-Dec-2021\'];
    
enlistcell_new = cell(1,max(ip.exp.sampleList)) ;

for n_sample = 1:3

    fprintf(1,['Sample ' num2str(n_sample) ' of ' num2str(max(ip.exp.sampleList)) '.' sprintf('\n')]);
    
    cells_in_group = [] ;
    
    progress_2 = [];
    progress_3 = [];
    progress_4 = [];
    
    for n_image = ip.sample{n_sample}.idx
        
        n_frame = ip.exp.splimg2frm(n_sample,n_image);
        
        for d_=1:1:size(progress_2,2)+size(progress_3,2)+size(progress_4,2) ; fprintf(1,'\b') ; end
        progress_2 = [sprintf('\t') 'Image ' num2str(n_image) ' in ' ...
            num2str(min(ip.sample{n_sample}.idx)) '-' num2str(max(ip.sample{n_sample}.idx)) '.'] ;
        fprintf(1,progress_2) ;
        
        % determine the channel number
        DAPIChannel = spl2cnum.DAPI(n_sample);
        SpotChannel = spl2cnum.(whichchanel)(n_sample) ;
        sr = InitializeSpotRecognitionParameters_c234(ip,n_frame,SpotChannel,datafolder);
        DAPI_sr = InitializeSpotRecognitionParameters_c234(ip,n_frame,DAPIChannel,datafolder);
        % already in the gated data folder (no need for thresh or gate)
        load([datafolder '\data\' ...
            'clusterdata_low' num2str(n_frame,'%03d') '.mat'],'clusterdata_low');
        load([datafolder '\data\' ...
            'clusterdata_high' num2str(n_frame,'%03d') '.mat'],'clusterdata_high');
        
        load([sr.seg.dir sr.seg.name],'LcFull') ;

%         cell_DAPI = imread([sr.image.fullname 'z' ...
%             num2str(median(sr.image.zrange),'%01d') 'c' ...
%             num2str(DAPIChannel,'%01d') '.tif']) ;
        
        progress_3 = [sprintf('\n\t\t') 'Calculating total fluorescence and lengths of cells...'] ;
        fprintf(1,progress_3) ;
        [CellArea,AvgCellFluor,CellLength,backgroundfluor] = FluorMeasure_indivZ(sr,LcFull) ;

        progress_4 = [sprintf('\n\t\t') 'Calculating total DAPI signal of cells...'] ;
        fprintf(1,progress_4) ;
        [~, DAPIpp, ~, ~] = FluorMeasure_indivZ(DAPI_sr,LcFull) ;           
        
        % imaginary spots will be discarded
        
        Num1 = max(LcFull(:)) ;
        for n_cell = 1:1:Num1
            
            cells_in_group = [cells_in_group ; [n_cell n_frame zeros(1,11)]] ;

            % Number of spots recognized in the cell.
            cells_in_group(end,3) = sum(clusterdata_low(:,1)==n_cell) ;
            cells_in_group(end,4) = sum(clusterdata_high(:,1)==n_cell) ;
        
            % Total spots intensity of the cell.
            cells_in_group(end,5) = sum(clusterdata_low(clusterdata_low(:,1)==n_cell,3)) ;
            cells_in_group(end,6) = sum(clusterdata_high(clusterdata_high(:,1)==n_cell,3)) ;
            if isnan(cells_in_group(end,6))
                disp('stop !');
            end
            
            % Total cell fluorescence intensity per pixel.
            cells_in_group(end,7) = AvgCellFluor(n_cell) ;
            % cell's background subtracted fluor
            cells_in_group(end,8) = backgroundfluor(n_cell) ;
            
            % Cell area in pixels.
            cells_in_group(end,9) = CellArea(n_cell) ;
            
            % Cell length in pixels.
            cells_in_group(end,10) = CellLength(n_cell) ;
            
            % Total DAPI intensity per pixel.
            cells_in_group(end,11) = DAPIpp(n_cell) ; 
            
            % Total spot intensity/pix in the cell
            try
                cells_in_group(end,12) = sum(clusterdata_low(clusterdata_low(:,1)==n_cell,3)) ...
                    ./ sum(clusterdata_low(clusterdata_low(:,1)==n_cell,4));
            catch E
                cells_in_group(end,12) = -1;
%                 disp(E);
%                 disp('Setting to 0');
            end
            try
                cells_in_group(end,13) = sum(clusterdata_high(clusterdata_high(:,1)==n_cell,3)) ...
                    ./ sum(clusterdata_high(clusterdata_high(:,1)==n_cell,4));
            catch E
                cells_in_group(end,13) = -1;
%                 disp(E);
%                 disp('Setting to 0');
            end

        end 
    end
    
    fprintf(1,sprintf('\n')) ;
    enlistcell_new{n_sample} = cells_in_group ;
end

save([datafolder '\data\' 'spots_data.mat'],...
    'enlistcell_new','-append');
end

%% 1D histogram of cell int/pix
%{
figure('Units','normalized',...
    'color','w',...
    'OuterPosition',[0.0    0.045    0.3*2    .95],...
    'NumberTitle','off','Name','1D spot intensity') ;
binN = 100 ;
ha = tight_subplot(ceil(length(whichsamples)/2), 2, [.02 .05], [.05 .01], [.05 .01]);

iplot = 1 ;
for isample = whichsamples
    
    if ~isempty(clusterlistlow_new{isample})
        
        axes(ha(iplot)) ;
        set(gca, 'xscale','log','fontsize',6)
        
        box on ; hold all ;
        %%%
        minP = min(enlistcell_new{isample}(:,12)) ;
        maxP = max(enlistcell_new{isample}(:,12)) ;
        nbin = logspace(log10(minP),log10(maxP),binN) ;
        %
        [ya,xa] = hist(enlistcell_new{isample}(:,12),nbin) ;
        %
        plot(xa,ya/sum(ya),'-',...
            'Color','r',...
            'LineWidth',1,...
            'DisplayName',[ip.sample{isample}.name '(low threshold)']) ;
        %%%
        y = enlistcell_new{isample}(:,13);
        minP = min(y(~isnan(y))) ;
        maxP = max(y(~isnan(y))) ;
        if isempty(minP) || isempty(maxP)
            %
        else
        nbin = logspace(log10(minP),log10(maxP),binN) ;
        %
        [ya,xa] = hist(y(~isnan(y)),nbin) ;
        %
        plot(xa,ya/sum(ya),'-',...
            'Color','b',...
            'LineWidth',1,...
            'DisplayName',[ip.sample{isample}.name '(high threshold)']) ;
        %%%
        end
        minP = min(enlistcell_new{isample}(:,7)) ;
        maxP = max(enlistcell_new{isample}(:,7)) ;
        nbin = logspace(log10(minP),log10(maxP),binN) ;
        %
        [ya,xa] = hist(enlistcell_new{isample}(:,7),nbin) ;
        %
        plot(xa,ya/sum(ya),'-',...
            'Color','k',...
            'LineWidth',1,...
            'DisplayName',[ip.sample{isample}.name '(whole Cell)']) ;
        %%%
        
        lh=legend('show');
        set(lh,'location','northeast');
        lh.FontSize = 8;

        xlim([1e2 1e4]) ;        
%         ylabel('Probability') ;
        set(gca,'XTickLabel',[]) ;
        iplot = iplot + 1 ;
        set(gca,'XTickLabelMode','auto') ;
        xlabel('Cell intensity/pixel') ;
    end
end
set(gca,'XTickLabelMode','auto') ;
xlabel('Cell intensity/pixel') ;
%}
figure('Units','normalized',...
    'color','w',...
    'OuterPosition',[0.0    0.045    0.3    .6],...
    'NumberTitle','off','Name','1D spot intensity') ;
ha = tight_subplot(3, ceil(length(whichsamples)/3), [.02 .05], [.08 .01], [.05 .01]);
    
for combo = 1:5
    
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
    
    datafolder = [ip.exp.path '\cluster_quantify\Run-GFP-' num2str(low_threshold) ...
        'low' num2str(high_threshold) 'high-cellnormed-thresh-16-Nov-2021\'];
    
    if exist([datafolder '\data\' 'spots_data.mat'],'file')
        load([datafolder '\data\' 'spots_data.mat'],'clusterlistlow_new');
    end
    if exist([datafolder '\data\' 'spots_data.mat'],'file')
        load([datafolder '\data\' 'spots_data.mat'],'clusterlisthigh_new');
    end
    if exist([datafolder '\data\' 'spots_data.mat'],'file')
        load([datafolder '\data\' 'spots_data.mat'],'enlistcell_new');
    end 
    
%     figure('Units','normalized',...
%         'color','w',...
%         'OuterPosition',[0.0    0.045    0.3    .6],...
%         'NumberTitle','off','Name','1D spot intensity') ;
    binN = 100 ;
%     ha = tight_subplot(3, ceil(length(whichsamples)/3), [.02 .05], [.08 .01], [.05 .01]);

    iplot = 1 ;
    for isample = whichsamples

        if ~isempty(clusterlistlow_new{isample})

            axes(ha(iplot)) ;
            set(gca, 'xscale','log','fontsize',6)

            box on ; hold all ;
            
            if combo == 1
                minP = min(enlistcell_new{isample}(:,7)) ;
                maxP = max(enlistcell_new{isample}(:,7)) ;
                nbin = logspace(log10(minP),log10(maxP),binN) ;
                %
                [ya,xa] = hist(enlistcell_new{isample}(:,7),nbin) ;
                %
                plot(xa,ya/sum(ya),'-',...
                    'Color','k',...
                    'LineWidth',1,...
                    'DisplayName',[short_name{isample} ' (whole Cell)']) ;
            end

            minP = min(enlistcell_new{isample}(:,12)) ;
            maxP = max(enlistcell_new{isample}(:,12)) ;
            nbin = logspace(log10(minP),log10(maxP),binN) ;
            %
            [ya,xa] = hist(enlistcell_new{isample}(:,12),nbin) ;
            %
            plot(xa,ya/sum(ya),'-',...
                'Color',[1 0 0] * (combo + 5 - 2) / (2*5 - 2),...
                'LineWidth',1,...
                'DisplayName',[short_name{isample} '(low threshold:' num2str(low_threshold) ')']) ;
            %%%
            y = enlistcell_new{isample}(:,13);
            minP = min(y(~isnan(y))) ;
            maxP = max(y(~isnan(y))) ;
            if isempty(minP) || isempty(maxP)
                iplot = iplot + 1 ;
                continue;
            end
            nbin = logspace(log10(minP),log10(maxP),binN) ;
            %
            [ya,xa] = hist(y(~isnan(y)),nbin) ;
            %
            plot(xa,ya/sum(ya),'-',...
                'Color',[0 0 1] * (combo + 5 - 2) / (2*5 - 2),...
                'LineWidth',1,...
                'DisplayName',[short_name{isample} '(high threshold:' num2str(high_threshold) ')']) ;

            lh=legend('show');
            set(lh,'location','northwest');
            lh.FontSize = 8;

            xlim([6e2 2e3]) ;        
    %         ylabel('Probability') ;
            set(gca,'XTickLabel',[]) ;
            iplot = iplot + 1 ;
%             set(gca,'XTickLabelMode','auto') ;
%             xlabel('Spot intensity/pixel') ;
        end
    end
    set(gca,'XTickLabelMode','auto') ;
    xlabel('Cell intensity/pixel') ;
    
%     savefig(gcf, [datafolder '\fit\1D_spot_intensity.fig']);
%     saveas(gcf, [datafolder '\fit\1D_spot_intensity.png']);
end
    savefig(gcf, [pwd '\cluster_rec\1D_cell_intensity_per_pix.fig']);
    saveas(gcf, [pwd '\cluster_rec\1D_cell_intensity_per_pix.png']);
    
%% 1D histogram of spot/cell
%{
figure('Units','normalized',...
    'color','w',...
    'OuterPosition',[0.0    0.045    0.3*2    .95],...
    'NumberTitle','off','Name','1D spot intensity') ;
binN = 100 ;
ha = tight_subplot(ceil(length(whichsamples)/2), 2, [.02 .05], [.05 .01], [.05 .01]);

iplot = 1 ;
for isample = whichsamples
    
    if ~isempty(clusterlistlow_new{isample})
        
        axes(ha(iplot)) ;
        set(gca, 'fontsize',6)

%         subplot(length(whichsamples),1,iplot,...
%             'XScale','log','FontSize',8) ;
        
        box on ; hold all ;
        
        minP = min(enlistcell_new{isample}(:,3)) ;
        maxP = max(enlistcell_new{isample}(:,3)) ;
        if isempty(minP) || isempty(maxP)
            %
        else
        nbin = linspace(minP,maxP,round(maxP-minP)+1) ;
        
        [ya,xa] = hist(enlistcell_new{isample}(:,3),nbin) ;
        
        bar(round(xa)+.15,ya/sum(ya),...
            'FaceColor','r',...
            'LineWidth',1,...
            'barwidth', .5,...
            'DisplayName',[ip.sample{isample}.name '(low threshold)']) ;
        %%%
        end
        minP = min(enlistcell_new{isample}(:,4)) ;
        maxP = max(enlistcell_new{isample}(:,4)) ;
        if isempty(minP) || isempty(maxP)
            %
        else
        nbin = linspace(minP,maxP,binN) ;
        
        [ya,xa] = hist(enlistcell_new{isample}(:,4),round(maxP-minP)+1) ;
        
        bar(round(xa)-.15,ya/sum(ya),...
            'FaceColor','b',...
            'LineWidth',1,...
            'barwidth', .5,...
            'DisplayName',[ip.sample{isample}.name '(high threshold)']) ;
        %%%
        end
        lh=legend('show');
        set(lh,'location','northeast');
        lh.FontSize = 8;

        xlim([-1 9]) ;        
%         ylabel('Probability') ;
        set(gca,'XTickLabel',[]) ;
        iplot = iplot + 1 ;
        set(gca,'XTickLabelMode','auto') ;
        xlabel('Clusters recognized per cell') ;
    end
end
set(gca,'XTickLabelMode','auto') ;
xlabel('Clusters recognized per cell') ;
%}
% figure('Units','normalized',...
%     'color','w',...
%     'OuterPosition',[0.0    0.045    0.95    .7],...
%     'NumberTitle','off','Name','1D spot intensity') ;
fig1 = figure('Units','normalized',...
    'color','w',...
    'OuterPosition',[0.0    0.045    0.5    .95/2]); 
% ha = subplot(3,1,1);%ha = tight_subplot(3, 1, [.02 .05], [.08 .01], [.05 .01]);
    
for combo = 3%1:5
    
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
    
    datafolder = [ip.exp.path '\cluster_quantify\Run-GFP-' num2str(low_threshold) ...
        'low' num2str(high_threshold) 'high-cellnormed-thresh-07-Dec-2021\'];
    
    if exist([datafolder '\data\' 'spots_data.mat'],'file')
        load([datafolder '\data\' 'spots_data.mat'],'clusterlistlow_new');
    end
    if exist([datafolder '\data\' 'spots_data.mat'],'file')
        load([datafolder '\data\' 'spots_data.mat'],'clusterlisthigh_new');
    end
    if exist([datafolder '\data\' 'spots_data.mat'],'file')
        load([datafolder '\data\' 'spots_data.mat'],'enlistcell_new');
    end 
    
%     figure('Units','normalized',...
%         'color','w',...
%         'OuterPosition',[0.0    0.045    0.3    .6],...
%         'NumberTitle','off','Name','1D spot intensity') ;
    binN = 100 ;
%     ha = tight_subplot(3, ceil(length(whichsamples)/3), [.02 .05], [.08 .01], [.05 .01]);

    iplot = 1 ;
    for isample = whichsamples

        if ~isempty(clusterlistlow_new{isample})

%             axes(ha((iplot)*1 + 0*combo)) ;
            subplot(1,3,iplot);
            title({['rpoc-GFP Clust/cell. Mean +/- SEM:'],...
                ['low: ' num2str(round(mean(enlistcell_new{isample}(:,3))*100)/100) ...
                ' +/- ' num2str(round(std(enlistcell_new{isample}(:,3))/sqrt(length(enlistcell_new{isample}(:,3)))*100)/100)],...
                ['high: ' num2str(round(mean(enlistcell_new{isample}(:,4))*100)/100) ...
                ' +/- ' num2str(round(std(enlistcell_new{isample}(:,3))/sqrt(length(enlistcell_new{isample}(:,4)))*100)/100)]});
%             set(gca, 'xscale','linear','fontsize',6)

            box on ; hold all ;
            
            minP = min(enlistcell_new{isample}(:,3)) ;
            maxP = max(enlistcell_new{isample}(:,3)) ;
            if isempty(minP) || isempty(maxP)
                %
            else
            nbin = linspace(minP,maxP,round(maxP-minP)+1) ;

            [ya,xa] = hist(enlistcell_new{isample}(:,3),nbin) ;

            bar(round(xa)+.15,ya/sum(ya),...
                'FaceColor','r',...
                'LineWidth',1,...
                'barwidth', .5,...
                'DisplayName',[short_name{isample} '(low thresh:' num2str(low_threshold) ...
                ')']) ;
            %%%
            end
            minP = min(enlistcell_new{isample}(:,4)) ;
            maxP = max(enlistcell_new{isample}(:,4)) ;
            if isempty(minP) || isempty(maxP)
                %
            else
            nbin = linspace(minP,maxP,binN) ;

            [ya,xa] = hist(enlistcell_new{isample}(:,4),round(maxP-minP)+1) ;

            bar(round(xa)-.15,ya/sum(ya),...
                'FaceColor','b',...
                'LineWidth',1,...
                'barwidth', .5,...
                'DisplayName',[short_name{isample} '(high thresh:' num2str(high_threshold) ...
                ')']) ;
            %%%
            end
            lh=legend('show');
            set(lh,'location','northeast');
%             lh.FontSize = 8;

            xlim([-0.5 6.5]) ;
            ylim([0 .7]);
            ylabel('Probability') ;
            set(gca,'XTickLabel',[]) ;
            iplot = iplot + 1 ;
%             set(gca,'XTickLabelMode','auto') ;
%             xlabel('Spot intensity/pixel') ;
            grid on;
                set(gca,'XTickLabelMode','auto') ;
    xlabel('clusters per cell') ;
        end
    end
%     set(gca,'XTickLabelMode','auto') ;
%     xlabel('clusters per cell') ;
    
%     savefig(gcf, [datafolder '\fit\1D_spot_intensity.fig']);
%     saveas(gcf, [datafolder '\fit\1D_spot_intensity.png']);
end
    savefig(gcf, [pwd '\cluster_rec\1D_cluster_per_cell.fig']);
    saveas(gcf, [pwd '\cluster_rec\1D_cluster_per_cell.png']);
    
%% 1D histogram of spot int
%{
figure('Units','normalized',...
    'color','w',...
    'OuterPosition',[0.0    0.045    0.3*2    .95],...
    'NumberTitle','off','Name','1D spot intensity') ;
binN = 100 ;
ha = tight_subplot(ceil(length(whichsamples)/2), 2, [.02 .05], [.05 .01], [.05 .01]);

iplot = 1 ;
for isample = whichsamples
    
    if ~isempty(clusterlistlow_new{isample})
        
        axes(ha(iplot)) ;
        set(gca, 'fontsize',6)

%         subplot(length(whichsamples),1,iplot,...
%             'XScale','log','FontSize',8) ;
        
        box on ; hold all ;
        
        minP = min(enlistcell_new{isample}(:,3)) ;
        maxP = max(enlistcell_new{isample}(:,3)) ;
        if isempty(minP) || isempty(maxP)
            %
        else
        nbin = linspace(minP,maxP,round(maxP-minP)+1) ;
        
        [ya,xa] = hist(enlistcell_new{isample}(:,3),nbin) ;
        
        bar(round(xa)+.15,ya/sum(ya),...
            'FaceColor','r',...
            'LineWidth',1,...
            'barwidth', .5,...
            'DisplayName',[ip.sample{isample}.name '(low threshold)']) ;
        %%%
        end
        minP = min(enlistcell_new{isample}(:,4)) ;
        maxP = max(enlistcell_new{isample}(:,4)) ;
        if isempty(minP) || isempty(maxP)
            %
        else
        nbin = linspace(minP,maxP,binN) ;
        
        [ya,xa] = hist(enlistcell_new{isample}(:,4),round(maxP-minP)+1) ;
        
        bar(round(xa)-.15,ya/sum(ya),...
            'FaceColor','b',...
            'LineWidth',1,...
            'barwidth', .5,...
            'DisplayName',[ip.sample{isample}.name '(high threshold)']) ;
        %%%
        end
        lh=legend('show');
        set(lh,'location','northeast');
        lh.FontSize = 8;

        xlim([-1 9]) ;        
%         ylabel('Probability') ;
        set(gca,'XTickLabel',[]) ;
        iplot = iplot + 1 ;
        set(gca,'XTickLabelMode','auto') ;
        xlabel('Clusters recognized per cell') ;
    end
end
set(gca,'XTickLabelMode','auto') ;
xlabel('Clusters recognized per cell') ;
%}
% figure('Units','normalized',...
%     'color','w',...
%     'OuterPosition',[0.0    0.045    0.95    .7],...
%     'NumberTitle','off','Name','1D spot intensity') ;
fig1 = figure('Units','normalized',...
    'color','w',...
    'OuterPosition',[0.0    0.045    0.5    .95/2]); 
% ha = subplot(3,1,1);%ha = tight_subplot(3, 1, [.02 .05], [.08 .01], [.05 .01]);
    
for combo = 3%1:5
    
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
    
    datafolder = [ip.exp.path '\cluster_quantify\Run-GFP-' num2str(low_threshold) ...
        'low' num2str(high_threshold) 'high-cellnormed-thresh-07-Dec-2021\'];
    
    if exist([datafolder '\data\' 'spots_data.mat'],'file')
        load([datafolder '\data\' 'spots_data.mat'],'clusterlistlow_new');
    end
    if exist([datafolder '\data\' 'spots_data.mat'],'file')
        load([datafolder '\data\' 'spots_data.mat'],'clusterlisthigh_new');
    end
    if exist([datafolder '\data\' 'spots_data.mat'],'file')
        load([datafolder '\data\' 'spots_data.mat'],'enlistcell_new');
    end 
    
%     figure('Units','normalized',...
%         'color','w',...
%         'OuterPosition',[0.0    0.045    0.3    .6],...
%         'NumberTitle','off','Name','1D spot intensity') ;
    binN = 20 ;
%     ha = tight_subplot(3, ceil(length(whichsamples)/3), [.02 .05], [.08 .01], [.05 .01]);

    iplot = 1 ;
    for isample = whichsamples

        if ~isempty(clusterlistlow_new{isample})

%             axes(ha((iplot)*1 + 0*combo)) ;
            subplot(1,3,iplot);
            title({['rpoc-GFP Clust/cell. Mean +/- SEM:'],...
                ['low: ' num2str(round(mean(clusterlistlow_new{isample}(:,3))/1000)*1000) ...
                ' +/- ' num2str(round(std(clusterlistlow_new{isample}(:,3))/sqrt(length(clusterlistlow_new{isample}(:,3)))/1000)*1000)],...
                ['high: ' num2str(round(mean(clusterlisthigh_new{isample}(:,3))/1000)*1000) ...
                ' +/- ' num2str(round(std(clusterlisthigh_new{isample}(:,3))/sqrt(length(clusterlisthigh_new{isample}(:,3)))/1000)*1000)]});
%             set(gca, 'xscale','linear','fontsize',6)

            box on ; hold all ;
            
            minP = min(clusterlistlow_new{isample}(:,3)) ;
            maxP = max(clusterlistlow_new{isample}(:,3)) ;
            if isempty(minP) || isempty(maxP)
                %
            else
%             nbin = linspace(minP,maxP,binN) ;
            nbin = linspace(0,10^6,binN) ;

            [ya,xa] = hist(clusterlistlow_new{isample}(:,3),nbin) ;

            bar(round(xa)+.15*(nbin(2)-nbin(1)),ya/sum(ya),...
                'FaceColor','r',...
                'LineWidth',1,...
                'barwidth', .5,...
                'DisplayName',[short_name{isample} '(low thresh:' num2str(low_threshold) ...
                ')']) ;
            %%%
            end
            minP = min(clusterlisthigh_new{isample}(:,3)) ;
            maxP = max(clusterlisthigh_new{isample}(:,3)) ;
            if isempty(minP) || isempty(maxP)
                %
            else
%             nbin = linspace(minP,maxP,binN) ;
            nbin = linspace(0,10^6,binN) ;

            [ya,xa] = hist(clusterlisthigh_new{isample}(:,3),nbin) ;

            bar(round(xa)-.15*(nbin(2)-nbin(1)),ya/sum(ya),...
                'FaceColor','b',...
                'LineWidth',1,...
                'barwidth', .5,...
                'DisplayName',[short_name{isample} '(high thresh:' num2str(high_threshold) ...
                ')']) ;
            %%%
            end
            lh=legend('show');
            set(lh,'location','northeast');
%             lh.FontSize = 8;

            xlim([0 1]*10^6) ;
            ylim([0 .25]);
            ylabel('Probability') ;
            set(gca,'XTickLabel',[]) ;
            iplot = iplot + 1 ;
%             set(gca,'XTickLabelMode','auto') ;
%             xlabel('Spot intensity/pixel') ;
            grid on;
                set(gca,'XTickLabelMode','auto') ;
    xlabel('clusters per cell') ;
        end
    end
%     set(gca,'XTickLabelMode','auto') ;
%     xlabel('clusters per cell') ;
    
%     savefig(gcf, [datafolder '\fit\1D_spot_intensity.fig']);
%     saveas(gcf, [datafolder '\fit\1D_spot_intensity.png']);
end
savefig(gcf, [pwd '\cluster_rec\1D_cluster_int.fig']);
saveas(gcf, [pwd '\cluster_rec\1D_cluster_int.png']);

%% 1D histogram of cell int
%{
figure('Units','normalized',...
    'color','w',...
    'OuterPosition',[0.0    0.045    0.3*2    .95],...
    'NumberTitle','off','Name','1D spot intensity') ;
binN = 100 ;
ha = tight_subplot(ceil(length(whichsamples)/2), 2, [.02 .05], [.05 .01], [.05 .01]);

iplot = 1 ;
for isample = whichsamples
    
    if ~isempty(clusterlistlow_new{isample})
        
        axes(ha(iplot)) ;
        set(gca, 'xscale','log','fontsize',6)

        
        box on ; hold all ;
        
        %%%
        minP = min(enlistcell_new{isample}(:,5)) ;
        maxP = max(enlistcell_new{isample}(:,5)) ;
        if minP == 0
            minP = 1;
        end
        nbin = logspace(log10(minP),log10(maxP),binN) ;
        
        [ya,xa] = hist(enlistcell_new{isample}(:,5),nbin) ;
        
        plot(xa,ya/sum(ya),'-',...
            'Color','r',...
            'LineWidth',1,...
            'DisplayName',[ip.sample{isample}.name '(low threshold)']) ;
        %%%
        minP = min(enlistcell_new{isample}(:,6)) ;
        if minP == 0
            minP = 1;
        end
        maxP = max(enlistcell_new{isample}(:,6)) ;
        if maxP ~= 0
        nbin = logspace(log10(minP),log10(maxP),binN) ;
        
        [ya,xa] = hist(enlistcell_new{isample}(:,6),nbin) ;
        
        plot(xa,ya/sum(ya),'-',...
            'Color','b',...
            'LineWidth',1,...
            'DisplayName',[ip.sample{isample}.name '(high threshold)']) ;
        end
        %%%
        
%         [ya,xa] = hist(spotlist_new{isample}(~inside_,5),nbin) ;
% %         
%         plot(xa,ya/sum(ya),'-',...
%             'Color',0.7*[1 1 1],...
%             'LineWidth',1,...
%             'DisplayName',[ip.sample{isample}.name '(outside)']) ;
        lh=legend('show');
        set(lh,'location','northeast');
        lh.FontSize = 8;

        xlim([1e4 1e7]) ;        
%         ylabel('Probability') ;
        set(gca,'XTickLabel',[]) ;
        iplot = iplot + 1 ;
        set(gca,'XTickLabelMode','auto') ;
        xlabel('All spot intensity within the cell') ;
    end
end
set(gca,'XTickLabelMode','auto') ;
xlabel('All spot intensity within the cell') ;
%}
figure('Units','normalized',...
    'color','w',...
    'OuterPosition',[0.0    0.045    0.3    .6],...
    'NumberTitle','off','Name','1D spot intensity') ;
ha = tight_subplot(3, ceil(length(whichsamples)/3), [.02 .05], [.08 .01], [.05 .01]);
    
for combo = 1:5
    
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
    
    datafolder = [ip.exp.path '\cluster_quantify\Run-GFP-' num2str(low_threshold) ...
        'low' num2str(high_threshold) 'high-cellnormed-thresh-16-Nov-2021\'];
    
    if exist([datafolder '\data\' 'spots_data.mat'],'file')
        load([datafolder '\data\' 'spots_data.mat'],'clusterlistlow_new');
    end
    if exist([datafolder '\data\' 'spots_data.mat'],'file')
        load([datafolder '\data\' 'spots_data.mat'],'clusterlisthigh_new');
    end
    if exist([datafolder '\data\' 'spots_data.mat'],'file')
        load([datafolder '\data\' 'spots_data.mat'],'enlistcell_new');
    end 
    
%     figure('Units','normalized',...
%         'color','w',...
%         'OuterPosition',[0.0    0.045    0.3    .6],...
%         'NumberTitle','off','Name','1D spot intensity') ;
    binN = 100 ;
%     ha = tight_subplot(3, ceil(length(whichsamples)/3), [.02 .05], [.08 .01], [.05 .01]);

    iplot = 1 ;
    for isample = whichsamples

        if ~isempty(clusterlistlow_new{isample})
            axes(ha(iplot)) ;
            set(gca, 'xscale','log','fontsize',6)

            box on ; hold all ;
            
            if combo == 1
                minP = min(enlistcell_new{isample}(:,7) .* enlistcell_new{isample}(:,9)) ;
                maxP = max(enlistcell_new{isample}(:,7) .* enlistcell_new{isample}(:,9)) ;
                nbin = logspace(log10(minP),log10(maxP),binN) ;
                %
                [ya,xa] = hist(enlistcell_new{isample}(:,7) .* enlistcell_new{isample}(:,9),nbin) ;
                %
                plot(xa,ya/sum(ya),'-',...
                    'Color','k',...
                    'LineWidth',1,...
                    'DisplayName',[short_name{isample} ' (whole Cell)']) ;
            end

            minP = min(max(enlistcell_new{isample}(:,5),1)) ;
            maxP = max(enlistcell_new{isample}(:,5)) ;
            nbin = logspace(log10(minP),log10(maxP),binN) ;
            %
            [ya,xa] = hist(enlistcell_new{isample}(:,5),nbin) ;
            %
            plot(xa,ya/sum(ya),'-',...
                'Color',[1 0 0] * (combo + 5 - 2) / (2*5 - 2),...
                'LineWidth',1,...
                'DisplayName',[short_name{isample} '(low threshold:' num2str(low_threshold) ')']) ;
            %%%
            y = enlistcell_new{isample}(:,6);
            minP = min(max(y(~isnan(y)),1)) ;
            maxP = max(y(~isnan(y))) ;
            if isempty(minP) || isempty(maxP)
                iplot = iplot + 1 ;
                continue;
            end
            nbin = logspace(log10(minP),log10(maxP),binN) ;
            %
            [ya,xa] = hist(y(~isnan(y)),nbin) ;
            %
            plot(xa,ya/sum(ya),'-',...
                'Color',[0 0 1] * (combo + 5 - 2) / (2*5 - 2),...
                'LineWidth',1,...
                'DisplayName',[short_name{isample} '(high threshold:' num2str(high_threshold) ')']) ;

            lh=legend('show');
            set(lh,'location','northwest');
            lh.FontSize = 8;

            xlim([1e3 1e7]) ;        
    %         ylabel('Probability') ;
            set(gca,'XTickLabel',[]) ;
            iplot = iplot + 1 ;
%             set(gca,'XTickLabelMode','auto') ;
%             xlabel('Spot intensity/pixel') ;
        end
    end
    set(gca,'XTickLabelMode','auto') ;
    xlabel('Cell intensity') ;
    
%     savefig(gcf, [datafolder '\fit\1D_spot_intensity.fig']);
%     saveas(gcf, [datafolder '\fit\1D_spot_intensity.png']);
end
    savefig(gcf, [pwd '\cluster_rec\1D_cell_intensity.fig']);
    saveas(gcf, [pwd '\cluster_rec\1D_cell_intensity.png']);
=======
% Initialize variables

clc
clear

ip = exp_20211029_InitializeExp();

spl2cnum = struct() ; 

spl2cnum.GFP = 4*ones(1,max(ip.exp.sampleList)) ;
spl2cnum.DAPI = 5*ones(1,max(ip.exp.sampleList)) ;

whichchanel = 'GFP' ;  

whichspl2cnum = spl2cnum.(whichchanel) ; 

% datafolder = [pwd '\cluster_quantify\Run-TXRED-.4low.7high-cellnormed-thresh-02-Nov-2021'] ;

whichsamples = ip.exp.sampleList ;

plot_col = 5 ;
plot_row = 2 ;

% if exist([datafolder '\data\' 'spots_data.mat'],'file')
%     load([datafolder '\data\' 'spots_data.mat'],'clusterlistlow_new');
% end
% if exist([datafolder '\data\' 'spots_data.mat'],'file')
%     load([datafolder '\data\' 'spots_data.mat'],'clusterlisthigh_new');
% end 
short_name = {'1mM IPTG', '10 \muM IPTG', '0mM IPTG'};

%% group clusterdata_X.mat of individual frames
for combo = 3%1:5
    
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
    
    datafolder = [ip.exp.path '\cluster_quantify\Run-GFP-' num2str(low_threshold) ...
        'low' num2str(high_threshold) 'high-cellnormed-thresh-07-Dec-2021\'];

whichsamples = ip.exp.sampleList ;

plot_col = 5 ;
plot_row = 3 ;

if exist([datafolder '\data\' 'spots_data.mat'],'file')
    load([datafolder '\data\' 'spots_data.mat'],'clusterlistlow_new');
end
if exist([datafolder '\data\' 'spots_data.mat'],'file')
    load([datafolder '\data\' 'spots_data.mat'],'clusterlisthigh_new');
end 



% only run for frames I segmented.
%ip.sample{2}.idx = 1:8 ;

clusterlistlow_new = cell(1,max(ip.exp.sampleList)) ;
clusterlisthigh_new = cell(1,max(ip.exp.sampleList)) ;

for n_sample = whichsamples
    
    fprintf(1,['Group ' num2str(n_sample) ' of ' num2str(max(ip.exp.sampleList)) '.' sprintf('\n')]);
    
    spots_in_group1 = [] ;
    spots_in_group2 = [] ;
    
    r1 = 0 ;
    r2 = 0 ;
    
    progress_2 = [];
    for n_image = ip.sample{n_sample}.idx

        n_frame = ip.exp.splimg2frm(n_sample,n_image);
        
        % get the channel number from the sample numebr
        SpotChannel = whichspl2cnum(n_sample) ; 
        sr = InitializeSpotRecognitionParameters_c234(ip,n_image,SpotChannel,datafolder); %% update the file name here. 
        
        load([sr.output 'clusterdata_low' num2str(n_frame,'%03d') '.mat']);
        
        if ~isempty(clusterdata_low)
            spots_in_group1 = [spots_in_group1 ; ...
                clusterdata_low(:,2) ...            % Z-slice  (1)
                [1:size(clusterdata_low,1)]' ...    % vector with a list of spot numbers for the frame (2)
                clusterdata_low(:,[3 4 5]) ...      % SpotInt area AvgSpotInt(3 4 5)
                clusterdata_low(:,1) ...    % cellnum (6)
                clusterdata_low(:,[6]) ... % frame (7)
                clusterdata_low(:,[9]) ... % low cluster number (8)
                clusterdata_low(:,[7 8])] ; % average x y position (9 10)
        end
        
        load([sr.output 'clusterdata_high' num2str(n_frame,'%03d') '.mat'],'clusterdata_high');
        
        if ~isempty(clusterdata_high)
            spots_in_group2 = [spots_in_group2 ; ...
                clusterdata_high(:,2) ...         % Z-slice  (1)
                [1:size(clusterdata_high,1)]' ...  % vector with a list of spot numbers for the frame (2)
                clusterdata_high(:,[3 4 5]) ...   % SpotInt area AvgSpotInt(3 4 5)
                clusterdata_high(:,1) ...  % cellnum (6)
                clusterdata_high(:, 6) ... % frame (7)
                clusterdata_high(:, 9) ... % low cluster number (8)
                clusterdata_high(:,[7 8])] ; % average x y position (9 10)
        end
    end
    
    if ~isempty(spots_in_group1)
        clusterlistlow_new{n_sample} = spots_in_group1 ;
    else
        clusterlistlow_new{n_sample} = zeros(1,7) ;
    end
    
    if ~isempty(spots_in_group2)
        clusterlisthigh_new{n_sample} = spots_in_group2 ;
    else
        clusterlisthigh_new{n_sample} = zeros(1,8) ;
    end
    
end

if exist([datafolder '\data\' 'spots_data.mat'],'file')
    save([datafolder '\data\' 'spots_data.mat'],'clusterlistlow_new','-append');
else
    save([datafolder '\data\' 'spots_data.mat'],'clusterlistlow_new');
end
if exist([datafolder '\data\' 'spots_data.mat'],'file')
    save([datafolder '\data\' 'spots_data.mat'],'clusterlisthigh_new','-append');
else
    save([datafolder '\data\' 'spots_data.mat'],'clusterlisthigh_new');
end
end

%% 1D histogram of spot int/pix
figure('Units','normalized',...
    'color','w',...
    'OuterPosition',[0.0    0.045    0.3    .6],...
    'NumberTitle','off','Name','1D spot intensity') ;
ha = tight_subplot(3, ceil(length(whichsamples)/3), [.02 .05], [.08 .01], [.05 .01]);
    
for combo = 1:5
    
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
    
    datafolder = [ip.exp.path '\cluster_quantify\Run-GFP-' num2str(low_threshold) ...
        'low' num2str(high_threshold) 'high-cellnormed-thresh-16-Nov-2021\'];
    
    if exist([datafolder '\data\' 'spots_data.mat'],'file')
        load([datafolder '\data\' 'spots_data.mat'],'clusterlistlow_new');
    end
    if exist([datafolder '\data\' 'spots_data.mat'],'file')
        load([datafolder '\data\' 'spots_data.mat'],'clusterlisthigh_new');
    end 
    
%     figure('Units','normalized',...
%         'color','w',...
%         'OuterPosition',[0.0    0.045    0.3    .6],...
%         'NumberTitle','off','Name','1D spot intensity') ;
    binN = 100 ;
%     ha = tight_subplot(3, ceil(length(whichsamples)/3), [.02 .05], [.08 .01], [.05 .01]);

    iplot = 1 ;
    for isample = whichsamples

        if ~isempty(clusterlistlow_new{isample})

            axes(ha(iplot)) ;
            set(gca, 'xscale','log','fontsize',6)

            box on ; hold all ;

            minP = min(clusterlistlow_new{isample}(:,3) ./ clusterlistlow_new{isample}(:,4)) ;
            maxP = max(clusterlistlow_new{isample}(:,3) ./ clusterlistlow_new{isample}(:,4)) ;
            nbin = logspace(log10(minP),log10(maxP),binN) ;
            %
            [ya,xa] = hist(clusterlistlow_new{isample}(:,3) ./ clusterlistlow_new{isample}(:,4),nbin) ;
            %
            plot(xa,ya/sum(ya),'-',...
                'Color',[1 0 0] * (combo + 5 - 2) / (2*5 - 2),...
                'LineWidth',1,...
                'DisplayName',[short_name{isample} '(low threshold:' num2str(low_threshold) ')']) ;
            %%%
            y = clusterlisthigh_new{isample}(:,3) ./ clusterlisthigh_new{isample}(:,4);
            minP = min(y(~isnan(y))) ;
            maxP = max(y(~isnan(y))) ;
            if isempty(minP) || isempty(maxP)
                iplot = iplot + 1 ;
                continue;
            end
            nbin = logspace(log10(minP),log10(maxP),binN) ;
            %
            [ya,xa] = hist(y(~isnan(y)),nbin) ;
            %
            plot(xa,ya/sum(ya),'-',...
                'Color',[0 0 1] * (combo + 5 - 2) / (2*5 - 2),...
                'LineWidth',1,...
                'DisplayName',[short_name{isample} '(high threshold:' num2str(high_threshold) ')']) ;

            lh=legend('show');
            set(lh,'location','northwest');
            lh.FontSize = 8;

            xlim([.7e3 2e3]) ;        
    %         ylabel('Probability') ;
            set(gca,'XTickLabel',[]) ;
            iplot = iplot + 1 ;
%             set(gca,'XTickLabelMode','auto') ;
%             xlabel('Spot intensity/pixel') ;
        end
    end
    set(gca,'XTickLabelMode','auto') ;
    xlabel('Spot intensity/pixel') ;
    
%     savefig(gcf, [datafolder '\fit\1D_spot_intensity.fig']);
%     saveas(gcf, [datafolder '\fit\1D_spot_intensity.png']);
end
    savefig(gcf, [pwd '\cluster_rec\1D_spot_intensity.fig']);
    saveas(gcf, [pwd '\cluster_rec\1D_spot_intensity.png']);

%% 1D histogram of spot area
%{
figure('Units','normalized',...
    'color','w',...
    'OuterPosition',[0.0    0.045    0.3*2    .95],...
    'NumberTitle','off','Name','1D spot intensity') ;
binN = 100 ;
ha = tight_subplot(ceil(length(whichsamples)/2), 2, [.02 .05], [.05 .01], [.05 .01]);

iplot = 1 ;
for isample = whichsamples
    
    if ~isempty(clusterlistlow_new{isample})
        
        axes(ha(iplot)) ;
        set(gca, 'xscale','log','fontsize',6)

%         subplot(length(whichsamples),1,iplot,...
%             'XScale','log','FontSize',8) ;
        
        box on ; hold all ;
        
        minP = min(clusterlistlow_new{isample}(:,4)) ;
        maxP = max(clusterlistlow_new{isample}(:,4)) ;
        nbin = logspace(log10(minP),log10(maxP),binN) ;
        
        [ya,xa] = hist(clusterlistlow_new{isample}(:,4),nbin) ;
        
        plot(xa,ya/sum(ya),'-',...
            'Color','r',...
            'LineWidth',1,...
            'DisplayName',[ip.sample{isample}.name '(low threshold)']) ;
        %%%
        minP = min(clusterlisthigh_new{isample}(:,4)) ;
        maxP = max(clusterlisthigh_new{isample}(:,4)) ;
        nbin = logspace(log10(minP),log10(maxP),binN) ;
        
        [ya,xa] = hist(clusterlisthigh_new{isample}(:,4),nbin) ;
        
        plot(xa,ya/sum(ya),'-',...
            'Color','b',...
            'LineWidth',1,...
            'DisplayName',[ip.sample{isample}.name '(high threshold)']) ;
        %
        lh=legend('show');
        set(lh,'location','northeast');
        lh.FontSize = 8;

        xlim([1e1 1e4]) ;        
%         ylabel('Probability') ;
        set(gca,'XTickLabel',[]) ;
        iplot = iplot + 1 ;
        set(gca,'XTickLabelMode','auto') ;
        xlabel('Spot area') ;
    end
end
set(gca,'XTickLabelMode','auto') ;
xlabel('Spot area') ;
%}
figure('Units','normalized',...
    'color','w',...
    'OuterPosition',[0.0    0.045    0.3    .6],...
    'NumberTitle','off','Name','1D spot intensity') ;
ha = tight_subplot(3, ceil(length(whichsamples)/3), [.02 .05], [.08 .01], [.05 .01]);
    
for combo = 1:5
    
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
    
    datafolder = [ip.exp.path '\cluster_quantify\Run-GFP-' num2str(low_threshold) ...
        'low' num2str(high_threshold) 'high-cellnormed-thresh-16-Nov-2021\'];
    
    if exist([datafolder '\data\' 'spots_data.mat'],'file')
        load([datafolder '\data\' 'spots_data.mat'],'clusterlistlow_new');
    end
    if exist([datafolder '\data\' 'spots_data.mat'],'file')
        load([datafolder '\data\' 'spots_data.mat'],'clusterlisthigh_new');
    end 
    
%     figure('Units','normalized',...
%         'color','w',...
%         'OuterPosition',[0.0    0.045    0.3    .6],...
%         'NumberTitle','off','Name','1D spot intensity') ;
    binN = 100 ;
%     ha = tight_subplot(3, ceil(length(whichsamples)/3), [.02 .05], [.08 .01], [.05 .01]);

    iplot = 1 ;
    for isample = whichsamples

        if ~isempty(clusterlistlow_new{isample})

            axes(ha(iplot)) ;
            set(gca, 'xscale','log','fontsize',6)

            box on ; hold all ;

            minP = min(clusterlistlow_new{isample}(:,4)) ;
            maxP = max(clusterlistlow_new{isample}(:,4)) ;
            nbin = logspace(log10(minP),log10(maxP),binN) ;
            %
            [ya,xa] = hist(clusterlistlow_new{isample}(:,4),nbin) ;
            %
            plot(xa,ya/sum(ya),'-',...
                'Color',[1 0 0] * (combo + 5 - 2) / (2*5 - 2),...
                'LineWidth',1,...
                'DisplayName',[short_name{isample} '(low threshold:' num2str(low_threshold) ')']) ;
            %%%
            y = clusterlisthigh_new{isample}(:,4);
            minP = min(max(y(~isnan(y)),1)) ;
            maxP = max(y(~isnan(y))) ;
            if isempty(minP) || isempty(maxP)
                iplot = iplot + 1 ;
                continue;
            end
            nbin = logspace(log10(minP),log10(maxP),binN) ;
            %
            [ya,xa] = hist(y(~isnan(y)),nbin) ;
            %
            plot(xa,ya/sum(ya),'-',...
                'Color',[0 0 1] * (combo + 5 - 2) / (2*5 - 2),...
                'LineWidth',1,...
                'DisplayName',[short_name{isample} '(high threshold:' num2str(high_threshold) ')']) ;

            lh=legend('show');
            set(lh,'location','northwest');
            lh.FontSize = 8;

            xlim([1e0 5e3]) ; 
            ylim([0 .1])
    %         ylabel('Probability') ;
            set(gca,'XTickLabel',[]) ;
            iplot = iplot + 1 ;
%             set(gca,'XTickLabelMode','auto') ;
%             xlabel('Spot intensity/pixel') ;
        end
    end
    set(gca,'XTickLabelMode','auto') ;
    xlabel('Spot pixel area') ;
    
%     savefig(gcf, [datafolder '\fit\1D_spot_intensity.fig']);
%     saveas(gcf, [datafolder '\fit\1D_spot_intensity.png']);
end
    savefig(gcf, [pwd '\cluster_rec\1D_spot_area.fig']);
    saveas(gcf, [pwd '\cluster_rec\1D_spot_area.png']);
    
%% 1D histogram of spot int
%{
figure('Units','normalized',...
    'color','w',...
    'OuterPosition',[0.0    0.045    0.3*2    .95],...
    'NumberTitle','off','Name','1D spot intensity') ;
binN = 100 ;
ha = tight_subplot(ceil(length(whichsamples)/2), 2, [.02 .05], [.05 .01], [.05 .01]);

iplot = 1 ;
for isample = whichsamples
    
    if ~isempty(clusterlistlow_new{isample})
        
        axes(ha(iplot)) ;
        set(gca, 'xscale','log','fontsize',6)

        
        box on ; hold all ;
        
        minP = min(clusterlistlow_new{isample}(:,3)) ;
        maxP = max(clusterlistlow_new{isample}(:,3)) ;
        nbin = logspace(log10(minP),log10(maxP),binN) ;
        
        [ya,xa] = hist(clusterlistlow_new{isample}(:,3),nbin) ;
        
        plot(xa,ya/sum(ya),'-',...
            'Color','r',...
            'LineWidth',1,...
            'DisplayName',[ip.sample{isample}.name '(low threshold)']) ;
        %%%
        minP = min(clusterlisthigh_new{isample}(:,3)) ;
        maxP = max(clusterlisthigh_new{isample}(:,3)) ;
        nbin = logspace(log10(minP),log10(maxP),binN) ;
        
        [ya,xa] = hist(clusterlisthigh_new{isample}(:,3),nbin) ;
        
        plot(xa,ya/sum(ya),'-',...
            'Color','b',...
            'LineWidth',1,...
            'DisplayName',[ip.sample{isample}.name '(high threshold)']) ;
        
%         [ya,xa] = hist(spotlist_new{isample}(~inside_,5),nbin) ;
% %         
%         plot(xa,ya/sum(ya),'-',...
%             'Color',0.7*[1 1 1],...
%             'LineWidth',1,...
%             'DisplayName',[ip.sample{isample}.name '(outside)']) ;
        lh=legend('show');
        set(lh,'location','northeast');
        lh.FontSize = 8;

        xlim([1e4 1e7]) ;        
%         ylabel('Probability') ;
        set(gca,'XTickLabel',[]) ;
        iplot = iplot + 1 ;
        set(gca,'XTickLabelMode','auto') ;
        xlabel('Spot intensity') ;
    end
end
set(gca,'XTickLabelMode','auto') ;
xlabel('Spot intensity') ;
%}
figure('Units','normalized',...
    'color','w',...
    'OuterPosition',[0.0    0.045    0.3    .6],...
    'NumberTitle','off','Name','1D spot intensity') ;
ha = tight_subplot(3, ceil(length(whichsamples)/3), [.02 .05], [.08 .01], [.05 .01]);
    
for combo = 1:5
    
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
    
    datafolder = [ip.exp.path '\cluster_quantify\Run-GFP-' num2str(low_threshold) ...
        'low' num2str(high_threshold) 'high-cellnormed-thresh-16-Nov-2021\'];
    
    if exist([datafolder '\data\' 'spots_data.mat'],'file')
        load([datafolder '\data\' 'spots_data.mat'],'clusterlistlow_new');
    end
    if exist([datafolder '\data\' 'spots_data.mat'],'file')
        load([datafolder '\data\' 'spots_data.mat'],'clusterlisthigh_new');
    end 
    
%     figure('Units','normalized',...
%         'color','w',...
%         'OuterPosition',[0.0    0.045    0.3    .6],...
%         'NumberTitle','off','Name','1D spot intensity') ;
    binN = 100 ;
%     ha = tight_subplot(3, ceil(length(whichsamples)/3), [.02 .05], [.08 .01], [.05 .01]);

    iplot = 1 ;
    for isample = whichsamples

        if ~isempty(clusterlistlow_new{isample})

            axes(ha(iplot)) ;
            set(gca, 'xscale','log','fontsize',6)

            box on ; hold all ;

            minP = min(clusterlistlow_new{isample}(:,3));
            maxP = max(clusterlistlow_new{isample}(:,3));
            nbin = logspace(log10(minP),log10(maxP),binN) ;
            %
            [ya,xa] = hist(clusterlistlow_new{isample}(:,3),nbin) ;
            %
            plot(xa,ya/sum(ya),'-',...
                'Color',[1 0 0] * (combo + 5 - 2) / (2*5 - 2),...
                'LineWidth',1,...
                'DisplayName',[short_name{isample} '(low threshold:' num2str(low_threshold) ')']) ;
            %%%
            y = clusterlisthigh_new{isample}(:,3);
            minP = min(max(y(~isnan(y)),1)) ;
            maxP = max(y(~isnan(y))) ;
            if isempty(minP) || isempty(maxP)
                iplot = iplot + 1 ;
                continue;
            end
            nbin = logspace(log10(minP),log10(maxP),binN) ;
            %
            [ya,xa] = hist(y(~isnan(y)),nbin) ;
            %
            plot(xa,ya/sum(ya),'-',...
                'Color',[0 0 1] * (combo + 5 - 2) / (2*5 - 2),...
                'LineWidth',1,...
                'DisplayName',[short_name{isample} '(high threshold:' num2str(high_threshold) ')']) ;

            lh=legend('show');
            set(lh,'location','northwest');
            lh.FontSize = 8;

            xlim([1e3 1e7]) ;        
    %         ylabel('Probability') ;
            set(gca,'XTickLabel',[]) ;
            iplot = iplot + 1 ;
%             set(gca,'XTickLabelMode','auto') ;
%             xlabel('Spot intensity/pixel') ;
        end
    end
    set(gca,'XTickLabelMode','auto') ;
    xlabel('Spot intensity') ;
    
%     savefig(gcf, [datafolder '\fit\1D_spot_intensity.fig']);
%     saveas(gcf, [datafolder '\fit\1D_spot_intensity.png']);
end
    savefig(gcf, [pwd '\cluster_rec\1D_spot_intensity.fig']);
    saveas(gcf, [pwd '\cluster_rec\1D_spot_intensity.png']);
    
%% 2D plot only spot area vs peak height (not changed)
%{
markersize = 2;
colors = [[255 0 0 ];[200 200 200];[99 184 255];]/255;

figure('Units','normalized',...
    'color','w',...
    'OuterPosition',[ 0 0 1 1],...
    'NumberTitle','off','Name','spot area VS peakheight scatter') ;
% figure;
iplot = 1 ;
for isample = whichsamples
    
    % spotlist fields
    % 1  Z-slice
    % 2  spot number in the frame
    % 3  BG
    % 4  peakInt
    % 5  SpotInt
    % 6  pi*major*minor (area of spot)
    % 7  X
    % 8  Y
    % 9  cellnum
    % 10 pos_long_axis (normalized)
    % 11 pos_across (normalized)
    % 12 pos_long_axis_pix
    % 13 pos_across_axis_pix
    % 14 frame
    % 15 Major_axis
    % 16 Minor_axis
    % 17 imaginary spot flag
    
    % 3. spot area vs peak height
    subplot(plot_row , plot_col , iplot); hold on;
    
    inside_ = spotlist_new{isample}(:,9) ~= 0 ;
    plot(spotlist_new{isample}(~inside_,6),spotlist_new{isample}(~inside_,4), ...
        '.','color',0.7*[1 1 1],'markersize',markersize);
    hold on ;
    plot(spotlist_new{isample}(inside_,6),spotlist_new{isample}(inside_,4), ...
        '.','color','r','markersize',markersize);
    
%     if ~strcmp(whichchanel,'CFP')
%         plot([1e0 1e4] , [1 1]*NegPeakThreshold{1},'r--')
%     end
    
    xlabel('Spot area');
    ylabel('Peak height');
    ylim([1e1 1e6]);
    xlim([1e0 1e4]);
    set(gca , 'ytick', [1e1 1e2 1e3 1e4 1e5]) ;
    set(gca , 'xtick', [1e0 1e1 1e2 1e3 1e4]) ;
    xl = get(gca, 'xlim') ;
    
    hold off
    set(gca,'XScale','log','YScale','log')
    title(ip.sample{isample}.name) ;
%     legend([num2str(length(spotlist_new{isample}(:,6))) 'spots'],...
%         'location','northeast')
%     legend('outside','inside') ;
    grid on
    grid minor
    
    iplot = iplot + 1 ;
end
%}

%% heatmap to show the histograms (not changed)

% close all

markersize = 2;
colors = [[255 0 0 ];[200 200 200];[99 184 255];]/255;

figure('Units','normalized',...
    'color','w',...
    'OuterPosition',[0 0 1 1],...
    'NumberTitle','off','Name','spot area VS peakheight heatmap') ;

iplot = 1 ;
for isample = whichsamples
    
    % spotlist fields
    % 1  Z-slice
    % 2  spot number in the frame
    % 3  BG
    % 4  peakInt
    % 5  SpotInt
    % 6  pi*major*minor (area of spot)
    % 7  X
    % 8  Y
    % 9  cellnum
    % 10 pos_long_axis (normalized)
    % 11 pos_across (normalized)
    % 12 pos_long_axis_pix
    % 13 pos_across_axis_pix
    % 14 frame
    % 15 Major_axis
    % 16 Minor_axis
    % 17 imaginary spot flag
    
    
    % 3. spot area vs peak height
    subplot(plot_row , plot_col , iplot); hold on;
    xedges = logspace(0,4,100);
    yedges = logspace(1,5,100);
    edges = {xedges,yedges};
    H = hist3(spotlist_new{isample}(:,[6 4]),'Edges',edges);
    imagesc(H'); colormap jet;
    
    xlabel('Spot area');
    ylabel('Peak height');
    set(gca,'XTick',[0 25 50 75 100],'XTickLabel',{'0','1','2','3','4'},...
        'YTick',[0 25 50 75 100],'YTickLabel',{'1','2','3','4','5'});
    xlabel('log_{10}( Spot area )');
    ylabel('log_{10}( Peak height )');
    title(ip.sample{isample}.name) ;
    xlim([0 100]);
    ylim([0 100]);
    iplot = iplot + 1 ;
end

%% % 1. Re-group cluster data into cell data
for combo = 3%1:5
    
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
    
    datafolder = [ip.exp.path '\cluster_quantify\Run-GFP-' num2str(low_threshold) ...
        'low' num2str(high_threshold) 'high-cellnormed-thresh-07-Dec-2021\'];
    
enlistcell_new = cell(1,max(ip.exp.sampleList)) ;

for n_sample = 1:3

    fprintf(1,['Sample ' num2str(n_sample) ' of ' num2str(max(ip.exp.sampleList)) '.' sprintf('\n')]);
    
    cells_in_group = [] ;
    
    progress_2 = [];
    progress_3 = [];
    progress_4 = [];
    
    for n_image = ip.sample{n_sample}.idx
        
        n_frame = ip.exp.splimg2frm(n_sample,n_image);
        
        for d_=1:1:size(progress_2,2)+size(progress_3,2)+size(progress_4,2) ; fprintf(1,'\b') ; end
        progress_2 = [sprintf('\t') 'Image ' num2str(n_image) ' in ' ...
            num2str(min(ip.sample{n_sample}.idx)) '-' num2str(max(ip.sample{n_sample}.idx)) '.'] ;
        fprintf(1,progress_2) ;
        
        % determine the channel number
        DAPIChannel = spl2cnum.DAPI(n_sample);
        SpotChannel = spl2cnum.(whichchanel)(n_sample) ;
        sr = InitializeSpotRecognitionParameters_c234(ip,n_frame,SpotChannel,datafolder);
        DAPI_sr = InitializeSpotRecognitionParameters_c234(ip,n_frame,DAPIChannel,datafolder);
        % already in the gated data folder (no need for thresh or gate)
        load([datafolder '\data\' ...
            'clusterdata_low' num2str(n_frame,'%03d') '.mat'],'clusterdata_low');
        load([datafolder '\data\' ...
            'clusterdata_high' num2str(n_frame,'%03d') '.mat'],'clusterdata_high');
        
        load([sr.seg.dir sr.seg.name],'LcFull') ;

%         cell_DAPI = imread([sr.image.fullname 'z' ...
%             num2str(median(sr.image.zrange),'%01d') 'c' ...
%             num2str(DAPIChannel,'%01d') '.tif']) ;
        
        progress_3 = [sprintf('\n\t\t') 'Calculating total fluorescence and lengths of cells...'] ;
        fprintf(1,progress_3) ;
        [CellArea,AvgCellFluor,CellLength,backgroundfluor] = FluorMeasure_indivZ(sr,LcFull) ;

        progress_4 = [sprintf('\n\t\t') 'Calculating total DAPI signal of cells...'] ;
        fprintf(1,progress_4) ;
        [~, DAPIpp, ~, ~] = FluorMeasure_indivZ(DAPI_sr,LcFull) ;           
        
        % imaginary spots will be discarded
        
        Num1 = max(LcFull(:)) ;
        for n_cell = 1:1:Num1
            
            cells_in_group = [cells_in_group ; [n_cell n_frame zeros(1,11)]] ;

            % Number of spots recognized in the cell.
            cells_in_group(end,3) = sum(clusterdata_low(:,1)==n_cell) ;
            cells_in_group(end,4) = sum(clusterdata_high(:,1)==n_cell) ;
        
            % Total spots intensity of the cell.
            cells_in_group(end,5) = sum(clusterdata_low(clusterdata_low(:,1)==n_cell,3)) ;
            cells_in_group(end,6) = sum(clusterdata_high(clusterdata_high(:,1)==n_cell,3)) ;
            if isnan(cells_in_group(end,6))
                disp('stop !');
            end
            
            % Total cell fluorescence intensity per pixel.
            cells_in_group(end,7) = AvgCellFluor(n_cell) ;
            % cell's background subtracted fluor
            cells_in_group(end,8) = backgroundfluor(n_cell) ;
            
            % Cell area in pixels.
            cells_in_group(end,9) = CellArea(n_cell) ;
            
            % Cell length in pixels.
            cells_in_group(end,10) = CellLength(n_cell) ;
            
            % Total DAPI intensity per pixel.
            cells_in_group(end,11) = DAPIpp(n_cell) ; 
            
            % Total spot intensity/pix in the cell
            try
                cells_in_group(end,12) = sum(clusterdata_low(clusterdata_low(:,1)==n_cell,3)) ...
                    ./ sum(clusterdata_low(clusterdata_low(:,1)==n_cell,4));
            catch E
                cells_in_group(end,12) = -1;
%                 disp(E);
%                 disp('Setting to 0');
            end
            try
                cells_in_group(end,13) = sum(clusterdata_high(clusterdata_high(:,1)==n_cell,3)) ...
                    ./ sum(clusterdata_high(clusterdata_high(:,1)==n_cell,4));
            catch E
                cells_in_group(end,13) = -1;
%                 disp(E);
%                 disp('Setting to 0');
            end

        end 
    end
    
    fprintf(1,sprintf('\n')) ;
    enlistcell_new{n_sample} = cells_in_group ;
end

save([datafolder '\data\' 'spots_data.mat'],...
    'enlistcell_new','-append');
end

%% 1D histogram of cell int/pix
%{
figure('Units','normalized',...
    'color','w',...
    'OuterPosition',[0.0    0.045    0.3*2    .95],...
    'NumberTitle','off','Name','1D spot intensity') ;
binN = 100 ;
ha = tight_subplot(ceil(length(whichsamples)/2), 2, [.02 .05], [.05 .01], [.05 .01]);

iplot = 1 ;
for isample = whichsamples
    
    if ~isempty(clusterlistlow_new{isample})
        
        axes(ha(iplot)) ;
        set(gca, 'xscale','log','fontsize',6)
        
        box on ; hold all ;
        %%%
        minP = min(enlistcell_new{isample}(:,12)) ;
        maxP = max(enlistcell_new{isample}(:,12)) ;
        nbin = logspace(log10(minP),log10(maxP),binN) ;
        %
        [ya,xa] = hist(enlistcell_new{isample}(:,12),nbin) ;
        %
        plot(xa,ya/sum(ya),'-',...
            'Color','r',...
            'LineWidth',1,...
            'DisplayName',[ip.sample{isample}.name '(low threshold)']) ;
        %%%
        y = enlistcell_new{isample}(:,13);
        minP = min(y(~isnan(y))) ;
        maxP = max(y(~isnan(y))) ;
        if isempty(minP) || isempty(maxP)
            %
        else
        nbin = logspace(log10(minP),log10(maxP),binN) ;
        %
        [ya,xa] = hist(y(~isnan(y)),nbin) ;
        %
        plot(xa,ya/sum(ya),'-',...
            'Color','b',...
            'LineWidth',1,...
            'DisplayName',[ip.sample{isample}.name '(high threshold)']) ;
        %%%
        end
        minP = min(enlistcell_new{isample}(:,7)) ;
        maxP = max(enlistcell_new{isample}(:,7)) ;
        nbin = logspace(log10(minP),log10(maxP),binN) ;
        %
        [ya,xa] = hist(enlistcell_new{isample}(:,7),nbin) ;
        %
        plot(xa,ya/sum(ya),'-',...
            'Color','k',...
            'LineWidth',1,...
            'DisplayName',[ip.sample{isample}.name '(whole Cell)']) ;
        %%%
        
        lh=legend('show');
        set(lh,'location','northeast');
        lh.FontSize = 8;

        xlim([1e2 1e4]) ;        
%         ylabel('Probability') ;
        set(gca,'XTickLabel',[]) ;
        iplot = iplot + 1 ;
        set(gca,'XTickLabelMode','auto') ;
        xlabel('Cell intensity/pixel') ;
    end
end
set(gca,'XTickLabelMode','auto') ;
xlabel('Cell intensity/pixel') ;
%}
figure('Units','normalized',...
    'color','w',...
    'OuterPosition',[0.0    0.045    0.3    .6],...
    'NumberTitle','off','Name','1D spot intensity') ;
ha = tight_subplot(3, ceil(length(whichsamples)/3), [.02 .05], [.08 .01], [.05 .01]);
    
for combo = 1:5
    
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
    
    datafolder = [ip.exp.path '\cluster_quantify\Run-GFP-' num2str(low_threshold) ...
        'low' num2str(high_threshold) 'high-cellnormed-thresh-16-Nov-2021\'];
    
    if exist([datafolder '\data\' 'spots_data.mat'],'file')
        load([datafolder '\data\' 'spots_data.mat'],'clusterlistlow_new');
    end
    if exist([datafolder '\data\' 'spots_data.mat'],'file')
        load([datafolder '\data\' 'spots_data.mat'],'clusterlisthigh_new');
    end
    if exist([datafolder '\data\' 'spots_data.mat'],'file')
        load([datafolder '\data\' 'spots_data.mat'],'enlistcell_new');
    end 
    
%     figure('Units','normalized',...
%         'color','w',...
%         'OuterPosition',[0.0    0.045    0.3    .6],...
%         'NumberTitle','off','Name','1D spot intensity') ;
    binN = 100 ;
%     ha = tight_subplot(3, ceil(length(whichsamples)/3), [.02 .05], [.08 .01], [.05 .01]);

    iplot = 1 ;
    for isample = whichsamples

        if ~isempty(clusterlistlow_new{isample})

            axes(ha(iplot)) ;
            set(gca, 'xscale','log','fontsize',6)

            box on ; hold all ;
            
            if combo == 1
                minP = min(enlistcell_new{isample}(:,7)) ;
                maxP = max(enlistcell_new{isample}(:,7)) ;
                nbin = logspace(log10(minP),log10(maxP),binN) ;
                %
                [ya,xa] = hist(enlistcell_new{isample}(:,7),nbin) ;
                %
                plot(xa,ya/sum(ya),'-',...
                    'Color','k',...
                    'LineWidth',1,...
                    'DisplayName',[short_name{isample} ' (whole Cell)']) ;
            end

            minP = min(enlistcell_new{isample}(:,12)) ;
            maxP = max(enlistcell_new{isample}(:,12)) ;
            nbin = logspace(log10(minP),log10(maxP),binN) ;
            %
            [ya,xa] = hist(enlistcell_new{isample}(:,12),nbin) ;
            %
            plot(xa,ya/sum(ya),'-',...
                'Color',[1 0 0] * (combo + 5 - 2) / (2*5 - 2),...
                'LineWidth',1,...
                'DisplayName',[short_name{isample} '(low threshold:' num2str(low_threshold) ')']) ;
            %%%
            y = enlistcell_new{isample}(:,13);
            minP = min(y(~isnan(y))) ;
            maxP = max(y(~isnan(y))) ;
            if isempty(minP) || isempty(maxP)
                iplot = iplot + 1 ;
                continue;
            end
            nbin = logspace(log10(minP),log10(maxP),binN) ;
            %
            [ya,xa] = hist(y(~isnan(y)),nbin) ;
            %
            plot(xa,ya/sum(ya),'-',...
                'Color',[0 0 1] * (combo + 5 - 2) / (2*5 - 2),...
                'LineWidth',1,...
                'DisplayName',[short_name{isample} '(high threshold:' num2str(high_threshold) ')']) ;

            lh=legend('show');
            set(lh,'location','northwest');
            lh.FontSize = 8;

            xlim([6e2 2e3]) ;        
    %         ylabel('Probability') ;
            set(gca,'XTickLabel',[]) ;
            iplot = iplot + 1 ;
%             set(gca,'XTickLabelMode','auto') ;
%             xlabel('Spot intensity/pixel') ;
        end
    end
    set(gca,'XTickLabelMode','auto') ;
    xlabel('Cell intensity/pixel') ;
    
%     savefig(gcf, [datafolder '\fit\1D_spot_intensity.fig']);
%     saveas(gcf, [datafolder '\fit\1D_spot_intensity.png']);
end
    savefig(gcf, [pwd '\cluster_rec\1D_cell_intensity_per_pix.fig']);
    saveas(gcf, [pwd '\cluster_rec\1D_cell_intensity_per_pix.png']);
    
%% 1D histogram of spot/cell
%{
figure('Units','normalized',...
    'color','w',...
    'OuterPosition',[0.0    0.045    0.3*2    .95],...
    'NumberTitle','off','Name','1D spot intensity') ;
binN = 100 ;
ha = tight_subplot(ceil(length(whichsamples)/2), 2, [.02 .05], [.05 .01], [.05 .01]);

iplot = 1 ;
for isample = whichsamples
    
    if ~isempty(clusterlistlow_new{isample})
        
        axes(ha(iplot)) ;
        set(gca, 'fontsize',6)

%         subplot(length(whichsamples),1,iplot,...
%             'XScale','log','FontSize',8) ;
        
        box on ; hold all ;
        
        minP = min(enlistcell_new{isample}(:,3)) ;
        maxP = max(enlistcell_new{isample}(:,3)) ;
        if isempty(minP) || isempty(maxP)
            %
        else
        nbin = linspace(minP,maxP,round(maxP-minP)+1) ;
        
        [ya,xa] = hist(enlistcell_new{isample}(:,3),nbin) ;
        
        bar(round(xa)+.15,ya/sum(ya),...
            'FaceColor','r',...
            'LineWidth',1,...
            'barwidth', .5,...
            'DisplayName',[ip.sample{isample}.name '(low threshold)']) ;
        %%%
        end
        minP = min(enlistcell_new{isample}(:,4)) ;
        maxP = max(enlistcell_new{isample}(:,4)) ;
        if isempty(minP) || isempty(maxP)
            %
        else
        nbin = linspace(minP,maxP,binN) ;
        
        [ya,xa] = hist(enlistcell_new{isample}(:,4),round(maxP-minP)+1) ;
        
        bar(round(xa)-.15,ya/sum(ya),...
            'FaceColor','b',...
            'LineWidth',1,...
            'barwidth', .5,...
            'DisplayName',[ip.sample{isample}.name '(high threshold)']) ;
        %%%
        end
        lh=legend('show');
        set(lh,'location','northeast');
        lh.FontSize = 8;

        xlim([-1 9]) ;        
%         ylabel('Probability') ;
        set(gca,'XTickLabel',[]) ;
        iplot = iplot + 1 ;
        set(gca,'XTickLabelMode','auto') ;
        xlabel('Clusters recognized per cell') ;
    end
end
set(gca,'XTickLabelMode','auto') ;
xlabel('Clusters recognized per cell') ;
%}
% figure('Units','normalized',...
%     'color','w',...
%     'OuterPosition',[0.0    0.045    0.95    .7],...
%     'NumberTitle','off','Name','1D spot intensity') ;
fig1 = figure('Units','normalized',...
    'color','w',...
    'OuterPosition',[0.0    0.045    0.5    .95/2]); 
% ha = subplot(3,1,1);%ha = tight_subplot(3, 1, [.02 .05], [.08 .01], [.05 .01]);
    
for combo = 3%1:5
    
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
    
    datafolder = [ip.exp.path '\cluster_quantify\Run-GFP-' num2str(low_threshold) ...
        'low' num2str(high_threshold) 'high-cellnormed-thresh-07-Dec-2021\'];
    
    if exist([datafolder '\data\' 'spots_data.mat'],'file')
        load([datafolder '\data\' 'spots_data.mat'],'clusterlistlow_new');
    end
    if exist([datafolder '\data\' 'spots_data.mat'],'file')
        load([datafolder '\data\' 'spots_data.mat'],'clusterlisthigh_new');
    end
    if exist([datafolder '\data\' 'spots_data.mat'],'file')
        load([datafolder '\data\' 'spots_data.mat'],'enlistcell_new');
    end 
    
%     figure('Units','normalized',...
%         'color','w',...
%         'OuterPosition',[0.0    0.045    0.3    .6],...
%         'NumberTitle','off','Name','1D spot intensity') ;
    binN = 100 ;
%     ha = tight_subplot(3, ceil(length(whichsamples)/3), [.02 .05], [.08 .01], [.05 .01]);

    iplot = 1 ;
    for isample = whichsamples

        if ~isempty(clusterlistlow_new{isample})

%             axes(ha((iplot)*1 + 0*combo)) ;
            subplot(1,3,iplot);
            title({['rpoc-GFP Clust/cell. Mean +/- SEM:'],...
                ['low: ' num2str(round(mean(enlistcell_new{isample}(:,3))*100)/100) ...
                ' +/- ' num2str(round(std(enlistcell_new{isample}(:,3))/sqrt(length(enlistcell_new{isample}(:,3)))*100)/100)],...
                ['high: ' num2str(round(mean(enlistcell_new{isample}(:,4))*100)/100) ...
                ' +/- ' num2str(round(std(enlistcell_new{isample}(:,3))/sqrt(length(enlistcell_new{isample}(:,4)))*100)/100)]});
%             set(gca, 'xscale','linear','fontsize',6)

            box on ; hold all ;
            
            minP = min(enlistcell_new{isample}(:,3)) ;
            maxP = max(enlistcell_new{isample}(:,3)) ;
            if isempty(minP) || isempty(maxP)
                %
            else
            nbin = linspace(minP,maxP,round(maxP-minP)+1) ;

            [ya,xa] = hist(enlistcell_new{isample}(:,3),nbin) ;

            bar(round(xa)+.15,ya/sum(ya),...
                'FaceColor','r',...
                'LineWidth',1,...
                'barwidth', .5,...
                'DisplayName',[short_name{isample} '(low thresh:' num2str(low_threshold) ...
                ')']) ;
            %%%
            end
            minP = min(enlistcell_new{isample}(:,4)) ;
            maxP = max(enlistcell_new{isample}(:,4)) ;
            if isempty(minP) || isempty(maxP)
                %
            else
            nbin = linspace(minP,maxP,binN) ;

            [ya,xa] = hist(enlistcell_new{isample}(:,4),round(maxP-minP)+1) ;

            bar(round(xa)-.15,ya/sum(ya),...
                'FaceColor','b',...
                'LineWidth',1,...
                'barwidth', .5,...
                'DisplayName',[short_name{isample} '(high thresh:' num2str(high_threshold) ...
                ')']) ;
            %%%
            end
            lh=legend('show');
            set(lh,'location','northeast');
%             lh.FontSize = 8;

            xlim([-0.5 6.5]) ;
            ylim([0 .7]);
            ylabel('Probability') ;
            set(gca,'XTickLabel',[]) ;
            iplot = iplot + 1 ;
%             set(gca,'XTickLabelMode','auto') ;
%             xlabel('Spot intensity/pixel') ;
            grid on;
                set(gca,'XTickLabelMode','auto') ;
    xlabel('clusters per cell') ;
        end
    end
%     set(gca,'XTickLabelMode','auto') ;
%     xlabel('clusters per cell') ;
    
%     savefig(gcf, [datafolder '\fit\1D_spot_intensity.fig']);
%     saveas(gcf, [datafolder '\fit\1D_spot_intensity.png']);
end
    savefig(gcf, [pwd '\cluster_rec\1D_cluster_per_cell.fig']);
    saveas(gcf, [pwd '\cluster_rec\1D_cluster_per_cell.png']);
    
%% 1D histogram of spot int
%{
figure('Units','normalized',...
    'color','w',...
    'OuterPosition',[0.0    0.045    0.3*2    .95],...
    'NumberTitle','off','Name','1D spot intensity') ;
binN = 100 ;
ha = tight_subplot(ceil(length(whichsamples)/2), 2, [.02 .05], [.05 .01], [.05 .01]);

iplot = 1 ;
for isample = whichsamples
    
    if ~isempty(clusterlistlow_new{isample})
        
        axes(ha(iplot)) ;
        set(gca, 'fontsize',6)

%         subplot(length(whichsamples),1,iplot,...
%             'XScale','log','FontSize',8) ;
        
        box on ; hold all ;
        
        minP = min(enlistcell_new{isample}(:,3)) ;
        maxP = max(enlistcell_new{isample}(:,3)) ;
        if isempty(minP) || isempty(maxP)
            %
        else
        nbin = linspace(minP,maxP,round(maxP-minP)+1) ;
        
        [ya,xa] = hist(enlistcell_new{isample}(:,3),nbin) ;
        
        bar(round(xa)+.15,ya/sum(ya),...
            'FaceColor','r',...
            'LineWidth',1,...
            'barwidth', .5,...
            'DisplayName',[ip.sample{isample}.name '(low threshold)']) ;
        %%%
        end
        minP = min(enlistcell_new{isample}(:,4)) ;
        maxP = max(enlistcell_new{isample}(:,4)) ;
        if isempty(minP) || isempty(maxP)
            %
        else
        nbin = linspace(minP,maxP,binN) ;
        
        [ya,xa] = hist(enlistcell_new{isample}(:,4),round(maxP-minP)+1) ;
        
        bar(round(xa)-.15,ya/sum(ya),...
            'FaceColor','b',...
            'LineWidth',1,...
            'barwidth', .5,...
            'DisplayName',[ip.sample{isample}.name '(high threshold)']) ;
        %%%
        end
        lh=legend('show');
        set(lh,'location','northeast');
        lh.FontSize = 8;

        xlim([-1 9]) ;        
%         ylabel('Probability') ;
        set(gca,'XTickLabel',[]) ;
        iplot = iplot + 1 ;
        set(gca,'XTickLabelMode','auto') ;
        xlabel('Clusters recognized per cell') ;
    end
end
set(gca,'XTickLabelMode','auto') ;
xlabel('Clusters recognized per cell') ;
%}
% figure('Units','normalized',...
%     'color','w',...
%     'OuterPosition',[0.0    0.045    0.95    .7],...
%     'NumberTitle','off','Name','1D spot intensity') ;
fig1 = figure('Units','normalized',...
    'color','w',...
    'OuterPosition',[0.0    0.045    0.5    .95/2]); 
% ha = subplot(3,1,1);%ha = tight_subplot(3, 1, [.02 .05], [.08 .01], [.05 .01]);
    
for combo = 3%1:5
    
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
    
    datafolder = [ip.exp.path '\cluster_quantify\Run-GFP-' num2str(low_threshold) ...
        'low' num2str(high_threshold) 'high-cellnormed-thresh-07-Dec-2021\'];
    
    if exist([datafolder '\data\' 'spots_data.mat'],'file')
        load([datafolder '\data\' 'spots_data.mat'],'clusterlistlow_new');
    end
    if exist([datafolder '\data\' 'spots_data.mat'],'file')
        load([datafolder '\data\' 'spots_data.mat'],'clusterlisthigh_new');
    end
    if exist([datafolder '\data\' 'spots_data.mat'],'file')
        load([datafolder '\data\' 'spots_data.mat'],'enlistcell_new');
    end 
    
%     figure('Units','normalized',...
%         'color','w',...
%         'OuterPosition',[0.0    0.045    0.3    .6],...
%         'NumberTitle','off','Name','1D spot intensity') ;
    binN = 20 ;
%     ha = tight_subplot(3, ceil(length(whichsamples)/3), [.02 .05], [.08 .01], [.05 .01]);

    iplot = 1 ;
    for isample = whichsamples

        if ~isempty(clusterlistlow_new{isample})

%             axes(ha((iplot)*1 + 0*combo)) ;
            subplot(1,3,iplot);
            title({['rpoc-GFP Clust/cell. Mean +/- SEM:'],...
                ['low: ' num2str(round(mean(clusterlistlow_new{isample}(:,3))/1000)*1000) ...
                ' +/- ' num2str(round(std(clusterlistlow_new{isample}(:,3))/sqrt(length(clusterlistlow_new{isample}(:,3)))/1000)*1000)],...
                ['high: ' num2str(round(mean(clusterlisthigh_new{isample}(:,3))/1000)*1000) ...
                ' +/- ' num2str(round(std(clusterlisthigh_new{isample}(:,3))/sqrt(length(clusterlisthigh_new{isample}(:,3)))/1000)*1000)]});
%             set(gca, 'xscale','linear','fontsize',6)

            box on ; hold all ;
            
            minP = min(clusterlistlow_new{isample}(:,3)) ;
            maxP = max(clusterlistlow_new{isample}(:,3)) ;
            if isempty(minP) || isempty(maxP)
                %
            else
%             nbin = linspace(minP,maxP,binN) ;
            nbin = linspace(0,10^6,binN) ;

            [ya,xa] = hist(clusterlistlow_new{isample}(:,3),nbin) ;

            bar(round(xa)+.15*(nbin(2)-nbin(1)),ya/sum(ya),...
                'FaceColor','r',...
                'LineWidth',1,...
                'barwidth', .5,...
                'DisplayName',[short_name{isample} '(low thresh:' num2str(low_threshold) ...
                ')']) ;
            %%%
            end
            minP = min(clusterlisthigh_new{isample}(:,3)) ;
            maxP = max(clusterlisthigh_new{isample}(:,3)) ;
            if isempty(minP) || isempty(maxP)
                %
            else
%             nbin = linspace(minP,maxP,binN) ;
            nbin = linspace(0,10^6,binN) ;

            [ya,xa] = hist(clusterlisthigh_new{isample}(:,3),nbin) ;

            bar(round(xa)-.15*(nbin(2)-nbin(1)),ya/sum(ya),...
                'FaceColor','b',...
                'LineWidth',1,...
                'barwidth', .5,...
                'DisplayName',[short_name{isample} '(high thresh:' num2str(high_threshold) ...
                ')']) ;
            %%%
            end
            lh=legend('show');
            set(lh,'location','northeast');
%             lh.FontSize = 8;

            xlim([0 1]*10^6) ;
            ylim([0 .25]);
            ylabel('Probability') ;
            set(gca,'XTickLabel',[]) ;
            iplot = iplot + 1 ;
%             set(gca,'XTickLabelMode','auto') ;
%             xlabel('Spot intensity/pixel') ;
            grid on;
                set(gca,'XTickLabelMode','auto') ;
    xlabel('clusters per cell') ;
        end
    end
%     set(gca,'XTickLabelMode','auto') ;
%     xlabel('clusters per cell') ;
    
%     savefig(gcf, [datafolder '\fit\1D_spot_intensity.fig']);
%     saveas(gcf, [datafolder '\fit\1D_spot_intensity.png']);
end
savefig(gcf, [pwd '\cluster_rec\1D_cluster_int.fig']);
saveas(gcf, [pwd '\cluster_rec\1D_cluster_int.png']);

%% 1D histogram of cell int
%{
figure('Units','normalized',...
    'color','w',...
    'OuterPosition',[0.0    0.045    0.3*2    .95],...
    'NumberTitle','off','Name','1D spot intensity') ;
binN = 100 ;
ha = tight_subplot(ceil(length(whichsamples)/2), 2, [.02 .05], [.05 .01], [.05 .01]);

iplot = 1 ;
for isample = whichsamples
    
    if ~isempty(clusterlistlow_new{isample})
        
        axes(ha(iplot)) ;
        set(gca, 'xscale','log','fontsize',6)

        
        box on ; hold all ;
        
        %%%
        minP = min(enlistcell_new{isample}(:,5)) ;
        maxP = max(enlistcell_new{isample}(:,5)) ;
        if minP == 0
            minP = 1;
        end
        nbin = logspace(log10(minP),log10(maxP),binN) ;
        
        [ya,xa] = hist(enlistcell_new{isample}(:,5),nbin) ;
        
        plot(xa,ya/sum(ya),'-',...
            'Color','r',...
            'LineWidth',1,...
            'DisplayName',[ip.sample{isample}.name '(low threshold)']) ;
        %%%
        minP = min(enlistcell_new{isample}(:,6)) ;
        if minP == 0
            minP = 1;
        end
        maxP = max(enlistcell_new{isample}(:,6)) ;
        if maxP ~= 0
        nbin = logspace(log10(minP),log10(maxP),binN) ;
        
        [ya,xa] = hist(enlistcell_new{isample}(:,6),nbin) ;
        
        plot(xa,ya/sum(ya),'-',...
            'Color','b',...
            'LineWidth',1,...
            'DisplayName',[ip.sample{isample}.name '(high threshold)']) ;
        end
        %%%
        
%         [ya,xa] = hist(spotlist_new{isample}(~inside_,5),nbin) ;
% %         
%         plot(xa,ya/sum(ya),'-',...
%             'Color',0.7*[1 1 1],...
%             'LineWidth',1,...
%             'DisplayName',[ip.sample{isample}.name '(outside)']) ;
        lh=legend('show');
        set(lh,'location','northeast');
        lh.FontSize = 8;

        xlim([1e4 1e7]) ;        
%         ylabel('Probability') ;
        set(gca,'XTickLabel',[]) ;
        iplot = iplot + 1 ;
        set(gca,'XTickLabelMode','auto') ;
        xlabel('All spot intensity within the cell') ;
    end
end
set(gca,'XTickLabelMode','auto') ;
xlabel('All spot intensity within the cell') ;
%}
figure('Units','normalized',...
    'color','w',...
    'OuterPosition',[0.0    0.045    0.3    .6],...
    'NumberTitle','off','Name','1D spot intensity') ;
ha = tight_subplot(3, ceil(length(whichsamples)/3), [.02 .05], [.08 .01], [.05 .01]);
    
for combo = 1:5
    
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
    
    datafolder = [ip.exp.path '\cluster_quantify\Run-GFP-' num2str(low_threshold) ...
        'low' num2str(high_threshold) 'high-cellnormed-thresh-16-Nov-2021\'];
    
    if exist([datafolder '\data\' 'spots_data.mat'],'file')
        load([datafolder '\data\' 'spots_data.mat'],'clusterlistlow_new');
    end
    if exist([datafolder '\data\' 'spots_data.mat'],'file')
        load([datafolder '\data\' 'spots_data.mat'],'clusterlisthigh_new');
    end
    if exist([datafolder '\data\' 'spots_data.mat'],'file')
        load([datafolder '\data\' 'spots_data.mat'],'enlistcell_new');
    end 
    
%     figure('Units','normalized',...
%         'color','w',...
%         'OuterPosition',[0.0    0.045    0.3    .6],...
%         'NumberTitle','off','Name','1D spot intensity') ;
    binN = 100 ;
%     ha = tight_subplot(3, ceil(length(whichsamples)/3), [.02 .05], [.08 .01], [.05 .01]);

    iplot = 1 ;
    for isample = whichsamples

        if ~isempty(clusterlistlow_new{isample})
            axes(ha(iplot)) ;
            set(gca, 'xscale','log','fontsize',6)

            box on ; hold all ;
            
            if combo == 1
                minP = min(enlistcell_new{isample}(:,7) .* enlistcell_new{isample}(:,9)) ;
                maxP = max(enlistcell_new{isample}(:,7) .* enlistcell_new{isample}(:,9)) ;
                nbin = logspace(log10(minP),log10(maxP),binN) ;
                %
                [ya,xa] = hist(enlistcell_new{isample}(:,7) .* enlistcell_new{isample}(:,9),nbin) ;
                %
                plot(xa,ya/sum(ya),'-',...
                    'Color','k',...
                    'LineWidth',1,...
                    'DisplayName',[short_name{isample} ' (whole Cell)']) ;
            end

            minP = min(max(enlistcell_new{isample}(:,5),1)) ;
            maxP = max(enlistcell_new{isample}(:,5)) ;
            nbin = logspace(log10(minP),log10(maxP),binN) ;
            %
            [ya,xa] = hist(enlistcell_new{isample}(:,5),nbin) ;
            %
            plot(xa,ya/sum(ya),'-',...
                'Color',[1 0 0] * (combo + 5 - 2) / (2*5 - 2),...
                'LineWidth',1,...
                'DisplayName',[short_name{isample} '(low threshold:' num2str(low_threshold) ')']) ;
            %%%
            y = enlistcell_new{isample}(:,6);
            minP = min(max(y(~isnan(y)),1)) ;
            maxP = max(y(~isnan(y))) ;
            if isempty(minP) || isempty(maxP)
                iplot = iplot + 1 ;
                continue;
            end
            nbin = logspace(log10(minP),log10(maxP),binN) ;
            %
            [ya,xa] = hist(y(~isnan(y)),nbin) ;
            %
            plot(xa,ya/sum(ya),'-',...
                'Color',[0 0 1] * (combo + 5 - 2) / (2*5 - 2),...
                'LineWidth',1,...
                'DisplayName',[short_name{isample} '(high threshold:' num2str(high_threshold) ')']) ;

            lh=legend('show');
            set(lh,'location','northwest');
            lh.FontSize = 8;

            xlim([1e3 1e7]) ;        
    %         ylabel('Probability') ;
            set(gca,'XTickLabel',[]) ;
            iplot = iplot + 1 ;
%             set(gca,'XTickLabelMode','auto') ;
%             xlabel('Spot intensity/pixel') ;
        end
    end
    set(gca,'XTickLabelMode','auto') ;
    xlabel('Cell intensity') ;
    
%     savefig(gcf, [datafolder '\fit\1D_spot_intensity.fig']);
%     saveas(gcf, [datafolder '\fit\1D_spot_intensity.png']);
end
    savefig(gcf, [pwd '\cluster_rec\1D_cell_intensity.fig']);
    saveas(gcf, [pwd '\cluster_rec\1D_cell_intensity.png']);
>>>>>>> 88dd4f4ac7eb6f5ae758552b04c47df241263be0
