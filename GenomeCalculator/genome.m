<<<<<<< HEAD
function [geno, gene, gene_time] = genome(t,tau,X,OC,L,gene_name,show_fig,C,D)
    if nargin == 7
        exp_C = [42,42,54,52,91,105];
        exp_D = [33,33,34,53,78,60];
        exp_tau = [0,60,96,137,190,259];
        
        f = fittype( @(b, c, x) log(1 + b.^(c*x-c*exp_tau(2)))./log(b) + exp_C(1));
        fitobject_C = fit(transpose(exp_tau), transpose(exp_C), f, 'Startpoint', [5,.3], 'Lower', [1.1, .1], 'Upper', [100, 5]);
        f = fittype( @(b, c, x) log(1 + b.^(c*x-c*exp_tau(2)))./log(b) + exp_D(1));
        fitobject_D = fit(transpose(exp_tau), transpose(exp_D), f, 'Startpoint', [5,.3], 'Lower', [1.1, .1], 'Upper', [100, 5]);
    
        C = fitobject_C(tau);
        D = fitobject_D(tau);
    elseif nargin == 9
        %
    else
        error('not enough input arguements')
    end
    
    t1 = t(t<=tau);
    t2 = t(t>tau) - tau;
    
    [geno1, gene1, gene_time1] = genome2(t1,tau,X,OC,L,gene_name,show_fig,C,D);
    if ~isempty(t2) || length(t2) == 1
        [geno2, gene2, gene_time2] = genome2(t2,tau,X,OC,L,gene_name,show_fig,C,D);
        geno = [geno1, geno2*2];
        gene = [gene1, gene2*2];
        gene_time = [gene_time1];
    end
    geno = geno1;
    gene = gene1;
    gene_time = gene_time1;
end


function [geno, gene, gene_time] = genome2(t,tau,X,OC,L,gene_name,show_fig,C,D)
    if nargin == 7
        exp_C = [42,42,54,52,91,105];
        exp_D = [33,33,34,53,78,60];
        exp_tau = [0,60,96,137,190,259];
        
        f = fittype( @(b, c, x) log(1 + b.^(c*x-c*exp_tau(2)))./log(b) + exp_C(1));
        fitobject_C = fit(transpose(exp_tau), transpose(exp_C), f, 'Startpoint', [5,.3], 'Lower', [1.1, .1], 'Upper', [100, 5]);
        f = fittype( @(b, c, x) log(1 + b.^(c*x-c*exp_tau(2)))./log(b) + exp_D(1));
        fitobject_D = fit(transpose(exp_tau), transpose(exp_D), f, 'Startpoint', [5,.3], 'Lower', [1.1, .1], 'Upper', [100, 5]);
    
        C = fitobject_C(tau);
        D = fitobject_D(tau);
    elseif nargin == 9
        %
    else
        error('not enough input arguements')
    end

    %% Genome Algorithm
    generations_crossed = floor((C+D)/tau) + 1;
    
    hanging_previous = (C+D) - floor((C+D)/tau)*tau;
    hanging_next = (C-min(C, hanging_previous)) - floor((C-min(C, hanging_previous))/tau)*tau;
    completely_covered = (C - min(C, hanging_previous) - hanging_next)/tau;
    
    %%% % Linear growth: t/C + 1
    %base_slope = 1/C;
    slope = sum(2.^[max(generations_crossed-2-completely_covered, 0):...
        generations_crossed-2] ./ C);
    overlap = min(tau - (hanging_previous + hanging_next), 0);
    if overlap == 0
        t_1 = hanging_next;
        t_2 = tau - hanging_previous;
        if 0 <= (generations_crossed-2-completely_covered)
            %Slope = [slope*(1+completely_covered), slope*completely_covered, slope*(1+completely_covered)];
            Slope = [slope, slope - 2^(generations_crossed-2-completely_covered)/C, ...
                slope - 2^(generations_crossed-2-completely_covered)/C + 2^(generations_crossed-1)/C];
            t_div = [min(hanging_next, t); min(max(0, tau-hanging_previous-hanging_next), max(0, t-hanging_next)); ...
                min(hanging_previous, max(0, t-(tau-hanging_previous)))];
            t_max = [min(hanging_next, tau); min(max(0, tau-hanging_previous-hanging_next), max(0, tau-hanging_next)); ...
                min(hanging_previous, max(0, tau-(tau-hanging_previous)))];
        else
            Slope = [slope, slope, slope + 2^(generations_crossed-1)/C];
            t_div = [min(hanging_next, t); min(max(0, tau-hanging_previous-hanging_next), max(0, t-hanging_next)); ...
                min(hanging_previous-D, max(0, t-(tau-hanging_previous)))];
            t_max = [min(hanging_next, tau); min(max(0, tau-hanging_previous-hanging_next), max(0, tau-hanging_next)); ...
                min(hanging_previous-D, max(0, tau-(tau-hanging_previous)))];
        end
    else
        %Slope = [slope*(1+completely_covered), slope*(0+completely_covered), slope*(1+completely_covered)];
        Slope = [slope, slope + 2^(generations_crossed-1)/C, ...
            slope - 2^(generations_crossed-2-completely_covered)/C + 2^(generations_crossed-1)/C];
        t_1 = tau-hanging_previous;
        t_2 = hanging_next;
        t_div = [min(tau-hanging_previous, t); -max(overlap, min(0, -t+(hanging_next+overlap))); ...
            min(hanging_previous, max(0, t-(tau-hanging_previous-overlap)))];
        t_max = [min(tau-hanging_previous, tau); -max(overlap, min(0, -tau+(hanging_next+overlap))); ...
            min(hanging_previous, max(0, tau-(tau-hanging_previous-overlap)))];
    end
    %disp(Slope .* C);
    %disp([0, t_1, t_2]);
    geno = (Slope * t_div) + (Slope * t_max);
    
%     %%% % Exponential growth: 2^(t/C)
%     
%     geno = 0;
%     %time_in_current_cycles = t+hanging_previous+(i-1)*tau
%     for i = 1:(completely_covered+1)
%         cycle_time = min(t+(i-1)*tau+min(hanging_previous, C), C);
%         geno = geno + (2.^(cycle_time ./ C) - 0);
%         cycle_time = min(0+(i-1)*tau+min(hanging_previous, C), C);
%         geno = geno - (2.^(cycle_time ./ C) - 0);
%     end
%     %time_in_next_cycle = min(max(0, t-time_until_hanging_previous), C)
%     cycle_time = min(max(0, t-(tau-hanging_previous)), C);
%     geno = geno + (2.^(cycle_time ./ C) - 1);
%     %add base amount
%     geno = geno + 1;
    

    %% Gene Algorithm
    speed = (L/2)/C;
    distance = min(abs(OC - X), abs(L-OC+X));
    time = distance/speed;
    C_prime = C - time+D;
    hanging_previous = (C_prime) - floor((C_prime)/tau)*tau;
    hanging_next = (C_prime-hanging_previous) - floor((C_prime-hanging_previous)/tau)*tau;
    completely_covered = (C_prime - hanging_previous - hanging_next)/tau;
    gene_time = tau-hanging_previous;
    gene = 2.^(completely_covered + ceil(min(1, max(0, t-tau+hanging_previous))));
    
    
    %% Plot of cell cycle
    if show_fig
        figure();
        for i = 1:generations_crossed
            D_time = [max(0, (generations_crossed - (i-1)) * tau - D),...
                (generations_crossed - (i-1)) * tau];
            C_time = [max(0, (generations_crossed - (i-1)) * tau - D - C),...
                max(0, (generations_crossed - (i-1)) * tau - D)];

    %         hold on;
    %         rectangle('Position', [D_time(1), i-1, D_time(end)-D_time(1), 1], 'FaceColor', 'b');
    %         hold on;
    %         rectangle('Position', [C_time(1), i-1, C_time(end)-C_time(1), 1], 'FaceColor', 'r');
            hold on;
            rectangle('Position', [D_time(1), (generations_crossed-i), D_time(end)-D_time(1), 1], 'FaceColor', 'b');
            hold on;
            rectangle('Position', [C_time(1), (generations_crossed-i), C_time(end)-C_time(1), 1], 'FaceColor', 'r');

            hold on;
            plot([D_time(end), D_time(end)], [0, generations_crossed], '--', ...
                'Color', 'k', 'LineWidth', 2, 'HandleVisibility', 'off');

            hold on;
            gene_time = (generations_crossed - (i-1)) * tau - C_prime;
    %         plot([gene_time, gene_time], [i-1, i], '--', 'Color', 'g', 'LineWidth', 2, 'HandleVisibility', 'off');
            plot([gene_time, gene_time], [(generations_crossed-i), (generations_crossed-i)+1], ...
                '--', 'Color', 'g', 'LineWidth', 2, 'HandleVisibility', 'off');
        end
        hold on;
        xlim([0, generations_crossed*tau]);
        xlabel('time (minutes)');
        %yticks([]);
        ylabel('Log_2(oriC in the cell)');
        yticks([0:1:100]);
        title(['C = ',num2str(C), 'min, D = ',num2str(D), 'min, tau = ',num2str(tau), 'min. K12-MG1655, ', gene_name]);
        bar([nan],[nan], 'b', 'DisplayName', 'D period');
        bar([nan],[nan], 'r', 'DisplayName', 'C period');
        plot([nan], [nan], '--', 'Color', 'g', 'DisplayName', 'Gene duplication');
        plot([nan], [nan], '--', 'Color', 'k', 'DisplayName', 'Cell division');
        legend();

        figure();
        subplot(3,1,1)
        for i = 1:generations_crossed
            D_time = [max(0, (generations_crossed - (i-1)) * tau - D),...
                (generations_crossed - (i-1)) * tau];
            C_time = [max(0, (generations_crossed - (i-1)) * tau - D - C),...
                max(0, (generations_crossed - (i-1)) * tau - D)];

    %         hold on;
    %         rectangle('Position', [D_time(1), i-1, D_time(end)-D_time(1), 1], 'FaceColor', 'b');
    %         hold on;
    %         rectangle('Position', [C_time(1), i-1, C_time(end)-C_time(1), 1], 'FaceColor', 'r');
            hold on;
            rectangle('Position', [D_time(1), (generations_crossed-i), D_time(end)-D_time(1), 1], 'FaceColor', 'b');
            hold on;
            rectangle('Position', [C_time(1), (generations_crossed-i), C_time(end)-C_time(1), 1], 'FaceColor', 'r');

            hold on;
            plot([D_time(end), D_time(end)], [0, generations_crossed], '--', ...
                'Color', 'k', 'LineWidth', 2, 'HandleVisibility', 'off');

            hold on;
            gene_time = (generations_crossed - (i-1)) * tau - C_prime;
    %         plot([gene_time, gene_time], [i-1, i], '--', 'Color', 'g', 'LineWidth', 2, 'HandleVisibility', 'off');
            plot([gene_time, gene_time], [(generations_crossed-i), (generations_crossed-i)+1], ...
                '--', 'Color', 'g', 'LineWidth', 2, 'HandleVisibility', 'off');
            if gene_time  > 0 && gene_time < tau
                gene_time_true = gene_time;
            end
        end
        hold on;
        xlim([0, tau]);
        xlabel('time (minutes)');
        %     yticks([]);
        ylabel('Log_2(oriC in the cell)');
        yticks([0:1:100]);
        title(['C = ',num2str(C), 'min, D = ',num2str(D), 'min, tau = ',num2str(tau), 'min. K12-MG1655, ', gene_name]);
        bar([nan],[nan], 'b', 'DisplayName', 'D period');
        bar([nan],[nan], 'r', 'DisplayName', 'C period');
        plot([nan], [nan], '--', 'Color', 'g', 'DisplayName', 'Gene duplication');
        legend();

        subplot(3,1,2)
        plot(t, geno);
        xlabel('time (minutes)');
        %xlim([0 tau]);
        ylabel('Genome copies (t)');
        title(['Genome copies over time (initial copies = ', num2str(Slope * t_max), 'x one chromosome)']);

        subplot(3,1,3)
        plot(t, gene);
        xlabel('time (minutes)');
        %xlim([0 tau]);
        ylim([0 inf]);
        ylabel('Gene copies (t)');
        title(['Gene copy number, duplicates at t = ', num2str(gene_time_true), 'minutes']);
    end
end
=======
function [geno, gene, gene_time] = genome(t,tau,X,OC,L,gene_name,show_fig,C,D)
    if nargin == 7
        exp_C = [42,42,54,52,91,105];
        exp_D = [33,33,34,53,78,60];
        exp_tau = [0,60,96,137,190,259];
        
        f = fittype( @(b, c, x) log(1 + b.^(c*x-c*exp_tau(2)))./log(b) + exp_C(1));
        fitobject_C = fit(transpose(exp_tau), transpose(exp_C), f, 'Startpoint', [5,.3], 'Lower', [1.1, .1], 'Upper', [100, 5]);
        f = fittype( @(b, c, x) log(1 + b.^(c*x-c*exp_tau(2)))./log(b) + exp_D(1));
        fitobject_D = fit(transpose(exp_tau), transpose(exp_D), f, 'Startpoint', [5,.3], 'Lower', [1.1, .1], 'Upper', [100, 5]);
    
        C = fitobject_C(tau);
        D = fitobject_D(tau);
    elseif nargin == 9
        %
    else
        error('not enough input arguements')
    end
    
    t1 = t(t<=tau);
    t2 = t(t>tau) - tau;
    
    [geno1, gene1, gene_time1] = genome2(t1,tau,X,OC,L,gene_name,show_fig,C,D);
    if ~isempty(t2) || length(t2) == 1
        [geno2, gene2, gene_time2] = genome2(t2,tau,X,OC,L,gene_name,show_fig,C,D);
        geno = [geno1, geno2*2];
        gene = [gene1, gene2*2];
        gene_time = [gene_time1];
    end
    geno = geno1;
    gene = gene1;
    gene_time = gene_time1;
end


function [geno, gene, gene_time] = genome2(t,tau,X,OC,L,gene_name,show_fig,C,D)
    if nargin == 7
        exp_C = [42,42,54,52,91,105];
        exp_D = [33,33,34,53,78,60];
        exp_tau = [0,60,96,137,190,259];
        
        f = fittype( @(b, c, x) log(1 + b.^(c*x-c*exp_tau(2)))./log(b) + exp_C(1));
        fitobject_C = fit(transpose(exp_tau), transpose(exp_C), f, 'Startpoint', [5,.3], 'Lower', [1.1, .1], 'Upper', [100, 5]);
        f = fittype( @(b, c, x) log(1 + b.^(c*x-c*exp_tau(2)))./log(b) + exp_D(1));
        fitobject_D = fit(transpose(exp_tau), transpose(exp_D), f, 'Startpoint', [5,.3], 'Lower', [1.1, .1], 'Upper', [100, 5]);
    
        C = fitobject_C(tau);
        D = fitobject_D(tau);
    elseif nargin == 9
        %
    else
        error('not enough input arguements')
    end

    %% Genome Algorithm
    generations_crossed = floor((C+D)/tau) + 1;
    
    hanging_previous = (C+D) - floor((C+D)/tau)*tau;
    hanging_next = (C-min(C, hanging_previous)) - floor((C-min(C, hanging_previous))/tau)*tau;
    completely_covered = (C - min(C, hanging_previous) - hanging_next)/tau;
    
    %%% % Linear growth: t/C + 1
    %base_slope = 1/C;
    slope = sum(2.^[max(generations_crossed-2-completely_covered, 0):...
        generations_crossed-2] ./ C);
    overlap = min(tau - (hanging_previous + hanging_next), 0);
    if overlap == 0
        t_1 = hanging_next;
        t_2 = tau - hanging_previous;
        if 0 <= (generations_crossed-2-completely_covered)
            %Slope = [slope*(1+completely_covered), slope*completely_covered, slope*(1+completely_covered)];
            Slope = [slope, slope - 2^(generations_crossed-2-completely_covered)/C, ...
                slope - 2^(generations_crossed-2-completely_covered)/C + 2^(generations_crossed-1)/C];
            t_div = [min(hanging_next, t); min(max(0, tau-hanging_previous-hanging_next), max(0, t-hanging_next)); ...
                min(hanging_previous, max(0, t-(tau-hanging_previous)))];
            t_max = [min(hanging_next, tau); min(max(0, tau-hanging_previous-hanging_next), max(0, tau-hanging_next)); ...
                min(hanging_previous, max(0, tau-(tau-hanging_previous)))];
        else
            Slope = [slope, slope, slope + 2^(generations_crossed-1)/C];
            t_div = [min(hanging_next, t); min(max(0, tau-hanging_previous-hanging_next), max(0, t-hanging_next)); ...
                min(hanging_previous-D, max(0, t-(tau-hanging_previous)))];
            t_max = [min(hanging_next, tau); min(max(0, tau-hanging_previous-hanging_next), max(0, tau-hanging_next)); ...
                min(hanging_previous-D, max(0, tau-(tau-hanging_previous)))];
        end
    else
        %Slope = [slope*(1+completely_covered), slope*(0+completely_covered), slope*(1+completely_covered)];
        Slope = [slope, slope + 2^(generations_crossed-1)/C, ...
            slope - 2^(generations_crossed-2-completely_covered)/C + 2^(generations_crossed-1)/C];
        t_1 = tau-hanging_previous;
        t_2 = hanging_next;
        t_div = [min(tau-hanging_previous, t); -max(overlap, min(0, -t+(hanging_next+overlap))); ...
            min(hanging_previous, max(0, t-(tau-hanging_previous-overlap)))];
        t_max = [min(tau-hanging_previous, tau); -max(overlap, min(0, -tau+(hanging_next+overlap))); ...
            min(hanging_previous, max(0, tau-(tau-hanging_previous-overlap)))];
    end
    %disp(Slope .* C);
    %disp([0, t_1, t_2]);
    geno = (Slope * t_div) + (Slope * t_max);
    
%     %%% % Exponential growth: 2^(t/C)
%     
%     geno = 0;
%     %time_in_current_cycles = t+hanging_previous+(i-1)*tau
%     for i = 1:(completely_covered+1)
%         cycle_time = min(t+(i-1)*tau+min(hanging_previous, C), C);
%         geno = geno + (2.^(cycle_time ./ C) - 0);
%         cycle_time = min(0+(i-1)*tau+min(hanging_previous, C), C);
%         geno = geno - (2.^(cycle_time ./ C) - 0);
%     end
%     %time_in_next_cycle = min(max(0, t-time_until_hanging_previous), C)
%     cycle_time = min(max(0, t-(tau-hanging_previous)), C);
%     geno = geno + (2.^(cycle_time ./ C) - 1);
%     %add base amount
%     geno = geno + 1;
    

    %% Gene Algorithm
    speed = (L/2)/C;
    distance = min(abs(OC - X), abs(L-OC+X));
    time = distance/speed;
    C_prime = C - time+D;
    hanging_previous = (C_prime) - floor((C_prime)/tau)*tau;
    hanging_next = (C_prime-hanging_previous) - floor((C_prime-hanging_previous)/tau)*tau;
    completely_covered = (C_prime - hanging_previous - hanging_next)/tau;
    gene_time = tau-hanging_previous;
    gene = 2.^(completely_covered + ceil(min(1, max(0, t-tau+hanging_previous))));
    
    
    %% Plot of cell cycle
    if show_fig
        figure();
        for i = 1:generations_crossed
            D_time = [max(0, (generations_crossed - (i-1)) * tau - D),...
                (generations_crossed - (i-1)) * tau];
            C_time = [max(0, (generations_crossed - (i-1)) * tau - D - C),...
                max(0, (generations_crossed - (i-1)) * tau - D)];

    %         hold on;
    %         rectangle('Position', [D_time(1), i-1, D_time(end)-D_time(1), 1], 'FaceColor', 'b');
    %         hold on;
    %         rectangle('Position', [C_time(1), i-1, C_time(end)-C_time(1), 1], 'FaceColor', 'r');
            hold on;
            rectangle('Position', [D_time(1), (generations_crossed-i), D_time(end)-D_time(1), 1], 'FaceColor', 'b');
            hold on;
            rectangle('Position', [C_time(1), (generations_crossed-i), C_time(end)-C_time(1), 1], 'FaceColor', 'r');

            hold on;
            plot([D_time(end), D_time(end)], [0, generations_crossed], '--', ...
                'Color', 'k', 'LineWidth', 2, 'HandleVisibility', 'off');

            hold on;
            gene_time = (generations_crossed - (i-1)) * tau - C_prime;
    %         plot([gene_time, gene_time], [i-1, i], '--', 'Color', 'g', 'LineWidth', 2, 'HandleVisibility', 'off');
            plot([gene_time, gene_time], [(generations_crossed-i), (generations_crossed-i)+1], ...
                '--', 'Color', 'g', 'LineWidth', 2, 'HandleVisibility', 'off');
        end
        hold on;
        xlim([0, generations_crossed*tau]);
        xlabel('time (minutes)');
        %yticks([]);
        ylabel('Log_2(oriC in the cell)');
        yticks([0:1:100]);
        title(['C = ',num2str(C), 'min, D = ',num2str(D), 'min, tau = ',num2str(tau), 'min. K12-MG1655, ', gene_name]);
        bar([nan],[nan], 'b', 'DisplayName', 'D period');
        bar([nan],[nan], 'r', 'DisplayName', 'C period');
        plot([nan], [nan], '--', 'Color', 'g', 'DisplayName', 'Gene duplication');
        plot([nan], [nan], '--', 'Color', 'k', 'DisplayName', 'Cell division');
        legend();

        figure();
        subplot(3,1,1)
        for i = 1:generations_crossed
            D_time = [max(0, (generations_crossed - (i-1)) * tau - D),...
                (generations_crossed - (i-1)) * tau];
            C_time = [max(0, (generations_crossed - (i-1)) * tau - D - C),...
                max(0, (generations_crossed - (i-1)) * tau - D)];

    %         hold on;
    %         rectangle('Position', [D_time(1), i-1, D_time(end)-D_time(1), 1], 'FaceColor', 'b');
    %         hold on;
    %         rectangle('Position', [C_time(1), i-1, C_time(end)-C_time(1), 1], 'FaceColor', 'r');
            hold on;
            rectangle('Position', [D_time(1), (generations_crossed-i), D_time(end)-D_time(1), 1], 'FaceColor', 'b');
            hold on;
            rectangle('Position', [C_time(1), (generations_crossed-i), C_time(end)-C_time(1), 1], 'FaceColor', 'r');

            hold on;
            plot([D_time(end), D_time(end)], [0, generations_crossed], '--', ...
                'Color', 'k', 'LineWidth', 2, 'HandleVisibility', 'off');

            hold on;
            gene_time = (generations_crossed - (i-1)) * tau - C_prime;
    %         plot([gene_time, gene_time], [i-1, i], '--', 'Color', 'g', 'LineWidth', 2, 'HandleVisibility', 'off');
            plot([gene_time, gene_time], [(generations_crossed-i), (generations_crossed-i)+1], ...
                '--', 'Color', 'g', 'LineWidth', 2, 'HandleVisibility', 'off');
            if gene_time  > 0 && gene_time < tau
                gene_time_true = gene_time;
            end
        end
        hold on;
        xlim([0, tau]);
        xlabel('time (minutes)');
        %     yticks([]);
        ylabel('Log_2(oriC in the cell)');
        yticks([0:1:100]);
        title(['C = ',num2str(C), 'min, D = ',num2str(D), 'min, tau = ',num2str(tau), 'min. K12-MG1655, ', gene_name]);
        bar([nan],[nan], 'b', 'DisplayName', 'D period');
        bar([nan],[nan], 'r', 'DisplayName', 'C period');
        plot([nan], [nan], '--', 'Color', 'g', 'DisplayName', 'Gene duplication');
        legend();

        subplot(3,1,2)
        plot(t, geno);
        xlabel('time (minutes)');
        %xlim([0 tau]);
        ylabel('Genome copies (t)');
        title(['Genome copies over time (initial copies = ', num2str(Slope * t_max), 'x one chromosome)']);

        subplot(3,1,3)
        plot(t, gene);
        xlabel('time (minutes)');
        %xlim([0 tau]);
        ylim([0 inf]);
        ylabel('Gene copies (t)');
        title(['Gene copy number, duplicates at t = ', num2str(gene_time_true), 'minutes']);
    end
end
>>>>>>> 88dd4f4ac7eb6f5ae758552b04c47df241263be0
