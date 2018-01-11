function dprime_unity( stat_struct, selected_current )

    SCH_drugoff_dprimes = [];
    SCH_drugon_dprimes  = [];
    SCH_drugoff_avgs = [];
    SCH_drugon_avgs = [];
    
    SKF_drugoff_dprimes = [];
    SKF_drugon_dprimes  = [];
    SKF_drugoff_avgs = [];
    SKF_drugon_avgs = [];
    
    for i = 1:length(stat_struct )
        
        if isfield(stat_struct(i).attend, 'dmat')
        
            disp(stat_struct(i).filename)
            
            if length(stat_struct(i).attend) > 1

                for j = 1:length(stat_struct(i).attend)
                    if (stat_struct(i).attend(j).current == selected_current)
                        dprime_subset = get_relAttIn_dprimes( stat_struct(i).attend(j).dmat );
                        curr_Off_dprimes = dprime_subset(:,1);
                        curr_On_dprimes  = dprime_subset(:,2);
                        curr_Off_avgs = dprime_subset(:,3);
                        curr_On_avgs = dprime_subset(:,4);
                    else
                        curr_Off_dprimes = [];
                        curr_On_dprimes = [];
                        curr_Off_avgs = [];
                        curr_On_avgs = [];
                    end
                end
            
            else
                if (stat_struct(i).attend.current == selected_current)
                    dprime_subset = get_relAttIn_dprimes( stat_struct(i).attend.dmat );
                    curr_Off_dprimes = dprime_subset(:,1);
                    curr_On_dprimes  = dprime_subset(:,2);
                    curr_Off_avgs = dprime_subset(:,3);
                    curr_On_avgs = dprime_subset(:,4);                   
                else
                    curr_Off_dprimes = [];
                    curr_On_dprimes = [];
                    curr_Off_avgs = [];
                    curr_On_avgs = [];
                end        
            end

            
            
            if strcmp(stat_struct(i).drug, 'SCH23390')
                SCH_drugoff_dprimes = [SCH_drugoff_dprimes curr_Off_dprimes'];
                SCH_drugon_dprimes  = [SCH_drugon_dprimes  curr_On_dprimes' ];
                SCH_drugoff_avgs = [SCH_drugoff_avgs curr_Off_avgs'];
                SCH_drugon_avgs  = [SCH_drugon_avgs curr_On_avgs'];               
            elseif strcmp(stat_struct(i).drug, 'SKF81297')
                SKF_drugoff_dprimes = [SKF_drugoff_dprimes curr_Off_dprimes'];
                SKF_drugon_dprimes  = [SKF_drugon_dprimes  curr_On_dprimes' ];
                SKF_drugoff_avgs = [SKF_drugoff_avgs curr_Off_avgs'];
                SKF_drugon_avgs  = [SKF_drugon_avgs curr_On_avgs'];
            else
                disp('Drug not recognized.')
            end
        end

    end
    
    % Get Rank Sum p-value
    if ~ isempty( SCH_drugoff_dprimes)
    SCH_rs = ranksum( SCH_drugoff_dprimes, SCH_drugon_dprimes );
    else
        SCH_rs = [];
    end

    
    SCH_drugoff_dprimes_FIXED = SCH_drugoff_dprimes(~isnan(SCH_drugoff_dprimes));
    SCH_drugon_dprimes_FIXED = SCH_drugon_dprimes(~isnan(SCH_drugon_dprimes));
    SCH_drugon_dprimes_FIXED = SCH_drugon_dprimes_FIXED(~isinf(SCH_drugon_dprimes_FIXED));
    
    [h1 p1] = ttest(SCH_drugoff_dprimes_FIXED );
    m1 = median(SCH_drugoff_dprimes_FIXED);
    %n1 = round(length(SCH_drugoff_dprimes_FIXED) / 3);
    [h2 p2] = ttest( SCH_drugon_dprimes_FIXED );
    m2 = median(SCH_drugon_dprimes_FIXED);
    %n2 = round(length(SCH_drugon_dprimes_FIXED) / 3);
    disp( {'SCH drug off ttest:' num2str(p1) 'SCH drug off median: ' num2str(m1) });%  'SCH Off Num Cells: ' num2str(n1) });
    disp( {'SCH drug on ttest: ' num2str(p2) 'SCH drug on median: ' num2str(m2) });%'SCH On Num Cells: ' num2str(n2) } );
    
    
    if ~ isempty( SKF_drugoff_dprimes)
    SKF_rs = ranksum( SKF_drugoff_dprimes, SKF_drugon_dprimes );
    else
        SKF_rs = [];
    end
    
    SKF_drugoff_dprimes_FIXED = SKF_drugoff_dprimes(~isnan(SKF_drugoff_dprimes));
    SKF_drugon_dprimes_FIXED = SKF_drugon_dprimes(~isnan(SKF_drugon_dprimes));
    SKF_drugon_dprimes_FIXED = SKF_drugon_dprimes_FIXED(~isinf(SKF_drugon_dprimes_FIXED));
    
    [h1 p1] = ttest(SKF_drugoff_dprimes_FIXED );
    m1 = median(SKF_drugoff_dprimes_FIXED);
    %n1 = round(length(SKF_drugoff_dprimes_FIXED) / 3);
    [h2 p2] = ttest( SKF_drugon_dprimes_FIXED );
    m2 = median(SKF_drugon_dprimes_FIXED);
    %n2 = round(length(SKF_drugoff_dprimes_FIXED) / 3);
    disp( {'SKF drug off ttest:' num2str(p1) 'SKF drug off median: ' num2str(m1)});% 'SKF Off Num Cells: ' num2str(n1) } );
    disp( {'SKF drug on ttest: ' num2str(p2) 'SKF drug on median: ' num2str(m2)});% 'SKF On Num Cells: ' num2str(n1)} );
    
    
    figure('units','normalized','position',[.1 .1 .3 .42]);
    %%% SCH %%%
    %subplot(2,2,1);
    hold on;
    plot( SCH_drugoff_dprimes, SCH_drugon_dprimes, 'ok', 'MarkerFaceColor', 'red' );
    xlabel('Control D''', 'FontSize', 16, 'FontWeight', 'bold'); ylabel( 'D1R Antagonist D''', 'FontSize', 16, 'FontWeight', 'bold' );
    xlim([-2 2]);
    ylim([-2 2]);
    
    ax = gca;
    ax.XTick = [-2 -1 0 1 2];
    ax.YTick = ax.XTick;
    set(gca,'FontSize',12, 'FontWeight', 'bold');
    line( [-3 3], [-3 3], 'Color', 'black');
    
    text(-1.5, 1.5, ['p = ' num2str(round(SCH_rs, 4))], 'FontSize', 16, 'FontWeight', 'bold');
    hold off;
    title(strcat( 'D1R Antagonist',{' '} , num2str(selected_current), 'nA' ), 'FontSize', 18, 'Color', 'r');
    
    % Get Histogram of distances from line of unity
%     distance = (SCH_drugoff_dprimes - SCH_drugon_dprimes);
% 
%     valid_vals1 = distance(~isnan(distance));
%     valid_vals2 = valid_vals1(~isinf(valid_vals1));
%     mean_dist = mean( valid_vals2 );
%     
%     subplot(2,2,2);
%     h = histogram(distance, 'FaceColor', 'red');
%     height = max(h.Values);
%     hold on;
%     line( [0 0], [0, round(height,-1)], 'LineStyle', '--', 'Color', 'black');
%     plot( mean_dist, height, 'v', 'MarkerFaceColor', 'red', 'MarkerSize', 8);
%     xlabel( 'Control - D1R Antagonist', 'FontSize', 16, 'FontWeight', 'bold' );
%     ylabel( 'Count', 'FontSize', 16, 'FontWeight', 'bold' );
%     ax = gca;
%     %ax.XTick = [-2 -1 0 1 2];
%     %ax.YTick = [0 10 20 30 40 50];
%     set(gca,'FontSize',12, 'FontWeight', 'bold');

%    suptitle(gen_title(stat_struct));
%     hold off;
    
    %%% SKF %%%
    figure('units','normalized','position',[.1 .1 .3 .42]);   % subplot(2,2,3);
    hold on;
    plot( SKF_drugoff_dprimes, SKF_drugon_dprimes, 'ok', 'MarkerFaceColor', 'blue' );
    xlabel('Control D''', 'FontSize', 16, 'FontWeight', 'bold'); ylabel( 'D1R Agonist D''', 'FontSize', 16, 'FontWeight', 'bold' );
    xlim([-2 2]);
    ylim([-2 2]);
    
    ax = gca;
    ax.XTick = [-2 -1 0 1 2];
    ax.YTick = ax.XTick;
    set(gca,'FontSize',12, 'FontWeight', 'bold');
    line( [-3 3], [-3 3], 'Color', 'black');
    
    text(-1.5, 1.5, ['p = ' num2str(round(SKF_rs, 4))], 'FontSize', 16, 'FontWeight', 'bold');
      
    hold off; 
    title(strcat( 'D1R Agonist',{' '} , num2str(selected_current), 'nA' ), 'FontSize', 18, 'Color', 'b');
    
    % Get Histogram of distances from line of unity
%     distance = (SKF_drugoff_dprimes - SKF_drugon_dprimes);
% 
%     valid_vals1 = distance(~isnan(distance));
%     valid_vals2 = valid_vals1(~isinf(valid_vals1));
%     mean_dist = mean( valid_vals2 );
%     
%     subplot(2,2,4);
%     h = histogram(distance, 'FaceColor', 'blue');
%     height = max(h.Values);
%     hold on;
%     line( [0 0], [0, round(height,-1)], 'LineStyle', '--', 'Color', 'black');
%     plot( mean_dist, height, 'v', 'MarkerFaceColor', 'red', 'MarkerSize', 8);
%     xlabel( 'Control - Drug', 'FontSize', 16, 'FontWeight', 'bold' );
%     ylabel( 'Count', 'FontSize', 16, 'FontWeight', 'bold' );
%     ax = gca;
%     %ax.XTick = [-2 -1 0 1 2];
%     %ax.YTick = [0 10 20 30 40 50];
%     set(gca,'FontSize',12, 'FontWeight', 'bold');
% 
%     suptitle(gen_title(stat_struct, selected_current));
%     hold off;
    
    
    
    %%% BAR PLOT %%%
    
%     figure()
%     
%     subplot(1,2,1);
%     
    badidx_SCHoff = find(isnan(SCH_drugoff_avgs));
    badidx_SCHon = find(isnan(SCH_drugon_avgs));
    badidx_SKFoff = find(isnan(SKF_drugoff_avgs));
    badidx_SKFon = find(isnan(SKF_drugon_avgs));
    if length(badidx_SCHoff)
        SCH_drugoff_avgs(badidx_SCHoff) = [];
        SCH_drugon_avgs(badidx_SCHoff) = [];
    end
    if length(badidx_SCHon)
        SCH_drugoff_avgs(badidx_SCHon) = [];
        SCH_drugon_avgs(badidx_SCHon) = [];
    end
    if length(badidx_SKFoff)
        SKF_drugoff_avgs(badidx_SKFoff) = [];
        SKF_drugon_avgs(badidx_SKFoff) = [];
    end
    if length(badidx_SKFon)
        SKF_drugoff_avgs(badidx_SKFon) = [];
        SKF_drugon_avgs(badidx_SKFon) = [];
    end
    
    SCH_drug_rs = ranksum( SCH_drugoff_avgs, SCH_drugon_avgs );
    
%     SCH_off_ste = std(SCH_drugoff_avgs) / (sqrt(length(SCH_drugoff_avgs)));
%     SCH_on_ste = std(SCH_drugon_avgs) / (sqrt(length(SCH_drugon_avgs)));
%    
%     bar( [mean(SCH_drugoff_avgs) mean(SCH_drugon_avgs)], 'k' );
%     hold on;
%     errorbar(1:2, [mean(SCH_drugoff_avgs) mean(SCH_drugon_avgs)], [SCH_off_ste SCH_on_ste], '.k', 'LineWidth', 2)
%     name = {'Drug Off';'Drug On'};
%     set(gca,'xticklabel',name)
%     ylabel( 'Firing Rate (Hz)', 'FontSize', 16, 'FontWeight', 'bold' );
%     set(gca,'FontSize',12, 'FontWeight', 'bold');
%     text(1.5, (max(ylim) - 1), ['p = ' num2str(round(SCH_drug_rs, 6))], 'FontSize', 16, 'FontWeight', 'bold');
% 
%     hold off;

    SKF_drug_rs = ranksum( SKF_drugoff_avgs, SKF_drugon_avgs );

%     SKF_off_ste = std(SKF_drugoff_avgs) / (sqrt(length(SKF_drugoff_avgs)));
%     SKF_on_ste = std(SKF_drugon_avgs) / (sqrt(length(SKF_drugon_avgs)));
%     
%     subplot(1,2,2);
%     
%     bar( [mean(SKF_drugoff_avgs) mean(SKF_drugon_avgs)], 'b' );
%     hold on;
%     errorbar(1:2, [mean(SKF_drugoff_avgs) mean(SKF_drugon_avgs)], [SKF_off_ste SKF_on_ste], '.k', 'LineWidth', 2)
%     name = {'Drug Off';'Drug On'};
%     set(gca,'xticklabel',name)
%     ylabel( 'Firing Rate (Hz)', 'FontSize', 16, 'FontWeight', 'bold' );
%     set(gca,'FontSize',12, 'FontWeight', 'bold');
%     
%     text(1.5, (max(ylim) - 1), ['p = ' num2str(round(SKF_drug_rs, 6))], 'FontSize', 16, 'FontWeight', 'bold');
%     hold off;
    
    
        
    %%%%%% ATTEND FIGURE %%%%%%
    figure('units','normalized','position',[.1 .1 .3 .42]);
    
    %%% SCH %%%
    %subplot(2,2,1);
    hold on;
    plot( SCH_drugoff_avgs, SCH_drugon_avgs, 'ok', 'MarkerFaceColor', 'red' );
    xlabel({'Control FR' '(spikes/s)'}, 'FontSize', 16, 'FontWeight', 'bold'); ylabel( {'D1R Antagonist FR' '(spikes/s)'}, 'FontSize', 16, 'FontWeight', 'bold' );
    xlim([0 50]);
    ylim([0 50]);
    
    ax = gca;
    ax.XTick = [0 10 20 30 40 50];
    ax.YTick = ax.XTick;
    set(gca,'FontSize',12, 'FontWeight', 'bold');
    line( [-50 50], [-50 50], 'Color', 'black');
    
    text(5, 45, ['p = ' num2str(round(SCH_drug_rs, 2))], 'FontSize', 16, 'FontWeight', 'bold');
    hold off;
    title(strcat( 'D1R Antagonist',{' '} , num2str(selected_current), 'nA' ), 'FontSize', 18, 'Color', 'r');
    
    % Get Histogram of distances from line of unity
%     distance = (SCH_drugoff_avgs - SCH_drugon_avgs);
% 
%     valid_vals1 = distance(~isnan(distance));
%     valid_vals2 = valid_vals1(~isinf(valid_vals1));
%     mean_dist = mean( valid_vals2 );
%     
%     subplot(2,2,2);
%     h = histogram(distance, 'FaceColor', 'red');
%     height = max(h.Values);
%     hold on;
%     line( [0 0], [0, round(height,-1)], 'LineStyle', '--', 'Color', 'black');
%     plot( mean_dist, height, 'v', 'MarkerFaceColor', 'black', 'MarkerSize', 8);
%     xlabel( 'Control - D1R Antagonist', 'FontSize', 16, 'FontWeight', 'bold' );
%     ylabel( 'Count', 'FontSize', 16, 'FontWeight', 'bold' );
%     ax = gca;
%     %ax.XTick = [-2 -1 0 1 2];
%     %ax.YTick = [0 10 20 30 40 50];
%     set(gca,'FontSize',12, 'FontWeight', 'bold');
% 
%     %suptitle(gen_title(stat_struct));
%     hold off;

    %%% SKF %%%
    figure('units','normalized','position',[.1 .1 .3 .42]); % subplot(2,2,3);
    hold on;
    plot( SKF_drugoff_avgs, SKF_drugon_avgs, 'ok', 'MarkerFaceColor', 'blue' );
    xlabel({'Control FR' '(spikes/s)'}, 'FontSize', 16, 'FontWeight', 'bold'); ylabel( {'D1R Agonist FR' '(spikes/s)'}, 'FontSize', 16, 'FontWeight', 'bold' );
    xlim([0 50]);
    ylim([0 50]);
    
    ax = gca;
    ax.XTick = [0 10 20 30 40 50];
    ax.YTick = ax.XTick;
    set(gca,'FontSize',12, 'FontWeight', 'bold');
    line( [-50 50], [-50 50], 'Color', 'black');
    
    text(5, 45, ['p = ' num2str(round(SKF_drug_rs, 2))], 'FontSize', 16, 'FontWeight', 'bold');
    hold off;
    title(strcat( 'D1R Agonist',{' '} , num2str(selected_current), 'nA' ), 'FontSize', 18, 'Color', 'b');
    
    
    % Get Histogram of distances from line of unity
%     distance = (SKF_drugoff_avgs - SKF_drugon_avgs);
% 
%     valid_vals1 = distance(~isnan(distance));
%     valid_vals2 = valid_vals1(~isinf(valid_vals1));
%     mean_dist = mean( valid_vals2 );
%     
%     subplot(2,2,4);
%     h = histogram(distance, 'FaceColor', 'blue');
%     height = max(h.Values);
%     hold on;
%     line( [0 0], [0, round(height,-1)], 'LineStyle', '--', 'Color', 'black');
%     plot( mean_dist, height, 'v', 'MarkerFaceColor', 'red', 'MarkerSize', 8);
%     xlabel( 'Control - D1R Agonist', 'FontSize', 16, 'FontWeight', 'bold' );
%     ylabel( 'Count', 'FontSize', 16, 'FontWeight', 'bold' );
%     ax = gca;
%     %ax.XTick = [-2 -1 0 1 2];
%     %ax.YTick = [0 10 20 30 40 50];
%     set(gca,'FontSize',12, 'FontWeight', 'bold');
% 
%     suptitle(gen_title(stat_struct, selected_current));
%     hold off;
        
end

% This function gets the subset of d's which are likely to be the 'attend
% in' d' condition and its two closest neighbours. One point is therefore
% being dropped.
function dprime_subset = get_relAttIn_dprimes( dmat )

    Offs = dmat(:,2);
    Ons  = dmat(:,3);
    OffAvgs = dmat(:,6);
    OnAvgs = dmat(:,7);
    
    maxOff_idx = find(Offs == max(Offs));
    if maxOff_idx == 8
        idxs = [7 8 1];
    elseif maxOff_idx == 1
        idxs = [8 1 2];
    else
        idxs = [ (maxOff_idx -1) maxOff_idx (maxOff_idx +1) ];
    end
    
    OffVals = Offs(idxs);
    OnVals  = Ons( idxs );
    OffAvgVals = OffAvgs(idxs);
    OnAvgVals = OnAvgs(idxs);
    
    dprime_subset = [OffVals, OnVals, OffAvgVals, OnAvgVals];

end


% Figures out what title to give this plot
function titlestring = gen_title(stat_struct, selected_current)
    %drug = stat_struct(1).drug;
    %current = num2str(stat_struct(1).attend.current);
    current = num2str(selected_current);
    titlestring = strcat( {'Late-Cue Modulation at'}, {' '}, {current}, {'nA'} );
end