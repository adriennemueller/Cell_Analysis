function dprime_unity_vis( stat_struct, selected_current )

    SCH_drugoff_dprimes = [];
    SCH_drugon_dprimes  = [];
    
    SKF_drugoff_dprimes = [];
    SKF_drugon_dprimes  = [];
    
    for i = 1:length(stat_struct )
        
        if isfield(stat_struct(i).attend, 'vis_pval')
        
            disp(stat_struct(i).filename)
            
            if length(stat_struct(i).attend) > 1
                for j = 1:length(stat_struct(i).attend)
                    if (stat_struct(i).attend(j).current == selected_current)  %&& (test_attend(stat_struct(i).attend(j).anova_mat))
                        %dprime_subset = get_dprime_subset( stat_struct(i).attend(j).vis_pval, stat_struct(i).attend(j).dmat );
                        dprime_subset = get_dprime_subset_single( stat_struct(i).attend(j).vis_pval, stat_struct(i).attend(j).dmat );
                        curr_Off_dprimes = dprime_subset(:,1);
                        curr_On_dprimes  = dprime_subset(:,2);
                    else
                        curr_Off_dprimes = [];
                        curr_On_dprimes = [];
                    end
                end
            
            else
                    if (stat_struct(i).attend.current == selected_current) %&& (test_attend(stat_struct(i).attend.anova_mat))
                        %dprime_subset = get_dprime_subset( stat_struct(i).attend.vis_pval, stat_struct(i).attend.dmat );
                        dprime_subset = get_dprime_subset_single( stat_struct(i).attend.vis_pval, stat_struct(i).attend.dmat );
                        curr_Off_dprimes = dprime_subset(:,1);
                        curr_On_dprimes  = dprime_subset(:,2);
                    else
                        curr_Off_dprimes = [];
                        curr_On_dprimes = [];
                    end
            
            end

            if strcmp(stat_struct(i).drug, 'SCH23390')
                SCH_drugoff_dprimes = [SCH_drugoff_dprimes curr_Off_dprimes'];
                SCH_drugon_dprimes  = [SCH_drugon_dprimes  curr_On_dprimes' ];
            elseif strcmp(stat_struct(i).drug, 'SKF81297')
                SKF_drugoff_dprimes = [SKF_drugoff_dprimes curr_Off_dprimes'];
                SKF_drugon_dprimes  = [SKF_drugon_dprimes  curr_On_dprimes' ];
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
    
    if ~ isempty( SKF_drugoff_dprimes)
        SKF_rs = ranksum( SKF_drugoff_dprimes, SKF_drugon_dprimes );
    else
        SKF_rs = [];
    end
    
    figure('units','normalized','position',[.1 .1 .3 .42]);
    %%% SCH %%%
    %subplot(2,2,1);
    hold on;
    plot( SCH_drugoff_dprimes, SCH_drugon_dprimes, 'ok', 'MarkerFaceColor', 'red' );
    xlabel({'Control FR' '(spikes/s)'}, 'FontSize', 16, 'FontWeight', 'bold'); ylabel( {'D1R Antagonist FR' '(spikes/s)'}, 'FontSize', 16, 'FontWeight', 'bold' );
    xlim([0 50]);
    ylim([0 50]);
    
    ax = gca;
    ax.XTick = [0 10 20 30 40 50];
    ax.YTick = ax.XTick;
    set(gca,'FontSize',12, 'FontWeight', 'bold');
    line( [-50 50], [-50 50], 'Color', 'black');
    
    SCH_rs_pval_str = get_pval_string( SCH_rs );
    text(5, 45, ['p ' SCH_rs_pval_str], 'FontSize', 16, 'FontWeight', 'bold');
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
%     plot( mean_dist, height, 'v', 'MarkerFaceColor', 'black', 'MarkerSize', 8);
%     xlabel( 'Drug Off - Drug On', 'FontSize', 16, 'FontWeight', 'bold' );
%     ylabel( 'Count', 'FontSize', 16, 'FontWeight', 'bold' );
%     ax = gca;
%     %ax.XTick = [-2 -1 0 1 2];
%     %ax.YTick = [0 10 20 30 40 50];
%     set(gca,'FontSize',12, 'FontWeight', 'bold');

    %suptitle(gen_title(stat_struct));
%     hold off;
    
    %%% SKF %%%
    figure('units','normalized','position',[.1 .1 .3 .42]);    % subplot(2,2,3);
    hold on;
    plot( SKF_drugoff_dprimes, SKF_drugon_dprimes, 'ok', 'MarkerFaceColor', 'blue' );
    xlabel({'Control FR' '(spikes/s)'}, 'FontSize', 16, 'FontWeight', 'bold'); ylabel( {'D1R Agonist FR' '(spikes/s)'}, 'FontSize', 16, 'FontWeight', 'bold' );
    xlim([0 50]);
    ylim([0 50]);
    
    ax = gca;
    ax.XTick = [0 10 20 30 40 50];
    ax.YTick = ax.XTick;
    set(gca,'FontSize',12, 'FontWeight', 'bold');
    line( [-50 50], [-50 50], 'Color', 'black');
    
    SKF_pval_str = get_pval_string( SKF_rs );
    text(5, 45, ['p ' SKF_pval_str], 'FontSize', 16, 'FontWeight', 'bold');
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
%     xlabel( 'Drug Off - Drug On', 'FontSize', 16, 'FontWeight', 'bold' );
%     ylabel( 'Count', 'FontSize', 16, 'FontWeight', 'bold' );
%     ax = gca;
%     %ax.XTick = [-2 -1 0 1 2];
%     %ax.YTick = [0 10 20 30 40 50];
%     set(gca,'FontSize',12, 'FontWeight', 'bold');
% 
%     suptitle(gen_title(stat_struct, selected_current));
%     hold off;
    

end

function pval_str = get_pval_string( pval )
    if isempty(pval)
        pval_str = '= NaN';
    elseif pval >= 0.01
        pval_str = [ '= ' num2str(round(pval, 2))];
    elseif (0.001 < pval) && (pval < 0.01)
        pval_str = ['< 0.01'];
    elseif pval < 0.000001
        pval_str = ['< 1*10-6' ];
    elseif pval < 0.001
        pval_str = ['< 0.001'];
    else
        pval_str = '= NaN';
    end
    % Make for smaller
end


function rslt = test_attend( anova_mat )
    attend_anova_idxs = [3, 5, 6, 7];
    pvals = anova_mat.p;
    
    rslt = pvals(attend_anova_idxs) <= 0.05;
    rslt = max(rslt);
    
end


% This function gets the subset of d's which are likely to be the 'attend
% in' d' condition and its two closest neighbours. One point is therefore
% being dropped.
function dprime_subset = get_dprime_subset_single( vispval, dmat )
    vis_Offs = vispval(:,4);
    vis_Ons  = vispval(:,5);

    Offs = dmat(:,2);
    Ons  = dmat(:,3);
    
    maxOff_idx = find(Offs == max(Offs));
    
    OffVal = vis_Offs(maxOff_idx);
    OnVal  = vis_Ons( maxOff_idx );

    dprime_subset = [OffVal, OnVal];

end


% This function gets the subset of d's which are likely to be the 'attend
% in' d' condition and its two closest neighbours. One point is therefore
% being dropped.
function dprime_subset = get_dprime_subset( vispval, dmat )


%%% If Taking All Four
%     Offs = vispval(:,4);
%     Ons  = vispval(:,5);
%     
%     idxs = 1:4;
%     
%     OffVals = Offs(idxs);
%     OnVals  = Ons( idxs );

    vis_Offs = vispval(:,4);
    vis_Ons  = vispval(:,5);

    Offs = dmat(:,2);
    Ons  = dmat(:,3);
    
    maxOff_idx = find(Offs == max(Offs));
    if maxOff_idx == 8
        idxs = [7 8 1];
    elseif maxOff_idx == 1
        idxs = [8 1 2];
    else
        idxs = [ (maxOff_idx -1) maxOff_idx (maxOff_idx +1) ];
    end
    
    OffVals = vis_Offs(idxs);
    OnVals  = vis_Ons( idxs );

    dprime_subset = [OffVals, OnVals];

end


% Figures out what title to give this plot
function titlestring = gen_title(stat_struct, selected_current)
    %drug = stat_struct(1).drug;
    current = num2str(selected_current);
    titlestring = strcat( {'Visual Response at'}, {' '}, {current}, {'nA'} );
end