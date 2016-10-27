function dprime_unity( stat_struct )

    drugoff_dprimes = [];
    drugon_dprimes  = [];
    
    for i = 1:length(stat_struct )
        
        dprime_subset = get_relAttIn_dprimes( stat_struct(i).attend.dmat );
        curr_Off_dprimes = dprime_subset(:,1);
        curr_On_dprimes  = dprime_subset(:,2);
        
        drugoff_dprimes = [drugoff_dprimes curr_Off_dprimes'];
        drugon_dprimes  = [drugon_dprimes  curr_On_dprimes' ];

    end
    
    % Get Rank Sum p-value
    rs = ranksum( drugoff_dprimes, drugon_dprimes );
    
    figure();
    subplot(1,2,1);
    plot( drugoff_dprimes, drugon_dprimes, 'ok', 'MarkerFaceColor', 'black' );
    xlabel('Drug Off D''', 'FontSize', 16, 'FontWeight', 'bold'); ylabel( 'Drug On D''', 'FontSize', 16, 'FontWeight', 'bold' );
    xlim([-2 2]);
    ylim([-2 2]);
    
    ax = gca;
    ax.XTick = [-2 -1 0 1 2];
    ax.YTick = ax.XTick;
    set(gca,'FontSize',12, 'FontWeight', 'bold');
    line( [-3 3], [-3 3], 'Color', 'black');
    
    text(-1.5, 1.5, ['p = ' num2str(round(rs, 3))], 'FontSize', 16, 'FontWeight', 'bold');
    
    
    % Get Histogram of distances from line of unity
    distance = (drugoff_dprimes - drugon_dprimes) ./ sqrt(2);

    valid_vals1 = distance(~isnan(distance));
    valid_vals2 = valid_vals1(~isinf(valid_vals1));
    mean_dist = mean( valid_vals2 );
    
    subplot(1,2,2);
    h = histogram(distance, 'FaceColor', 'black');
    height = max(h.Values);
    hold on;
    line( [0 0], [0, round(height,-1)], 'LineStyle', '--', 'Color', 'black');
    plot( mean_dist, height, 'v', 'MarkerFaceColor', 'red', 'MarkerSize', 8);
    xlabel( 'Distance from Line of Unity', 'FontSize', 16, 'FontWeight', 'bold' );
    ylabel( 'Count', 'FontSize', 16, 'FontWeight', 'bold' );
    ax = gca;
    ax.XTick = [-2 -1 0 1 2];
    ax.YTick = [0 10 20 30 40 50];
    set(gca,'FontSize',12, 'FontWeight', 'bold');

    suptitle(gen_title(stat_struct));
    hold off;

end

% This function gets the subset of d's which are likely to be the 'attend
% in' d' condition and its two closest neighbours. One point is therefore
% being dropped.
function dprime_subset = get_relAttIn_dprimes( dmat )

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
    
    OffVals = Offs(idxs);
    OnVals  = Ons( idxs );
    
    dprime_subset = [OffVals, OnVals];

end


% Figures out what title to give this plot
function titlestring = gen_title(stat_struct)
    drug = stat_struct(1).drug;
    current = num2str(stat_struct(1).attend.current);
    titlestring = strcat( {'10mM '}, drug, {' '}, {current}, {'nA'} );
end