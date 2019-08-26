function tmp_get_anova_counts( mfs, drug, current, paradigm )


    if nargin < 1, load master_file_struct; mfs = master_file_struct; end
    
    
    counts = [0 0 0];
    total = 0;
    
    % Get Substruct of Specified Drug
    substruct = mfs.session( find( strcmp( {mfs.session.drug}, drug ) ) );
    
    %%% TMP - JOSE ONLY CODE %%%
%     tmp_list = strfind( {substruct.event_file}, 'Jose' );
%     tmp_list( cellfun(@isempty, tmp_list) ) = {0};
%     tmp_list = cell2mat( tmp_list );
%     substruct = substruct( tmp_list ~= 0 );
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Loop through substruct and strip out actual data - not filtered for
    % current or paradigm yet
    
    stripped_struct = []; fix_means = []; vis_means = [];
    for i = 1:length(substruct)
        for j = 1:length( substruct(i).vis_signif )
            
            % If chosen paradigm is not among paradigms, skip this unit
            paradigms = substruct(i).paradigms{j};
            if isempty( strcmp( paradigms, paradigm ) )
                continue
            end
            
            if isempty( substruct(i).vis_signif{j} )
                continue
            end
           
            currents = [substruct(i).vis_signif{j}.current];
            for k = 1:length( currents )
                if currents(k) == current
                    
                    pvals = substruct(i).vis_signif{1,j}(k).ps;
                    total = total+1;
                    
                    fix_means = [fix_means, substruct(i).vis_signif{1,j}(k).average_fix_fr ];
                    vis_means = [vis_means, substruct(i).vis_signif{1,j}(k).average_vis_fr ];
                    
                    for l = 1:3
                        if pvals(l) < 0.05
                            counts(l) = counts(l)+1;
                        end
                    end
                end
            end
        end
    end
                        
    figure('units','normalized','position',[.1 .1 0.4 0.2]);
    m_color = get_marker_color( drug );

    bar( counts ./ total *100, 'FaceColor', m_color)
    
    xt_labs = {'Fix vs Vis', 'Position', 'Fix/Vis * Position'};

    ylabel( '% of Neurons' );
    TickLabel_FontSize = 12;
    set( gca, 'YTick', [0 10 20 30 40 50 60], 'XTick', [1, 2, 3], 'XTickLabel', xt_labs, 'FontSize', TickLabel_FontSize,  ... 
        'FontWeight', 'Bold' ); xlim([0 4]); ylim( [0 60] ); box( gca, 'off');

    text(3.5, 5, ['N = ' num2str(total)], 'FontSize', 16, 'FontWeight', 'bold');
    title(strcat( drug,{' '} , num2str(current), 'nA' ), 'FontSize', 18, 'Color', m_color);
    
    figure('units','normalized','position',[.1 .1 .3 .42]);
    hold on; 
    plot( fix_means, vis_means, 'ok');
%    xl = xlim; yl = ylim;
%    max_axis_vals = max(xl, yl);
    max_axis_vals = [0 60]
    xlim(max_axis_vals); ylim(max_axis_vals);
    line( max_axis_vals, max_axis_vals, 'Color', 'black');

end

function m_color = get_marker_color( drug )

    if strcmp( drug, 'SCH23390' )
        m_color = 'red';
    elseif strcmp( drug, 'SKF81297' )
        m_color = 'blue';
    end

end