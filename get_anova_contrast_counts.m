function get_anova_contrast_counts( mfs, drug, current, paradigm )


    if nargin < 1, load master_file_struct; mfs = master_file_struct; end
    
    
    counts = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ];
    total = 0;
    
    % Get Substruct of Specified Drug
    substruct = mfs.session( find( strcmp( {mfs.session.drug}, drug ) ) );
    
    % Loop through substruct and strip out actual data - not filtered for
    % current or paradigm yet
    
    stripped_struct = [];
    for i = 1:length(substruct)
        for j = 1:length( substruct(i).attendContrast_stats)
            
            % If chosen paradigm is not among paradigms, skip this unit
            paradigms = substruct(i).paradigms{j};
            if isempty( strcmp( paradigms, paradigm ) )
                continue
            end
            
            if isempty( substruct(i).attendContrast_stats{j} )
                continue
            end
           
            currents = [substruct(i).attendContrast_stats{j}.current];
            for k = 1:length( currents )
                if currents(k) == current
                    
                    pvals = substruct(i).attendContrast_stats{1,j}(k).anova_mat.tbl( 2:16 , 7 );
                    pvals = cell2mat(pvals);
                    total = total+1;
                    
                    for l = 1:15
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
    
        
    xt_labs = {'Pos', 'Drug', 'Att Dir', 'Contrast', 'Pos * Drug', 'Pos * Att Dir', 'Pos Contrast', ...
               'Drug * Att Dir', 'Drug * Contrast', 'Att Dir * Contrast',  'Pos * Drug * Att Dir', 'Pos * Drug * Contrast', ...
               'Pos * Att Dir * Contrast', 'Drug * Att Dir * Contrast', 'Pos * Drug * Att Dir * Contrast'};

    ylabel( '% of Neurons' );
    TickLabel_FontSize = 12;
    set( gca, 'YTick', [0 25 50 75 100], 'XTick', [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15], 'XTickLabel', xt_labs, 'FontSize', TickLabel_FontSize,  ... 
        'FontWeight', 'Bold' ); xlim([0 16]); ylim( [0 100] ); box( gca, 'off');
    xtickangle(45)

    text(15.5, 5, ['N = ' num2str(total)], 'FontSize', 16, 'FontWeight', 'bold');
    title(strcat( drug,{' '} , num2str(current), 'nA' ), 'FontSize', 18, 'Color', m_color);
    

end

function m_color = get_marker_color( drug )

    if strcmp( drug, 'SCH23390' )
        m_color = 'red';
    elseif strcmp( drug, 'SKF81297' )
        m_color = 'blue';
    end

end