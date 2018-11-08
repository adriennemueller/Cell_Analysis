function plot_contrast_dprimes()

    load( 'master_file_struct', 'master_file_struct' );
    mfs = master_file_struct;

    contrast_dprime_struct = struct;
    
    % Get set of attend d's for each unit, for each contrast. Should have a
    % cell array? of X cell, each with 8 dir by 4 contrast. - One for
    % Control, One for drug
    contrast_substruct = get_contrast_substruct( mfs );
   
    [ctrl_mat_full, drug_mat_full] = reshape_contrast_direc( contrast_substruct );
    
    
    % Plot what happens to the d's for each unit. (Four subplots - one for
    % each direction.
    for i = 1:length( contrast_substruct )
        
        % Take only 4 directions
        ctrl_mat = ctrl_mat_full( :, [7,8,1,2], i ); % ROWS ARE CONTRASTS, COLUMNS ARE DIRECTIONS
        drug_mat = drug_mat_full( :, [7,8,1,2], i ); % ROWS ARE CONTRASTS, COLUMNS ARE DIRECTIONS

       figure(i); hold on;
        set(gcf, 'Units', 'centimeters' );
        set(gcf, 'Position', [30, 30, 40, 10]);
         TickLabel_FontSize = 12;
        
        for j = 1:4 % 4 directions
            subplot( 1, 4, j )
            hold on;
            plot( ctrl_mat( :, j), 'k-');%j = each direction
            plot( drug_mat( :, j), 'r-'); %j = each direction
            xlabel( 'Contrast' );
            ylabel( 'Attentional Modulation (d'')' );
            set( gca, 'YTick', [-3 -2 -1 0 1 2 3], 'XTick', [1, 2, 3, 4], 'XTickLabel', {'10%', '15%', '20%', '25%'}, 'FontSize', TickLabel_FontSize,  ...
            'FontWeight', 'Bold' ); xlim([0 5]); ylim( [-3.5 3.5] ); box( gca, 'off');
            %legend( 'Control', 'SCH23390 50nA', 'Location', 'northwest');
            hold off;
        end
        tightfig( gcf );
        hold off;
    end
        
    % Make a summary plot showing how average d's (average within cells -
    % across directions - and then average across cells) change with drug -
    % for each contrast (4x2 bars = 4 contrasts, drug on, drug off );
    
    figure();
    set(gcf, 'Units', 'centimeters' );
    set(gcf, 'Position', [30, 30, 40, 10]);
        
    for i = 1:4 % NUM directions
        ctrl_mean_across_cells = nanmean( ctrl_mat_full( :, i, : ), 3 );
        ctrl_std_across_cells  = nanstd( ctrl_mat_full( :, i, : ), 0, 3 );
        
        drug_mean_across_cells = nanmean( drug_mat_full( :, i, : ), 3 );
        drug_std_across_cells  = nanstd( drug_mat_full( :, i, : ), 0, 3 );
        
        subplot( 1,4,i)
        hold on;
        plot( ctrl_mean_across_cells, 'k-' );
        plot( drug_mean_across_cells, 'r-' );
        errorbar( [1, 2, 3, 4], ctrl_mean_across_cells, ctrl_std_across_cells ./ sqrt(size(ctrl_mat_full, 3)), 'k' );
        errorbar( [1, 2, 3, 4], drug_mean_across_cells, drug_std_across_cells ./ sqrt(size(drug_mat_full, 3)), 'r' );
        
        xlabel( 'Contrast' );
        ylabel( 'Attentional Modulation (d'')' );
        set( gca, 'YTick', [-3 -2 -1 0 1 2 3], 'XTick', [1, 2, 3, 4], 'XTickLabel', {'10%', '15%', '20%', '25%'}, 'FontSize', TickLabel_FontSize,  ... 
            'FontWeight', 'Bold' ); xlim([0 5]); ylim( [-3.5 3.5] ); box( gca, 'off');
        hold off;
    end
    
    
    %%% Max Direction Only
    figure();
    set(gcf, 'Units', 'centimeters' );
    set(gcf, 'Position', [30, 30, 15, 15]);
        
    % Find direction with max control d' at 25% contrast
    
    numcells = size( ctrl_mat_full, 3 );
    
    for i = 1:numcells
        dprimes_at_25 = ctrl_mat_full( 4,:,i );
        max_dprimeat25_idx = find(dprimes_at_25 ==max(dprimes_at_25));
        
        ctrl_maxed25_dprime(i,:) = ctrl_mat_full( :, max_dprimeat25_idx, i );
        drug_maxed25_dprime(i,:) = drug_mat_full( :, max_dprimeat25_idx, i );
    end
        
    ctrl_mean25_across_cells = nanmean( ctrl_maxed25_dprime, 1 );
    ctrl_std25_across_cells  = nanstd( ctrl_maxed25_dprime, 1 );
        
    drug_mean25_across_cells = nanmean( drug_maxed25_dprime, 1 );
    drug_std25_across_cells  = nanstd( drug_maxed25_dprime, 1 );
        
    rs_list = [];
    for i = 1:4
        rs_val = ranksum( ctrl_maxed25_dprime(:,i), drug_maxed25_dprime(:,i)  );
        rs_list(i) = rs_val;
    end

        
    hold on;
    plot( ctrl_mean25_across_cells, 'k-' );
    plot( drug_mean25_across_cells, 'r-' );
    errorbar( [1, 2, 3, 4], ctrl_mean25_across_cells, ctrl_std25_across_cells ./ sqrt(length( ctrl_maxed25_dprime) ), 'k' );
    errorbar( [1, 2, 3, 4], drug_mean25_across_cells, drug_std25_across_cells ./ sqrt(length( drug_maxed25_dprime)), 'r' );
 
    
    
    xlabel( 'Contrast' );
    ylabel( 'Attentional Modulation (d'')' );
    set( gca, 'YTick', [0 1 2 3], 'XTick', [1, 2, 3, 4], 'XTickLabel', {'10%', '15%', '20%', '25%'}, 'FontSize', TickLabel_FontSize,  ... 
        'FontWeight', 'Bold' ); xlim([0 5]); ylim( [-0.5 4] ); box( gca, 'off');
    hold off;
    
end

function [ctrl_mat, drug_mat] = reshape_contrast_direc( contrast_substruct )
    for u = 1:length( contrast_substruct ) 
        for i = 1:length( contrast_substruct(u,:) )
            ctrl_dprimes = get_direc_dprimes( [contrast_substruct(u,i).ctrl_attend_stats.direction],  [contrast_substruct(u,i).ctrl_attend_stats.dprime_val] );
            drug_dprimes = get_direc_dprimes( [contrast_substruct(u,i).drug_attend_stats.direction],  [contrast_substruct(u,i).drug_attend_stats.dprime_val] );

            ctrl_dprimes( ~isfinite( ctrl_dprimes )) = NaN;
            drug_dprimes( ~isfinite( drug_dprimes )) = NaN;
            
            ctrl_mat( i,:,u ) = ctrl_dprimes;
            drug_mat( i,:,u ) = drug_dprimes;
        end
    end
end

function direc_vals = get_direc_dprimes( direcs, dprimes )

    full_direcs = [0 45 90 135 180 225 270 315];

    direc_vals = NaN( 1, 8 );
    
    for i = 1:length(direcs) 
        direc_idx = find( full_direcs ==  direcs(i) );
        direc_vals( direc_idx ) = dprimes( i );
    end

end


function c_substruct = get_contrast_substruct( mfs )

    mfs_indices = find( ~cellfun(@isempty, {mfs.session.attendContrast}) );
    mfs_substruct = mfs.session( mfs_indices );
    
    for i = 1:length(mfs_substruct)

        
        [numcells,numcontrasts] = size( (mfs_substruct(i).attendContrast) );
        
        for j = 1:numcells
            
            for k = 1:numcontrasts
            
                curr_struct(k).contrast = mfs_substruct(i).attendContrast{j,k};

                within_cell_idx = find( [mfs_substruct(i).attendContrast_stats{j,k}.current] == 50 );

                attend_struct = mfs_substruct(i).attendContrast_stats{j,k};

                curr_struct(k).ctrl_attend_stats = attend_struct(within_cell_idx).control_dmat.dmat;
                curr_struct(k).drug_attend_stats = attend_struct(within_cell_idx).drug_dmat.dmat;
            
            end
                 if ~exist( 'c_substruct' ), c_substruct = curr_struct;
                 else, c_substruct = vertcat(c_substruct, curr_struct);
                 end
                 
        end
    end

end