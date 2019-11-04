function epoch_anova_summary_fig = epoch_anova_barchart( mfs, drug, current, paradigm )

    if nargin < 1, load master_file_struct; mfs = master_file_struct; end
        
    
    if strcmp( paradigm, 'Attention' )
        epoch_list =  {'attend_fixation_stats', 'attend_visual_stats', 'attend_attend_stats', ...
                      'attend_blank_stats', 'attend_post_blank_stats', 'attend_reward_stats'};
                  
    elseif strcmp( paradigm, 'WM' )
        epoch_list =  {'wm_fixation_stats', 'wm_visual_stats', 'wm_delay_stats', ...
                     'wm_response_stats', 'wm_reward_stats'};
    end
    
   anova_struct = struct;     
    
    % Get Substruct of Specified Drug
    substruct = mfs.session( find( strcmp( {mfs.session.drug}, drug ) ) );
    
    % Loop through substruct and strip out actual data
    stripped_struct = [];
    for i = 1:length(substruct)
        for j = 1:length( substruct(i).stats)
            
            
            % If chosen paradigm is not among paradigms, skip this unit
            paradigms = substruct(i).paradigms{j};
            if isempty( strcmp( paradigms, paradigm ) ), continue, end
            
            if isempty( substruct(i).stats ), continue, end
            if isempty( substruct(i).stats{1,j}), continue, end
            
            % If there are no attend stats, skip this unit
            if  strcmp( paradigm, 'Attention' )
                if ~ isfield( substruct(i).stats{1,j}, 'attend_fixation_stats' ), continue, end
            elseif strcmp( paradigm, 'WM' )
                if ~ isfield( substruct(i).stats{1,j}, 'wm_fixation_stats' ), continue, end
            end
            
            %currents = [substruct(i).attend_stats{j}.current];
            currents = [substruct(i).currents{j}];
            currents = currents(2:end);
            for k = 1:length( currents )
            
                % Only grab data if current is the desired current
                if currents(k) == current
                   
                    % Loop through epoch list
                    for l = 1: length( epoch_list )

                        stats = substruct(i).stats{j};
                        factors = [stats.([epoch_list{l}])(k).anova_mat.tbl(2:end-2,1)];
                        p_vals = [stats.([epoch_list{l}])(k).anova_mat.tbl{2:end-2,7}];
                        
                        anova_struct.([epoch_list{l}]).factors = factors;
                        
                        if ~ isfield( anova_struct.([epoch_list{l}]), 'p_vals' )
                            anova_struct.([epoch_list{l}]).p_vals = p_vals; 
                        else
                            anova_struct.([epoch_list{l}]).p_vals = vertcat(anova_struct.([epoch_list{l}]).p_vals, p_vals ); 
                        end
                        
                    end
                end
            end
        end
    end

    % Loop through epochs and make bar plot of signif vs nonsignif for the
    % different anova factors
    
   epoch_anova_summary_fig = figure();
   hold on;
   
   fig_width = 1500;
   fig_height = 300;
   set(epoch_anova_summary_fig, 'Position', [10 10 fig_width 300]);
   
   num_epochs = length(epoch_list);
   
   for i = 1:num_epochs
      
       subplot( 1, num_epochs, i );
       
       signif_logical_mat = anova_struct.([epoch_list{i}]).p_vals <= 0.05 ;
       signif    = sum(signif_logical_mat,1);
       nonsignif = sum(~signif_logical_mat,1);
       
       factors = anova_struct.([epoch_list{i}]).factors;
       
       x = 1:length( factors );
       vals = [signif; nonsignif];
       b = bar( x, vals );
       
       xticklabels(factors);
       set(gca,'XTickLabelRotation',45)
       b(1).FaceColor = 'r';
       b(2).FaceColor = 'k';
       
       % Add epoch titles
       if strcmp(paradigm, 'Attention' )
           epoch_text = {'Fixation', 'Visual', 'Cue', 'Blank', 'Post-Blank', 'Reward'};
       elseif strcmp(paradigm, 'WM' )
           epoch_text = {'Fixation', 'Visual', 'Delay', 'Response', 'Reward'};
       end
   
       ylim_vals = ylim;
       new_ymax = ylim_vals(2) * 1.2;
       ylim( [0 new_ymax] );
       text(1, new_ymax*0.92, epoch_text{i}, 'FontSize', 14, 'FontWeight', 'bold');
  
       ylabel( 'Number of Units' )
   end

   %tightfig( epoch_anova_summary_fig );
   sgtitle( strcat( drug, {' '}, num2str(current), 'nA' ), 'FontWeight', 'bold' );
     
   %legend( {'p <= 0.05', 'p >= 0.05'}, 'Location', 'east' );
    legend( {'p <= 0.05', 'p >= 0.05'}, 'Position', [0.15 0.1, 0.05, 0.1] );

   
   save_name = strcat( 'tmp_figs/population_figs/', 'Anova_Summary', '_', paradigm,'_', drug, '_', num2str(current), 'nA', '.png' );
   print( epoch_anova_summary_fig, save_name, '-dpng'); % .png and .fig also posisble. % May Want to add paradigm to this eventually
     
   
end