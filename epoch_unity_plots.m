

function epoch_unity_fig = epoch_unity_plots( mfs, drug, current, paradigm, theta_invariant_flag )

    if nargin < 1, load master_file_struct; mfs = master_file_struct; end
        
    
    if strcmp( paradigm, 'Attention' )
        epoch_list =  {'attend_fixation_stats', 'attend_visual_stats', 'attend_attend_stats', ...
                      'attend_blank_stats', 'attend_post_blank_stats', 'attend_reward_stats'};
                  
    elseif strcmp( paradigm, 'WM' )
        epoch_list =  {'wm_fixation_stats', 'wm_visual_stats', 'wm_delay_stats', ...
                     'wm_response_stats', 'wm_reward_stats'};
    end
    
    ctrl_vals = [];
    drug_vals = [];
    curr_unit = 1;
        
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
                        
                        % Test if theta invariant
                        if theta_invariant_flag
                            % Discard epoch element if not theta_invariant
                            signif_theta = test_theta_invariance( stats.([epoch_list{l}])(k).anova_mat.tbl );
                            if signif_theta, continue, end
                        end
                        
                        ctrl_stats = [stats.([epoch_list{l}])(k).control_summ_stats];
                        ctrl_mean = gen_grand_av( [ctrl_stats.avg_fr], [ctrl_stats.num_trials] );
                        ctrl_vals(l,curr_unit) = ctrl_mean;
                        
                        drug_stats = [stats.([epoch_list{l}])(k).drug_summ_stats];
                        drug_mean = gen_grand_av( [drug_stats.avg_fr], [drug_stats.num_trials] );
                        drug_vals(l,curr_unit) = drug_mean;
                        
                    end
                    
                    curr_unit = curr_unit+1;
                end
            end
        end
    end
    
   epoch_unity_fig = figure();
   hold on;
   
   set(epoch_unity_fig, 'Position', [10 10 1500 200]);
    
   
   for i = 1:length(epoch_list)
      
       subplot( 1, length(epoch_list), i );
       
       plot(ctrl_vals(i,:), drug_vals(i,:), 'k.');
       
       % Fix dimensions and add unity line
       xlabel('Control FR', 'FontSize', 12, 'FontWeight', 'bold'); ylabel( [drug 'FR'], 'FontSize', 12, 'FontWeight', 'bold' );
       xlim([0 50]); ylim([0 50]);
       ax.XTick = [0 50];
       ax.YTick = ax.XTick;
       line( [-50 50], [-50 50], 'Color', 'black');
       p_x = l*2; p_y = 40;
       
       % Add P value
       rs_pval = signrank( ctrl_vals(i,:), drug_vals(i,:) );
       rs_str = get_pval_string(rs_pval);
       text(p_x, p_y, ['p ' rs_str], 'FontSize', 12, 'FontWeight', 'bold');

       
       % Add epoch titles
       if strcmp(paradigm, 'Attention' )
           epoch_text = {'Fixation', 'Visual', 'Cue', 'Blank', 'Post-Blank', 'Reward'};
       elseif strcmp(paradigm, 'WM' )
            epoch_text = {'Fixation', 'Visual', 'Delay', 'Response', 'Reward'};
       end
       text(p_x, p_y + 7, epoch_text{i}, 'FontSize', 14, 'FontWeight', 'bold');
   
   end
    
    
end


function signif_theta = test_theta_invariance( anova_tbl )


    theta_idx_cell = cellfun(@(x) length(strfind(x, 'theta')) > 0, anova_tbl(:,1), 'UniformOutput', false);
    theta_idxs = find( [theta_idx_cell{:}] == 1 );
    
    if isempty( theta_idxs )
            signif_theta = 0;
            return
    end

    p_vals = [anova_tbl{ theta_idxs, 7 }];    
    signif_theta = sum(find( p_vals <= 0.05 )); % Should consider making this smaller b/c multiple comparisons.
    
end

function grand_av = gen_grand_av( vals, ns )

    total_n = sum(ns);
    weights = ns ./ total_n;

    grand_av = nansum( vals .* weights ); 

end



function pval_str = get_pval_string( pval )
    if isempty(pval)
        pval_str = '= NaN';
    elseif pval >= 0.01
        pval_str = [ '= ' num2str(round(pval, 3))];
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
