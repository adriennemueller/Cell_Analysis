function summary_fig = gen_summary_fig( stats, currents )

    % 
    summary_fig = figure();
    

    % ATTEND
    
%     % If there are attend stats, create a 1x6 figure for each current
%     
%     if isfield( stats, 'attend_fixation_stats' )
%         
%         att_fig = figure();
%         
%         att_fields = {'attend_fixation_stats', 'attend_visual_stats', 'attend_attend_stats', ...
%                       'attend_blank_stats', 'attend_post_blank_stats', 'attend_reward_stats'};
%         
%          % For each current
%          currents = [stats.attend_fixation_stats.current];
%          for j = 1: length( currents )
% 
%              % Fixation plot
%              
%              % 
%              
%             for k = 1:6
%                 ctrl_means = [stats.([att_fields{k}]).control_summ_stats(j,:).avg_fr];
%                 ctrl_stes  = [stats.([att_fields{k}]).control_summ_stats(j,:).std_err];
%                 drug_means = [stats.([att_fields{k}]).drug_summ_stats(j,:).avg_fr];
%                 drug_stes  = [stats.([att_fields{k}]).drug_summ_stats(j,:).std_err];
%                 
%                 figure();
%                 errorbar( ctrl_means, ctrl_stes, 'k.' );
%                 hold on;
%                 errorbar( drug_means, drug_stes, 'r.' );
%                 
%                 
%                 % Add ANOVA RESULT
%             end
%             
%          end
%     end
         
    % WM
    % If there are wm stats, create a 1x5 figure for each current
    if isfield( stats, 'wm_fixation_stats' )
        
        wm_fig = figure(); hold on;
        wm_fields = {'wm_fixation_stats', 'wm_visual_stats', 'wm_delay_stats', ...
                     'wm_response_stats', 'wm_reward_stats'};

        % For each current
        currents = [stats.wm_fixation_stats.current];
        for j = 1: length( currents )

            % Fixation plot
            ctrl_stats = stats.wm_fixation_stats.control_summ_stats(j,:);
            ctrl_mean = gen_grand_av( [ctrl_stats.avg_fr], [ctrl_stats.num_trials] );
            ctrl_ste  = gen_grand_av( [ctrl_stats.std_err], [ctrl_stats.num_trials] );
             
            drug_stats = stats.wm_fixation_stats.drug_summ_stats(j,:);
            drug_mean = gen_grand_av( [drug_stats.avg_fr], [drug_stats.num_trials] );
            drug_ste  = gen_grand_av( [drug_stats.std_err], [drug_stats.num_trials] );
            
            subplot( 1, 5, 1)
            errorbar( ctrl_mean, ctrl_ste, 'k.' );
            hold on;
            errorbar( drug_mean, drug_ste, 'r.' );
            
            % Add anova results text
            stats_string = get_stats_string( stats.wm_fixation_stats.anova_mat );
            [text_x_loc, text_y_loc] = get_text_location( xlim, ylim );
            text( text_x_loc, text_y_loc, stats_string );
            
            for k = 2:5 % All Epochs after fixation
                ctrl_means = [stats.([wm_fields{k}]).control_summ_stats(j,:).avg_fr];
                ctrl_stes  = [stats.([wm_fields{k}]).control_summ_stats(j,:).std_err];
                drug_means = [stats.([wm_fields{k}]).drug_summ_stats(j,:).avg_fr];
                drug_stes  = [stats.([wm_fields{k}]).drug_summ_stats(j,:).std_err];
                
                subplot( 1, 5, k )
                errorbar( ctrl_means, ctrl_stes, 'k.' );
                hold on;
                errorbar( drug_means, drug_stes, 'r.' );

                % Add anova results text
                stats_string = get_stats_string( stats.([wm_fields{k}]).anova_mat );
                [text_x_loc, text_y_loc] = get_text_location( xlim, ylim );
                text( text_x_loc, text_y_loc, stats_string );
                
            end                
        end
    end
end

function [x_pos, y_pos] = get_text_location( xlimits, ylimits )
    x_pos = xlimits(1) + 0.2;
    y_pos = ylimits(2) - range(ylimits / 10);
end

function stats_string = get_stats_string( stats_struct )
    
    stats_string = [];

    factors_list = stats_struct.tbl( 2:end-2, 1 ); % Factors Column 
    p_vals       = stats_struct.tbl( 2:end-2, 7 ); % P values column
    
    for i = 1:length( factors_list )
        stats_string{i} = strcat( factors_list(i), {': '}, num2str( round(p_vals{i}, 3) ) );
    end
    
    stats_string = string(stats_string);
end        



function grand_av = gen_grand_av( vals, ns )

    total_n = sum(ns);
    weights = ns ./ total_n;

    grand_av = nansum( vals .* weights ); 

end
