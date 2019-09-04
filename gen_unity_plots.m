
% Variables to plot include:
% unity_type = 'd_prime', 'mean_fr', 'mod_idx'
% position_type = 'best_vis', 'best_cue'
%   position_type determines which position (of the 4 or 8 positions) to
%   choose to plot.
%   'best_vis' will plot the position with the highest activity
%       during the visual, targets on, period.
%   'best_cue' will plot the position with the highest activity during the
%       cue period.
% signif_att_type = 'none', 'visVfix', 'attend'
%   signif_att_type determines whether to only plot those units with
%   anova-determined significant activity
%   'none' will plot all units
%   'visVfix' will only plot units with significant activity comparing the
%       visual period to the fixation period, with position also a factor.
%   'attend' will only plot units with significant attentional modulation
%       during the cue period (drug and position also included as factors)

function [control_vals, drug_vals, rs_pval] = gen_unity_plots( mfs, drug, current, paradigm, unity_type, position_type, signif_att_flag )

    if nargin < 1, load master_file_struct; mfs = master_file_struct; end
        
    % Get Substruct of Specified Drug
    substruct = mfs.session( find( strcmp( {mfs.session.drug}, drug ) ) );
    
    % Loop through substruct and strip out actual data
    stripped_struct = [];
    for i = 1:length(substruct)
        for j = 1:length( substruct(i).attend_stats)
            
            % If chosen paradigm is not among paradigms, skip this unit
            paradigms = substruct(i).paradigms{j};
            if isempty( strcmp( paradigms, paradigm ) ), continue, end
            
            % If there are no attend stats, skip this unit
            if isempty( substruct(i).attend_stats{j} ), continue, end
            
            currents = [substruct(i).attend_stats{j}.current];
            for k = 1:length( currents )
            
                % Only grab data if current is the desired current
                if currents(k) == current
                    
                    break_me = 0;
                    if strcmp( 'visVfix', signif_att_flag )
                        pvals = substruct(i).vis_signif{1,j}(k).ps; % Using significant visual period activity compared to fixation as criterion
                        if ~sum(find(pvals <= 0.05))
                            break_me = 1;
                        end 
                    elseif strcmp( 'attend', signif_att_flag )
                        pvals = substruct(i).attend_visual_stats{1,j}(k).anova_mat.tbl( [2,5,6,8] , 7 );
                        pvals = cell2mat(pvals);
                        if ~sum(find(pvals < 0.05))
                            break_me = 1;
                        end
                    end
                    if break_me, continue; end
                    
                    
                    tmp_element.attend_stats = substruct(i).attend_stats{j};
                    tmp_element.attend_visual_stats = substruct(i).attend_visual_stats{j};
                    tmp_element.attend_fixation_stats = substruct(i).attend_fixation_stats{j};
                    tmp_element.vis_signif = substruct(i).vis_signif{j};
                    
                    new_element.attend_stats          = tmp_element.attend_stats(k);
                    new_element.attend_visual_stats   = tmp_element.attend_visual_stats(k);
                    new_element.attend_fixation_stats = tmp_element.attend_fixation_stats(k);
                    new_element.vis_signif = tmp_element.vis_signif(k); %% CHECK THIS! WHY K?
                    
                    if strfind( substruct(i).event_file, 'Jose' ), animal = 'Jose'; else animal = 'Garfunkel'; end
                    new_element.animal = animal;
            
                    
                    stripped_struct = [stripped_struct, new_element];
                end
            end
        end
    end
    
    
    %%% FOUR POINTS 
    % Loop through structs making a vector of drug on and drug off for
%     % specified period and unity_type (d_prime or mean_fr )
%     if strcmp( unity_type, 'd_prime' )
%         four_pt_ctrl_vals = [];
%         four_pt_drug_vals = [];
% 
%         for i = 1:length(stripped_struct)
%         
%             tmp_four_pt_ctrl_vals = [stripped_struct(i).attend_stats.control_dmat.dmat.dprime_val];
%             tmp_four_pt_drug_vals = [stripped_struct(i).attend_stats.drug_dmat.dmat.dprime_val];
%             
%             four_pt_ctrl_vals = horzcat( four_pt_ctrl_vals, tmp_four_pt_ctrl_vals(1, [8,1,2]) );
%             four_pt_drug_vals = horzcat( four_pt_drug_vals, tmp_four_pt_drug_vals(1, [8,1,2]) );
%         end
%         unity_plot( four_pt_ctrl_vals, four_pt_drug_vals, drug, current, 'd_prime' );
%     end

    
    %%% FIXATION FIRING RATE
    if strcmp( unity_type, 'mean_fr' )
        bestdir_ctrl_vals = [];
        bestdir_drug_vals = [];

        for i = 1:length(stripped_struct)
        
            best_idx = find([stripped_struct(i).attend_visual_stats.control_summ_stats.avg_fr] == max([stripped_struct(i).attend_visual_stats.control_summ_stats.avg_fr]));
           %  best_idx = best_idx(ismember( best_idx, [7,8,1,2] ));
%           
            best_idx = best_idx(ismember( best_idx, [1,2,3,4] )); %%% THIS ONE


            tmp_bestdir_ctrl_vals = [stripped_struct(i).attend_fixation_stats.control_summ_stats.avg_fr];
            tmp_bestdir_drug_vals = [stripped_struct(i).attend_fixation_stats.drug_summ_stats.avg_fr];
            
            bestdir_ctrl_vals = horzcat( bestdir_ctrl_vals, tmp_bestdir_ctrl_vals(1, best_idx) );
            bestdir_drug_vals = horzcat( bestdir_drug_vals, tmp_bestdir_drug_vals(1, best_idx) );
        end
       [control_vals, drug_vals, rs_pval] = unity_plot( bestdir_ctrl_vals, bestdir_drug_vals, drug, current, 'mean_fr' );
     end
    
    
    %%% Best PreCue Direction
    if strcmp( unity_type, 'd_prime' )
        bestdir_ctrl_vals = [];
        bestdir_drug_vals = [];

        for i = 1:length(stripped_struct)
        
            best_idx = find(round([stripped_struct(i).attend_visual_stats.control_summ_stats.avg_fr], 4) == max(round([stripped_struct(i).attend_visual_stats.control_summ_stats.avg_fr],4)));
            best_idx = best_idx(ismember( best_idx, [1,2,7,8] )); %Because two idential; pick the more contralateral one. Could try 3 instead of 8
            
            
%             best_idx = best_idx(ismember( best_idx, [7,8,1,2] ));
%             best_idx = best_idx(ismember( best_idx, [8,1,2,3] ));
%             best_idx = best_idx(ismember( best_idx, [2,3,4,5] ));
%             best_idx = best_idx(ismember( best_idx, [3,4,5,6] ));
%             best_idx = best_idx(ismember( best_idx, [4,5,6,7] ));
%             best_idx = best_idx(ismember( best_idx, [5,6,7,8] ));
%             best_idx = best_idx(ismember( best_idx, [6,7,8,1] ));
             
            % ORIGINAL WITH D PRIME
            %tmp_bestdir_ctrl_vals = [stripped_struct(i).attend_stats.control_dmat.dmat.dprime_val];
            %tmp_bestdir_drug_vals = [stripped_struct(i).attend_stats.drug_dmat.dmat.dprime_val];
            
            %bestdir_ctrl_vals = horzcat( bestdir_ctrl_vals, tmp_bestdir_ctrl_vals(1, best_idx) );
            %bestdir_drug_vals = horzcat( bestdir_drug_vals, tmp_bestdir_drug_vals(1, best_idx) );
            
            % NEW WITH MOD. INDEX
            tmp_bestdir_ctrl_vals = [stripped_struct(i).attend_stats.control_dmat.dmat.mod_idx];
            tmp_bestdir_drug_vals = [stripped_struct(i).attend_stats.drug_dmat.dmat.mod_idx];
            
            bestdir_ctrl_vals = horzcat( bestdir_ctrl_vals, tmp_bestdir_ctrl_vals(1, best_idx) );
            bestdir_drug_vals = horzcat( bestdir_drug_vals, tmp_bestdir_drug_vals(1, best_idx) );
        end
        [control_vals, drug_vals, rs_pval] = unity_plot( bestdir_ctrl_vals, bestdir_drug_vals, drug, current, 'd_prime' );
    end
    
   
    %%% ATTEND WINDOW FIRING RATE
    
    
    
end