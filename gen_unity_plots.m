
% unity_type = 'd_prime' or 'mean_fr'

function [control_vals, drug_vals, rs_pval] = gen_unity_plots( mfs, drug, current, paradigm, unity_type, signif_att_flag )

    if nargin < 1, load master_file_struct; mfs = master_file_struct; end
    
    
    % Get Substruct of Specified Drug
    substruct = mfs.session( find( strcmp( {mfs.session.drug}, drug ) ) );
    
    % Loop through substruct and strip out actual data - not filtered for
    % current or paradigm yet
    
    stripped_struct = [];
    for i = 1:length(substruct)
        for j = 1:length( substruct(i).attend_stats)
            
            % If chosen paradigm is not among paradigms, skip this unit
            paradigms = substruct(i).paradigms{j};
            if isempty( strcmp( paradigms, paradigm ) )
                continue
            end
            
            if isempty( substruct(i).attend_stats{j} )
                continue
            end
            
            
            currents = [substruct(i).attend_stats{j}.current];
            for k = 1:length( currents )
                if currents(k) == current
                    
                    break_me = 0;
                    if signif_att_flag
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
            
                    new_element.attend_stats          = tmp_element.attend_stats(k);
                    new_element.attend_visual_stats   = tmp_element.attend_visual_stats(k);
                    new_element.attend_fixation_stats = tmp_element.attend_fixation_stats(k);
                    
                    stripped_struct = [stripped_struct, new_element];
                end
            end
        end
    end
    
    
    %%% FOUR POINTS 
    % Loop through structs making a vector of drug on and drug off for
    % specified period and unity_type (d_prime or mean_fr )
    if strcmp( unity_type, 'd_prime' )
        four_pt_ctrl_vals = [];
        four_pt_drug_vals = [];

        for i = 1:length(stripped_struct)
        
            tmp_four_pt_ctrl_vals = [stripped_struct(i).attend_stats.control_dmat.dmat.dprime_val];
            tmp_four_pt_drug_vals = [stripped_struct(i).attend_stats.drug_dmat.dmat.dprime_val];
            
            four_pt_ctrl_vals = horzcat( four_pt_ctrl_vals, tmp_four_pt_ctrl_vals(1, [8,1,2]) );
            four_pt_drug_vals = horzcat( four_pt_drug_vals, tmp_four_pt_drug_vals(1, [8,1,2]) );
        end
        unity_plot( four_pt_ctrl_vals, four_pt_drug_vals, drug, current, 'd_prime' );
    end

    
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
        
            best_idx = find([stripped_struct(i).attend_visual_stats.control_summ_stats.avg_fr] == max([stripped_struct(i).attend_visual_stats.control_summ_stats.avg_fr]));
%             best_idx = best_idx(ismember( best_idx, [7,8,1,2] ));
%             best_idx = best_idx(ismember( best_idx, [8,1,2,3] ));
%             best_idx = best_idx(ismember( best_idx, [2,3,4,5] ));
%             best_idx = best_idx(ismember( best_idx, [3,4,5,6] ));
%             best_idx = best_idx(ismember( best_idx, [4,5,6,7] ));
%             best_idx = best_idx(ismember( best_idx, [5,6,7,8] ));
%             best_idx = best_idx(ismember( best_idx, [6,7,8,1] ));
                
            tmp_bestdir_ctrl_vals = [stripped_struct(i).attend_stats.control_dmat.dmat.dprime_val];
            tmp_bestdir_drug_vals = [stripped_struct(i).attend_stats.drug_dmat.dmat.dprime_val];
            
            bestdir_ctrl_vals = horzcat( bestdir_ctrl_vals, tmp_bestdir_ctrl_vals(1, best_idx) );
            bestdir_drug_vals = horzcat( bestdir_drug_vals, tmp_bestdir_drug_vals(1, best_idx) );
        end
        [control_vals, drug_vals, rs_pval] = unity_plot( bestdir_ctrl_vals, bestdir_drug_vals, drug, current, 'd_prime' );
    end
    
   
    %%% ATTEND WINDOW FIRING RATE
    
    
    
end