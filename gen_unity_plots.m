
% Variables to plot include:
% unity_type = 'd_prime', 'mean_fr', 'mod_idx'
% position_type = 'best_vis', 'best_cue'
%   position_type determines which position (of the 4 or 8 positions) to
%   choose to plot.
%   'best_vis' will plot the position with the highest activity
%       during the visual, targets on, period.
%   'best_cue' will plot the position with the highest activity during the
%       cue period.


function [control_vals, drug_vals, rs_pval] = gen_unity_plots( mfs, drug, current, paradigm, unity_type, position_type )

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
                    
                    signif_visfix = 0;
                    signif_attend = 0;
                    % Using significant visual period activity compared to fixation as criterion
                    visfix_pvals = substruct(i).vis_signif{1,j}(k).ps; 
                    if sum(find(visfix_pvals <= 0.05)), signif_visfix = 1; end 
                    
                    %%% CHECK WHAT THIS IS ACTUALLY CHECKING
                    attend_pvals = substruct(i).attend_visual_stats{1,j}(k).anova_mat.tbl( [2,5,6,8] , 7 );
                    attend_pvals = cell2mat(attend_pvals);
                    if sum(find(attend_pvals < 0.05)), signif_attend = 1; end
                    
                    tmp_element.attend_stats = substruct(i).attend_stats{j};
                    tmp_element.attend_visual_stats = substruct(i).attend_visual_stats{j};
                    tmp_element.attend_fixation_stats = substruct(i).attend_fixation_stats{j};
                    tmp_element.vis_signif = substruct(i).vis_signif{j};
                    
                    new_element.attend_stats          = tmp_element.attend_stats(k);
                    new_element.attend_visual_stats   = tmp_element.attend_visual_stats(k);
                    new_element.attend_fixation_stats = tmp_element.attend_fixation_stats(k);
                    new_element.vis_signif = tmp_element.vis_signif(k); %% CHECK THIS! WHY K?
                    
                    if strfind( substruct(i).event_file, 'Jose' ), animal = 'Jose'; else, animal = 'Garfunkel'; end
                    new_element.animal = animal;
                    new_element.signif_visfix = signif_visfix;
                    new_element.signif_attend = signif_attend;
            
                    stripped_struct = [stripped_struct, new_element];
                end
            end
        end
    end
    
    % Setup Control and Drug Vectors
    bestdir_ctrl_vals = [];
    bestdir_drug_vals = [];

    for i = 1:length(stripped_struct)
    
        if strcmp( position_type,  'best_vis' )
            best_idx = find(round([stripped_struct(i).attend_visual_stats.control_summ_stats.avg_fr], 4) == max(round([stripped_struct(i).attend_visual_stats.control_summ_stats.avg_fr],4)));
        elseif strcmp( position_type, 'best_cue' )
            best_idx = find([stripped_struct(i).attend_visual_stats.control_summ_stats.avg_fr] == max([stripped_struct(i).attend_visual_stats.control_summ_stats.avg_fr]));
        end
        
        best_idx = best_idx(ismember( best_idx, [7,8,1,2] )); 
        %best_idx = best_idx(ismember( best_idx, [8,1,2,3] )); %Because could choose either 90º or 270º position 
        
        % If two identical values for some reason, just pick one.
        if length(best_idx), best_idx = best_idx(1); end
        
        % Plot Firing Rate
        if strcmp( unity_type, 'mean_fr' )
            tmp_bestdir_ctrl_vals = [stripped_struct(i).attend_fixation_stats.control_summ_stats.avg_fr];
            tmp_bestdir_drug_vals = [stripped_struct(i).attend_fixation_stats.drug_summ_stats.avg_fr];
        % Plot d prime
        elseif strcmp( unity_type, 'd_prime' )
            tmp_bestdir_ctrl_vals = [stripped_struct(i).attend_stats.control_dmat.dmat.dprime_val];
            tmp_bestdir_drug_vals = [stripped_struct(i).attend_stats.drug_dmat.dmat.dprime_val];
        % Plot modulation index
        elseif strcmp( unity_type, 'mod_idx' )
            tmp_bestdir_ctrl_vals = [stripped_struct(i).attend_stats.control_dmat.dmat.mod_idx];
            tmp_bestdir_drug_vals = [stripped_struct(i).attend_stats.drug_dmat.dmat.mod_idx];
        end
            
        bestdir_ctrl_vals = horzcat( bestdir_ctrl_vals, tmp_bestdir_ctrl_vals(1, best_idx) );
        bestdir_drug_vals = horzcat( bestdir_drug_vals, tmp_bestdir_drug_vals(1, best_idx) );  
   
    end
    
    unity_struct.bestdir_ctrl_vals = bestdir_ctrl_vals;
    unity_struct.bestdir_drug_vals = bestdir_drug_vals;
    unity_struct.animal = {stripped_struct.animal};
    unity_struct.signif_visfix = [stripped_struct.signif_visfix];
    unity_struct.signif_attend = [stripped_struct.signif_attend];
    
    
    % Perform Plot on control vs drug values
    [control_vals, drug_vals, rs_pval] = unity_plot( unity_struct, drug, current, unity_type );
end