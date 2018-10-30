
% unity_type = 'd_prime' or 'mean_fr'

function gen_unity_plots( mfs, drug, current, paradigm, unity_type )

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
            
            four_pt_ctrl_vals = horzcat( four_pt_ctrl_vals, tmp_four_pt_ctrl_vals(1, [7,8,1,2]) );
            four_pt_drug_vals = horzcat( four_pt_drug_vals, tmp_four_pt_drug_vals(1, [7,8,1,2]) );
        end
        unity_plot( four_pt_ctrl_vals, four_pt_drug_vals, drug, current, 'd_prime' );
    end

    %%% SIGNIFICANT ATTEND ONLY 

    
    %%% FIXATION FIRING RATE
    if strcmp( unity_type, 'mean_fr' )
        four_pt_ctrl_vals = [];
        four_pt_drug_vals = [];

        for i = 1:length(stripped_struct)
        
            tmp_four_pt_ctrl_vals = [stripped_struct(i).attend_fixation_stats.control_summ_stats.avg_fr];
            tmp_four_pt_drug_vals = [stripped_struct(i).attend_fixation_stats.drug_summ_stats.avg_fr];
            
            four_pt_ctrl_vals = horzcat( four_pt_ctrl_vals, tmp_four_pt_ctrl_vals(1, [7,8,1,2]) );
            four_pt_drug_vals = horzcat( four_pt_drug_vals, tmp_four_pt_drug_vals(1, [7,8,1,2]) );
        end
        unity_plot( four_pt_ctrl_vals, four_pt_drug_vals, drug, current, 'mean_fr' );
    end
    
   
    %%% ATTEND WINDOW FIRING RATE
    
    
    
end