% Will run through all processed files and calculate d's, anovas, etc for
% different epochs of the trial and save them to master_file_struct
function mfs = process_stats( mfs )

    if nargin < 1, load master_file_struct; mfs = master_file_struct; end
    
    % If process_stats has not been run before, prepare the master file
    % struct to receive stats output by making appropriate fields.
    mfs = populate_fields( mfs );
    
    % Loop through all processed files
    for i = 1:length(mfs.session)
        
        % Debug Code to test a specific session
%         if ~ strcmp(mfs.session(i).sub_direc, '2016.08.26')
%             continue
%         end
        
        for j = 1:length(mfs.session(i).processed_files)
            
            % Identify the paradigms present in those files
            data_struct = load_processed_file( mfs.session(i).sub_direc, mfs.session(i).processed_files{j} );
            paradigm_list = mfs.session(i).paradigms{j};
            currents = mfs.session(i).currents{j};
            
            disp({'Processing', mfs.session(i).processed_files{j}})
            
            % Calculate stats for each paradigm (anova, d', etc) and 
            % save substruct of stats for each paradigm into mfs
            if contains( paradigm_list, 'Attention' )
                attend_trial_struct = data_struct( contains({data_struct.paradigm}, 'Attention' ) );
                mfs.session(i).stats{j}.attend_fixation_stats   = windowed_stats( attend_trial_struct, currents, 'fixation', 'Attention' );
                mfs.session(i).stats{j}.attend_visual_stats     = windowed_stats( attend_trial_struct, currents, 'visual', 'Attention' );
                mfs.session(i).stats{j}.attend_attend_stats     = windowed_stats( attend_trial_struct, currents, 'attend', 'Attention' );
                mfs.session(i).stats{j}.attend_blank_stats      = windowed_stats( attend_trial_struct, currents, 'blank', 'Attention' );
                mfs.session(i).stats{j}.attend_post_blank_stats = windowed_stats( attend_trial_struct, currents, 'post_blank', 'Attention' );
                mfs.session(i).stats{j}.attend_reward_stats     = windowed_stats( attend_trial_struct, currents, 'reward', 'Attention' );
                mfs.session(i).stats{j}.vis_signif = visual_significance( attend_trial_struct, currents );
            end

            if sum( strcmp( paradigm_list, 'WM' ) )
                wm_trial_struct = data_struct( contains({data_struct.paradigm}, 'WM' ) );
                mfs.session(i).stats{j}.wm_fixation_stats = windowed_stats( wm_trial_struct, currents, 'fixation', 'WM' );
                mfs.session(i).stats{j}.wm_visual_stats   = windowed_stats( wm_trial_struct, currents, 'visual', 'WM' );
                mfs.session(i).stats{j}.wm_delay_stats    = windowed_stats( wm_trial_struct, currents, 'wm_delay', 'WM' );
                mfs.session(i).stats{j}.wm_response_stats = windowed_stats( wm_trial_struct, currents, 'wm_response', 'WM' );
                mfs.session(i).stats{j}.wm_reward_stats   = windowed_stats( wm_trial_struct, currents, 'reward', 'WM' );
            end
            
            if sum( strcmp( paradigm_list, 'Attention_Contrast' ) )
                attContrast_trial_struct = data_struct( contains({data_struct.paradigm}, 'Attention_Contrast' ) );
                
                % Analyze irrespective of Contrast
                mfs.session(i).stats{j}.attend_fixation_stats = windowed_stats( attContrast_trial_struct, currents, 'fixation', 'Attention' );
                mfs.session(i).stats{j}.attend_visual_stats = windowed_stats( attContrast_trial_struct, currents, 'visual', 'Attention' );
                mfs.session(i).stats{j}.attend_stats = windowed_stats( attContrast_trial_struct, currents, 'attend', 'Attention' );
                
                % Analyze for each Contrast separately
                mfs.session(i).stats{j}.attendContrast_fixation_stats = windowed_stats( attContrast_trial_struct, currents, 'fixation', 'Attention_Contrast' );
                mfs.session(i).stats{j}.attendContrast_visual_stats = windowed_stats( attContrast_trial_struct, currents, 'visual', 'Attention_Contrast' );
                mfs.session(i).stats{j}.attendContrast_stats = windowed_stats( attContrast_trial_struct, currents, 'attContrast', 'Attention_Contrast' );
            end
    
        end
        
        % Save out master_file_struct
        master_file_struct = mfs;
        save( 'master_file_struct', 'master_file_struct' );
    end
end


function vis_signif = visual_significance( data_struct, currents )

    contrast_flag = 0;
    
    % Go through all different ejected currents in the file. (currents(1)
    % will be the retain current.
    for i = 2:length(currents)
        
        retain_current = currents(1);
        eject_current  = currents(i); 
        
        corr_idx = find( [data_struct.trial_error] == 0 );
        
        % Get fixation period spikemat
        fix_spikemat = get_directional_spikemat( data_struct, retain_current, 'fixation','', contrast_flag );
        
        % Get visual period spikemat
        vis_spikemat = get_directional_spikemat( data_struct, retain_current, 'visual', '', contrast_flag );
        
        fix_mat = []; vis_mat = []; vis_direcs = [];
        for j = 1:8 % Number of directions
            if j > 4, direc = j - 4;
            else, direc = j;
            end
            
            tmp_fixmat = (sum(fix_spikemat( j ).spikes) / size(fix_spikemat(j).spikes, 1)) * 1000;
            tmp_fixdirec = ones( 1, length(tmp_fixmat)) .* direc;
            fix_mat = horzcat( fix_mat, tmp_fixmat );
            
            tmp_vismat = (sum(vis_spikemat( j ).spikes) / size(vis_spikemat(j).spikes, 1)) * 1000;
            tmp_visdirec = ones( 1, length(tmp_vismat)) .* direc;
            vis_direcs = horzcat( vis_direcs, tmp_visdirec );
            vis_mat = horzcat( vis_mat, tmp_vismat );
        end
        
        % Get Anova for this
        data_vec = horzcat( fix_mat, vis_mat );
        fix_vis_idx = horzcat( ones(1,length(fix_mat)), ones(1, length(vis_mat) ) * 2 ); 
        vis_direcs = horzcat( vis_direcs, vis_direcs );
        [p,tbl] = anovan( data_vec, {fix_vis_idx, vis_direcs}, 'model','full','varnames',{'fix_vis', 'direc'}, 'display','off', 'sstype', 1 );
        
        % Return results for each current separately
        vis_signif(i-1).current = eject_current;
        vis_signif(i-1).average_fix_fr = mean( fix_mat );
        vis_signif(i-1).average_vis_fr = mean( vis_mat );
        vis_signif(i-1).ps = p;
        vis_signif(i-1).tbl = tbl;

    end
end

% Clear out old stats
function mfs = populate_fields( mfs )
    fnames = fieldnames( mfs.session );
    if sum( contains( fnames, 'stats' ) )
        stat_fields_idxs = contains( fnames, 'stats' );
        stat_fields = fnames( stat_fields_idxs );
        mfs.session = rmfield( mfs.session, stat_fields );
    end

end


function win_stats = windowed_stats( data_struct, currents, window_str, paradigm )
    
    % Set contrast_flag
    if strcmp(window_str, 'attContrast')
        contrast_flag = 1; else contrast_flag = 0;
    end

    % Go through all different ejected currents in the file. (currents(1)
    % will be the retain current.
    for i = 2:length(currents)
        
        retain_current = currents(1);
        eject_current  = currents(i); 
        
        if strcmp( paradigm, 'WM' )
            
        else
            control_spikemat_attin  = get_directional_spikemat( data_struct, retain_current, window_str, 'in', contrast_flag );
            control_spikemat_attout = get_directional_spikemat( data_struct, retain_current, window_str, 'out', contrast_flag );
            drug_spikemat_attin     = get_directional_spikemat( data_struct, eject_current, window_str, 'in', contrast_flag );
            drug_spikemat_attout    = get_directional_spikemat( data_struct, eject_current, window_str, 'out', contrast_flag );
        end    
        
        % Get D' for result of this
        control_dmat = gen_dprime_struct_wrapper( control_spikemat_attin, control_spikemat_attout );
        drug_dmat    = gen_dprime_struct_wrapper( drug_spikemat_attin, drug_spikemat_attout );
    
        % Get Anova for this
        anova_mat = create_anova_data_mat( data_struct, retain_current, eject_current, window_str, contrast_flag );
        
        % Get Summary Statistics for this window
        control_summ_stats = gen_summ_stats( control_spikemat_attin, control_spikemat_attout ); 
        drug_summ_stats    = gen_summ_stats( drug_spikemat_attin, drug_spikemat_attout );
        
        % Return results for each current separately
        win_stats(i-1).current = eject_current;
        win_stats(i-1).control_dmat = control_dmat;
        win_stats(i-1).drug_dmat = drug_dmat;
        win_stats(i-1).anova_mat = anova_mat;
        win_stats(i-1).control_summ_stats = control_summ_stats;
        win_stats(i-1).drug_summ_stats    = drug_summ_stats;
    end
end


function anova_data_mat = create_anova_data_mat( data_struct, retain_current, eject_current, window_str, contrast_flag )
    
    currents = [retain_current, eject_current];
    spikes = []; theta = []; target_change = []; drug = []; contrast = [];

    for i = 1:length(currents)
    
        % Loop through datastruct, creating a matrix based on the epoch 
        [ctrl_spike_mat, contrasts, events] = filtered_windowed_spikemat( data_struct, currents(i), window_str, [], [], contrast_flag );
 
        spikes        = horzcat( spikes, sum( ctrl_spike_mat, 1 ));
        theta         = horzcat( theta, [events.theta{:}] );
        target_change = horzcat( target_change, [events.target_change{:}] );
        drug          = horzcat( drug, repmat( currents(i), 1, length(events.theta)) );
        contrast      = horzcat( contrast, contrasts );
        %attend        =  
    end
    
    % Just drug as factor - 'Fixation' Window
    if strcmp( window_str, 'fixation' )
        merged_mat = {spikes, drug};
        factors_list = {'drug'};
    % If direction also a factor
    elseif sum( strcmp( window_str, {'visual', 'wm_delay', 'wm_response', 'reward'} ) )
        merged_mat = {spikes, drug, theta};
        factors_list = { 'drug', 'theta' };
    % IF ATTEND A FACTOR - DOESN'T WORK YET
    elseif sum( strcmp( window_str, {'attend', 'blank'} ) )
%         merged_mat = {spikes, drug, theta, attend}; FIX THIS
%         factors_list = { 'drug', 'theta', 'attend' }; FIX THIS
        merged_mat = {spikes, drug, theta};
        factors_list = { 'drug', 'theta' };
    % If whether target flips is a factor
    elseif strcmp( window_str, 'post_blank' )
%         merged_mat = {spikes, drug, theta, attend, target_change}; FIX THIS
%         factors_list = { 'drug', 'theta', 'attend', 'target_change' }; FIX ThiS
        merged_mat = {spikes, drug, theta, target_change}; 
        factors_list = {'drug', 'theta', 'target_change'};
    end        

    % Add Contrast Flag if applicable
    if contrast_flag
        merged_mat = {merge_mat, contrasts};
        factors_list = strcat( factors_list, 'contrast' );
    end
    
    anova_data_mat = get_anova_struct( merged_mat, factors_list );

end

function anova_mat = get_anova_struct( anova_data_mat, factors_list )
    
    data_vec = anova_data_mat{ 1, :};
    factors  = anova_data_mat( 2:end );
    
    [anova_mat.p, anova_mat.tbl] = anovan( data_vec, factors, 'model','full','varnames', factors_list, 'display','off', 'sstype', 1 );
    
end


function rslt = gen_summ_stats( attin_vals, attout_vals) 
    
    for i = 1:length( [attin_vals.direction] )
        
        rslt(i).direction = attin_vals(i).direction;
        
        spikes = horzcat( attin_vals(i).spikes, attout_vals(i).spikes );
        [ms, trials] = size( spikes );
        numspikes = sum( spikes, 1 );
        
        rslt(i).avg_fr = mean(numspikes ./ ms) * 1000;
        rslt(i).std_dev = std(numspikes ./ ms) * 1000;
        rslt(i).std_err = rslt(i).std_dev / sqrt(trials);
        rslt(i).num_trials = trials;
    end

end



% Wrapper for generating dprime statistics. Needed because of contrast
% paradigm - will return single struct when no contrasts present, or deeper
% struct if contrast info present.
function rslt = gen_dprime_struct_wrapper( groupA, groupB )

    if isfield( groupA, 'contrast' )
        contrasts = unique( [groupA.contrast] );
        for i = 1:length(contrasts)
            contrast = contrasts(i);
            rslt(i).contrast = contrast;
            
            groupAmod = cell2struct( num2cell([groupA.direction]), {'direction'} ); 
            Aspikes = cellfun(@(spks, ctrst) spks( :, ctrst == contrast ), {groupA.spikes}, {groupA.contrast}, 'UniformOutput', 0);
            [groupAmod.spikes] = Aspikes{:};
            groupAmod = groupAmod';
          
            groupBmod = cell2struct( num2cell([groupB.direction]), {'direction'} ); 
            Bspikes = cellfun(@(spks, ctrst) spks( :, ctrst == contrast ), {groupB.spikes}, {groupB.contrast}, 'UniformOutput', 0);
            [groupBmod.spikes] = Bspikes{:};
            groupBmod = groupBmod';

            rslt(i).dmat = gen_dprime_struct( groupAmod, groupBmod);
        end
    else
        rslt.dmat = gen_dprime_struct( groupA, groupB);
    end
end


% Get a vector of contrasts for the subset of data that is being processed.
% This function repeats a lot of the work done in
% filtered_windowed_spikemat. Might want to consider adding contrast as a
% filterable variable for filtered_windowed_spikemat instead.
function contrast_vec = get_contrasts( curr_data_mat, current, direction, inout )
    % Filter for correct trials
    correct_idxs = find( [curr_data_mat.trial_error] == 0 );

    % Filter for current
    corr_current_idxs = find( [curr_data_mat.drug] == current );
    valid_idxs = correct_idxs( ismember( correct_idxs, corr_current_idxs ) );
 
    % Filter for (attend) direction
    if strcmp( inout, 'out' )
        direction = reversed(direction);
    end
    if ~ isempty(direction)
         corr_direc_idxs = find( [curr_data_mat.theta] == direction );
         valid_idxs = valid_idxs( ismember( valid_idxs, corr_direc_idxs ) );
    end        
    
    contrast_vec = [curr_data_mat(valid_idxs).contrast];
end

% reversed gives you the opposite direction (1-8) to the one you input. eg
% 0->180, 270->90 etc.
function rslt = reversed(direction)
    rslt = mod((direction + 135), 360) + 45;
    if rslt == 360, rslt = 0; end
end
        
% gen_dprime_struct takes two spike_mats and generates the d' value for
% between respective elements. So, two spikemats with 8 directions each
% would yield a vector of length 8, with each element betin a d' value.
function dmat = gen_dprime_struct( groupA, groupB )   
    dmat = struct;

    for i = 1:length(unique([groupA.direction]))  % For the 8 directions
        
        dmat(i).direction = groupA(i).direction;
        
        groupA_n_spikes = sum(groupA(i).spikes);
        groupB_n_spikes = sum(groupB(i).spikes);
        
        if ~ isempty(groupA_n_spikes) &&  ~ isempty(groupB_n_spikes)
            dmat(i).ranksum_val = ranksum( groupA_n_spikes, groupB_n_spikes );
        else
           dmat(i).ranksum_val = nan;
        end
       
        dmat(i).dprime_val = d_prime( groupA_n_spikes, groupB_n_spikes );
        dmat(i).mod_idx = (mean(groupA_n_spikes) - mean(groupB_n_spikes)) / (mean(groupA_n_spikes) + mean(groupB_n_spikes));
       
    end
end

% 
% function rslt = gen_anova_struct( ctrl_attin, ctrl_attout, drug_attin, drug_attout, contrast_flag )
%     data_vec = []; 
%     direction = [];
%     drug = [];
%     attend = [];
%     
%     % Identify whether this is attend_Contrast paradigm data.
%     %contrast_flag = isfield( ctrl_attin, 'contrast' ); % NO SUCH DATA IN
%     %ATT MAT
%     if contrast_flag, 
%         contrast = [];
%     end    
%     
%     for i = 1:length(unique([ctrl_attin.direction])) 
%         
%         direc = ctrl_attin(i).direction;
%         if (direc >= 180), continue, end % Only 4 dir b/c copied att in/out.
%      
%         % Find the Drug Off, Attend In trials and make a vector of the
%         % number of spikes for them.
%         if ~isempty(ctrl_attin(i).spikes)
%         data_vec = [data_vec sum(ctrl_attin(i).spikes)]; 
%         ctrl_attin_length = size(ctrl_attin(i).spikes, 2);
%         direction = [ direction repmat(direc, 1, ctrl_attin_length )];
%         drug = [drug zeros(1, ctrl_attin_length)];
%         attend = [attend ones(1, ctrl_attin_length)];
%         if contrast_flag, contrast = [contrast ctrl_attin(i).contrasts]; end
%         end
%         
%         % Same for Drug Off, Attend Out
%         if ~isempty(ctrl_attout(i).spikes)
%         data_vec = [data_vec sum(ctrl_attout(i).spikes)]; 
%         ctrl_attout_length = size(ctrl_attout(i).spikes, 2);
%         direction = [ direction repmat(direc, 1, ctrl_attout_length )];
%         drug = [drug zeros(1, ctrl_attout_length)];
%         attend = [attend zeros(1, ctrl_attout_length)];
%         if contrast_flag, contrast = [contrast ctrl_attout(i).contrasts]; end
%         end
%         
%         % Find the Drug On, Attend In trials and make a vector of the
%         % number of spikes for them.        
%         if ~isempty(drug_attin(i).spikes)
%         data_vec = [data_vec sum(drug_attin(i).spikes)]; 
%         drug_attin_length = size(drug_attin(i).spikes, 2);
%         direction = [ direction repmat(direc, 1, drug_attin_length )];
%         drug = [drug ones(1, drug_attin_length)];
%         attend = [attend ones(1, drug_attin_length)]; 
%         if contrast_flag, contrast = [contrast drug_attin(i).contrasts]; end
%         end
%         
%         % Same for Drug On, Attend Out
%         if ~isempty(drug_attout(i).spikes)
%         data_vec = [data_vec sum(drug_attout(i).spikes)]; 
%         drug_attout_length = size(drug_attout(i).spikes, 2);
%         direction = [ direction repmat(direc, 1, drug_attout_length )];
%         drug = [drug ones(1, drug_attout_length)];
%         attend = [attend zeros(1, drug_attout_length)];
%         if contrast_flag, contrast = [contrast drug_attout(i).contrasts]; end
%         end
%         
%     end
%     
%     if contrast_flag
%         [p,tbl] = anovan( data_vec, {direction, drug, attend, cell2mat(contrast)}, 'model','full','varnames',{'direction','drug','attend', 'contrast'}, 'display','off', 'sstype', 1 );
%     else
%         [p,tbl] = anovan( data_vec, {direction, drug, attend}, 'model','full','varnames',{'direction','drug','attend'}, 'display','off', 'sstype', 1 );
%     end
%     
%     rslt.p = p;
%     rslt.tbl = tbl;
% end