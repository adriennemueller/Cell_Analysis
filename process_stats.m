% Will run through all processed files and calculate d's, anovas, etc for
% different epochs of the trial and save them to master_file_struct
function mfs = process_stats( mfs )

    if nargin < 1, load master_file_struct; mfs = master_file_struct; end
    
    % If process_stats has not been run before, prepare the master file
    % struct to receive stats output by making appropriate fields.
    mfs = populate_fields( mfs );
    
    % Loop through all processed files
    for i = 1:length(mfs.session)
        
        for j = 1:length(mfs.session(i).processed_files)
            
            % Identify the paradigms present in those files
            data_struct = load_processed_file( mfs.session(i).sub_direc, mfs.session(i).processed_files{j} );
             if length(unique([data_struct.theta])) > 9
                 data_struct = adjust_theta( data_struct );
             end
            paradigm_list = mfs.session(i).paradigms{j};
            currents = mfs.session(i).currents{j};
            
            % Calculate stats for each paradigm (anova, d', etc) and 
            % save substruct of stats for each paradigm into mfs
            if contains( paradigm_list, 'Attention' )
                attend_trial_struct = data_struct( contains({data_struct.paradigm}, 'Attention' ) );
                mfs.session(i).attend_stats{j} = windowed_stats( attend_trial_struct, currents, 'attend' );
                mfs.session(i).attend_visual_stats{j} = windowed_stats( attend_trial_struct, currents, 'visual' );
                mfs.session(i).attend_fixation_stats{j} = windowed_stats( attend_trial_struct, currents, 'fixation' );
            end

            if sum( strcmp( paradigm_list, 'WM' ) )
                wm_trial_struct = data_struct( contains({data_struct.paradigm}, 'WM' ) );
                mfs.session(i).wm_stats{j} = windowed_stats( wm_trial_struct, currents, 'wm' );
                mfs.session(i).wm_visual_stats{j} = windowed_stats( wm_trial_struct, currents, 'visual' );
                mfs.session(i).wm_fixation_stats{j} = windowed_stats( wm_trial_struct, currents, 'fixation' );
            end
            
            if sum( strcmp( paradigm_list, 'Attention_Contrast' ) )
                attContrast_trial_struct = data_struct( contains({data_struct.paradigm}, 'Attention_Contrast' ) );
                mfs.session(i).attContrast_stats{j} = windowed_stats( attContrast_trial_struct, currents, 'attContrast' );
                mfs.session(i).attContrast_visual_stats{j} = windowed_stats( attContrast_trial_struct, currents, 'visual' );
                mfs.session(i).attContrast_fixation_stats{j} = windowed_stats( attContrast_trial_struct, currents, 'fixation' );
            end
    
        end
        
        % Save out master_file_struct
        master_file_struct = mfs;
        save( 'master_file_struct', 'master_file_struct' );
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


function win_stats = windowed_stats( data_struct, currents, window_str )
    
    % Set contrast_flag
    if strcmp(window_str, 'attContrast')
        contrast_flag = 1; else contrast_flag = 0;
    end

    % Go through all different ejected currents in the file. (currents(1)
    % will be the retain current.
    for i = 2:length(currents)
        
        retain_current = currents(1);
        eject_current  = currents(i); 
        
        corr_idx = find( [data_struct.trial_error] == 0 );
        window = get_window( data_struct(corr_idx(1)), window_str );
        
        control_spikemat_attin  = get_directional_spikemat( data_struct, retain_current, window, 'in', contrast_flag );
        control_spikemat_attout = get_directional_spikemat( data_struct, retain_current, window, 'out', contrast_flag  );
        drug_spikemat_attin     = get_directional_spikemat( data_struct, eject_current, window, 'in', contrast_flag );
        drug_spikemat_attout    = get_directional_spikemat( data_struct, eject_current, window, 'out', contrast_flag );
                
        % Get D' for result of this
        control_dmat = gen_dprime_struct_wrapper( control_spikemat_attin, control_spikemat_attout );
        drug_dmat    = gen_dprime_struct_wrapper( drug_spikemat_attin, drug_spikemat_attout );
    
        % Get Anova for this
        anova_mat = gen_anova_struct( control_spikemat_attin, control_spikemat_attout, drug_spikemat_attin, drug_spikemat_attout );
 
        % Return results for each current separately
        win_stats(i-1).current = eject_current;
        win_stats(i-1).control_dmat = control_dmat;
        win_stats(i-1).drug_dmat = drug_dmat;
        win_stats(i-1).anova_mat = anova_mat;
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
       
    end
end


function rslt = gen_anova_struct( ctrl_attin, ctrl_attout, drug_attin, drug_attout )
    data_vec = []; 
    direction = [];
    drug = [];
    attend = [];
    
    % Identify whether this is attend_Contrast paradigm data.
    contrast_flag = isfield( ctrl_attin, 'contrast' );
    if contrast_flag, contrast = []; end    
    
    for i = 1:length(unique([ctrl_attin.direction])) 
        
        direc = ctrl_attin(i).direction;
        if (direc >= 180), continue, end % Only 4 dir b/c copied att in/out.
     
        % Find the Drug Off, Attend In trials and make a vector of the
        % number of spikes for them.
        if ~isempty(ctrl_attin(i).spikes)
        data_vec = [data_vec sum(ctrl_attin(i).spikes)]; 
        ctrl_attin_length = size(ctrl_attin(i).spikes, 2);
        direction = [ direction repmat(direc, 1, ctrl_attin_length )];
        drug = [drug zeros(1, ctrl_attin_length)];
        attend = [attend ones(1, ctrl_attin_length)];
        if contrast_flag, contrast = [contrast ctrl_attin(i).contrast]; end
        end
        
        % Same for Drug Off, Attend Out
        if ~isempty(ctrl_attout(i).spikes)
        data_vec = [data_vec sum(ctrl_attout(i).spikes)]; 
        ctrl_attout_length = size(ctrl_attout(i).spikes, 2);
        direction = [ direction repmat(direc, 1, ctrl_attout_length )];
        drug = [drug zeros(1, ctrl_attout_length)];
        attend = [attend zeros(1, ctrl_attout_length)];
        if contrast_flag, contrast = [contrast ctrl_attout(i).contrast]; end
        end
        
        % Find the Drug On, Attend In trials and make a vector of the
        % number of spikes for them.        
        if ~isempty(drug_attin(i).spikes)
        data_vec = [data_vec sum(drug_attin(i).spikes)]; 
        drug_attin_length = size(drug_attin(i).spikes, 2);
        direction = [ direction repmat(direc, 1, drug_attin_length )];
        drug = [drug ones(1, drug_attin_length)];
        attend = [attend ones(1, drug_attin_length)]; 
        if contrast_flag, contrast = [contrast drug_attin(i).contrast]; end
        end
        
        % Same for Drug On, Attend Out
        if ~isempty(drug_attout(i).spikes)
        data_vec = [data_vec sum(drug_attout(i).spikes)]; 
        drug_attout_length = size(drug_attout(i).spikes, 2);
        direction = [ direction repmat(direc, 1, drug_attout_length )];
        drug = [drug ones(1, drug_attout_length)];
        attend = [attend zeros(1, drug_attout_length)];
        if contrast_flag, contrast = [contrast drug_attout(i).contrast]; end
        end
        
    end
    
    if contrast_flag
        [p,tbl] = anovan( data_vec, {direction, drug, attend, contrast}, 'model','full','varnames',{'direction','drug','attend', 'contrast'}, 'display','off', 'sstype', 1 );
    else
        [p,tbl] = anovan( data_vec, {direction, drug, attend}, 'model','full','varnames',{'direction','drug','attend'}, 'display','off', 'sstype', 1 );
    end
    
    rslt.p = p;
    rslt.tbl = tbl;
end