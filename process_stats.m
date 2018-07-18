% Will run through all processed files and calculate d's, anovas, etc for
% different epochs of the trial and save them to master_file_struct
function mfs = process_stats( mfs )

    if nargin < 1, load master_file_struct; mfs = master_file_struct; end
    
    % If process_stats has not been run before, prepare the master file
    % struct to receive stats output by making appropriate fields.
    mfs = populate_fields( mfs );
    
    % Loop through all processed files
    for i = 1:length(mfs.session)
        
        if mfs.session(i).processed_stats == 1, continue; end % Skip this session if stats already processed.
        
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
            if contains( 'Attention', paradigm_list )
                attend_trial_struct = data_struct( contains({data_struct.paradigm}, 'Attention' ) );
                mfs.session(i).attend_stats = [ mfs.session(i).attend_stats windowed_stats( attend_trial_struct, currents, 'attend' ) ];
            end

            if contains( 'WM', paradigm_list )
                if ~isfield(mfs.session(i), 'wm_stats'), mfs.session(i).wm_stats = []; end
                wm_trial_struct = data_struct( contains({data_struct.paradigm}, 'WM' ) );
                mfs.session(i).wm_stats = [ mfs.session(i).wm_stats windowed_stats( wm_trial_struct, currents, 'wm' ) ];
            end
            
%             if contains( 'Attenion_Contrast', paradigm_list )
%                 attContrast_trial_struct = data_struct( contains({data_struct.paradigm}, 'Attention_Contrast' ) );
%                 mfs.session(i).processed_files(j).wm_stats = wm_stats( attContrast_trial_struct, currents );
%             end
            
            % Pretrial Also. %%% TODO %%%
            
          
    
        end
        mfs.session(i).processed_stats = 1;

        % Save out master_file_struct
        master_file_struct = mfs;
        save( 'master_file_struct', 'master_file_struct' );
    end
end

function mfs = populate_fields( mfs )
    if ~isfield( mfs.session, 'attend_stats' )
        mfs.session(1).attend_stats = [];
    end
    
    if ~isfield( mfs.session, 'wm_stats' )
        mfs.session(1).wm_stats = [];
    end 

end


function win_stats = windowed_stats( data_struct, currents, window_str )
    
    % Go through all different ejected currents in the file. (currents(1)
    % will be the retain current.
    for i = 2:length(currents)
        
        retain_current = currents(1);
        eject_current  = currents(i); 
        
        corr_idx = find( [data_struct.trial_error] == 0 );
        window = get_window( data_struct(corr_idx(1)), window_str );
        
        control_spikemat_attin  = get_directional_spikemat( data_struct, retain_current, window, 'in' );
        control_spikemat_attout = get_directional_spikemat( data_struct, retain_current, window, 'out'  );
        drug_spikemat_attin     = get_directional_spikemat( data_struct, eject_current, window, 'in' );
        drug_spikemat_attout    = get_directional_spikemat( data_struct, eject_current, window, 'out' );
                
        % Get D' for result of this
        control_dmat = gen_dprime_struct( control_spikemat_attin, control_spikemat_attout );
        drug_dmat    = gen_dprime_struct( drug_spikemat_attin, drug_spikemat_attout );
    
        % Get Anova for this
        anova_mat = gen_anova_struct( control_spikemat_attin, control_spikemat_attout, drug_spikemat_attin, drug_spikemat_attout );
 
        % Return results for each current separately
        win_stats(i-1).current = eject_current;
        win_stats(i-1).control_dmat = control_dmat;
        win_stats(i-1).drug_dmat = drug_dmat;
        win_stats(i-1).anova_mat = anova_mat;
    end
end


function spikemat_struct = get_directional_spikemat( spikemat, current, window, inout  )
    spikemat_struct = struct;
    directions = unique([spikemat.theta]);
    % Loop through directions and place into struct
    for i = 1:length(directions)
        spikemat_struct(i).direction = directions(i);
        spikemat_struct(i).spikes = filtered_windowed_spikemat( spikemat, current, window, directions(i), inout);
    end
end

 
function window = get_window( correct_trial, window_string )

    if strcmp( window_string, 'attend' )
        if find(correct_trial.event_codes == 121) % Need a function for this because I changed the event codes in Jan/Mar 2016 
            window = [121 126];
        else
            window = [133 126];
        end
    elseif strcmp(window_string, 'wm')
        window = [155 161]; 
    elseif strcmp( window_string, 'pre_trial' )
        %%% FILL THIS IN
    end
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
        attend = [attend ones(1, ctrl_attin_length)];end
        
        % Same for Drug Off, Attend Out
        if ~isempty(ctrl_attout(i).spikes)
        data_vec = [data_vec sum(ctrl_attout(i).spikes)]; 
        ctrl_attout_length = size(ctrl_attout(i).spikes, 2);
        direction = [ direction repmat(direc, 1, ctrl_attout_length )];
        drug = [drug zeros(1, ctrl_attout_length)];
        attend = [attend zeros(1, ctrl_attout_length)]; end
        
        
        % Find the Drug On, Attend In trials and make a vector of the
        % number of spikes for them.        
        if ~isempty(drug_attin(i).spikes)
        data_vec = [data_vec sum(drug_attin(i).spikes)]; 
        drug_attin_length = size(drug_attin(i).spikes, 2);
        direction = [ direction repmat(direc, 1, drug_attin_length )];
        drug = [drug ones(1, drug_attin_length)];
        attend = [attend ones(1, drug_attin_length)]; end
        
        % Same for Drug On, Attend Out
        if ~isempty(drug_attout(i).spikes)
        data_vec = [data_vec sum(drug_attout(i).spikes)]; 
        drug_attout_length = size(drug_attout(i).spikes, 2);
        direction = [ direction repmat(direc, 1, drug_attout_length )];
        drug = [drug ones(1, drug_attout_length)];
        attend = [attend zeros(1, drug_attout_length)]; end
        
    end
    
    [p,tbl] = anovan( data_vec, {direction, drug, attend}, 'model','full','varnames',{'direction','drug','attend'}, 'display','off', 'sstype', 1 );
    
    rslt.p = p;
    rslt.tbl = tbl;
end



% adjust_theta takes a trial-struct and changes the thetas to the closest
% '8-direction' theta. Non-ideal, but RFs are large enough that grouping
% them shouldn't change the firing too much. Also, collapsing 3-4
% directions down to 1 direction, if the RFs are small, will only make my
% signal worse.
% Crude method.
function adjusted_trials = adjust_theta( trials )

    for i = 1:length(trials)
        theta = round(trials(i).theta);
        if (theta <= 20) || (340 <= theta);
            trials(i).theta = 0;
        elseif (30 <= theta) && (theta <= 60)
            trials(i).theta = 45;
        elseif (70 <= theta) && (theta <= 110)
            trials(i).theta = 90;
        elseif (120 <= theta) && (theta <= 150)
            trials(i).theta = 135;
        elseif (160 <= theta) && (theta <= 200)
            trials(i).theta = 180;
        elseif (210 <= theta) && (theta <= 240)
            trials(i).theta = 225;
        elseif (250 <= theta) && (theta <= 290)
            trials(i).theta = 270;
        elseif (300 <= theta) && (theta <= 330)
            trials(i).theta = 315;
        end
    end
    
    adjusted_trials = trials;
end