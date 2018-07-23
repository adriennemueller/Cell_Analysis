% sanitize_structs makes a set of structs - one for each unit/bhv
% combination in the session. These structs are a parred-down version of
% the monkeylogic behavior information, and the aligned unit data. The
% output of this function is a cell array of structs; each struct
% containing the sanitized trial-data for a unit/bhv pairing. The times of
% the spikes and the events have been aligned on trial-onset (now 0).
% Because we take more unit data than we have bhv_data for (before the
% trial starts and after it ends) - there are two separate 'millis', one
% for behavioral data, and one for the unit data - that contain the
% relevant aligment times. They are coregistered.
function [session_struct valid_trial_struct_idx] = sanitize_structs( raw_struct, alignments )

    session_struct = {};
    valid_trial_struct_idx = [];

    % For each unit/bhv-file pair - Make a new trial-based struct, throwing
    % out trials without neural data.
    for i = 1:length(alignments)
        bhv  = raw_struct.bhvs(alignments(i).bhv_file);
        unit = raw_struct.units(alignments(i).unit_file);
        
        % Remove trials in from ML bhv data that have no associated PL
        % data. (Generated when ML is left running after PL has been turned
        % off).
        bhv = truncate_trials( bhv, length(alignments(i).PL_bhv_times(:,1)), alignments(i).PL_starting_index);
        
        trial_struct = struct;
        % Iterate through trials
        for j = 1:length(bhv.CodeTimes)
           
            %%% ADD AN OFFSET IN CASE PL STARTS AFTER ML HAS STARTED?
             %%%%%% ctd.animal = bhv.SubjectName; FOR DEGREE CONVERSION
            
            pl_beg = alignments(i).PL_bhv_times(j,1);
            pl_end = alignments(i).PL_bhv_times(j,2);
            unit_beg = alignments(i).unit_beg_time;
            unit_end = alignments(i).unit_end_time;
          
            % if unit activity is recorded for this trial
            if (unit_beg <= pl_beg ) && (unit_end >= pl_end)
                
                % Make a new struct with bhv_info for this trial
                bhv_code_times = bhv.CodeTimes{j};
                offset = bhv_code_times(1);
                % Subtract ~110ms offset caused by ML collecting ~110ms of analogue data before trial start.
                bhv_code_times = bhv_code_times - offset; 
                
                % Take only the relevant info from the original
                % monkey-logic bhv struct
                clean_trial_data = get_clean_trial_data( bhv, j, offset, bhv_code_times );
                
                % Get spike times aligned to Code Times, -500ms before and +500ms after
                % Aligned to TRIAL START
                padding = 500; % How much padding we want before and after a trial in ms
                spike_times = get_aligned_spike_times( unit.SpikeTimes, pl_beg, pl_end, padding );
                spike_time_series = make_ts( spike_times, pl_end - pl_beg + 1, padding );
                clean_trial_data.spike_millis = spike_time_series(2,:)';    
                clean_trial_data.spikes       = spike_time_series(1,:)'; 
                clean_trial_data.sden150      = spike_density( clean_trial_data.spikes, 50 ); %150ms std on Gaussian
              
                % Add them as a field for this struct and append this struct to trial_struct
                if isempty(fieldnames(trial_struct))
                    trial_struct = clean_trial_data;
                else
                    trial_struct = [trial_struct clean_trial_data];
                end               
                
            end
            

        end

        % Remove any trials with fewer than 5 spikes because too sparse.
        %%% CONSIDER ONLY DOING THIS IF CONSECUTIVE FOR STRETCH OF MORE
        %%% THAN 5 trials or so.
        trial_struct = remove_spikeless_trials( trial_struct, 5 );
        
        % Add the drug field to the trial-struct and remove inappropriate
        % drug trials. 
        trial_struct = adjust_drug( trial_struct );
        
        % Add the paradigm field to the trial-struct.
        trial_struct = identify_paradigm( trial_struct, bhv.ConditionsFile );

        % Make sure there are enough trials of each unique drug condition
        % in this trial_struct. If not empty the trial_struct
        [trial_struct, keep] = verify_drug_trialcounts( trial_struct, 20 );
        if keep, valid_trial_struct_idx = [valid_trial_struct_idx i]; end
        
        % Add the full trial-struct for this session to the session-struct.
        % (Each element a unit/bhv pairing).
        session_struct{i} = trial_struct;
    end
   
end


% Go through trial_struct and identify, for each trial, whether it is an
% attention trial, a working memory trial, a probe trial or an
% attention-contrast trial. Unique event codes are:
% Probe Trials:         153
% Attention:            126
% Attention Contrast:   126, need bhv_filename
% Working Memory:       155 or 161
%%% NEED TO DEAL WITH ALL THE UNKNOWNs, which are probably error trials...
%%% or maybe probe trials.
% Are probe trials tossed out in early sessions.
% Do probe rials and wm trials share the same code(wm used)?
function trials = identify_paradigm( trials, full_bhvfile )


    for i = 1:length(trials)
        if max([trials(i).event_codes] == 000) %%% NEED ACTUAL PROBE TRIAL IDENTIFIER
            trials(i).paradigm = 'Probe';
        elseif max([trials(i).event_codes] == 126) && contains(full_bhvfile, 'Contrast')
            trials(i).paradigm = 'Attention_Contrast';
        elseif max([trials(i).event_codes] == 126)
            trials(i).paradigm = 'Attention';
        elseif max([trials(i).event_codes] == 155) % NEED SEPARATE ONE FOR MGS ONLY FILES?
            trials(i).paradigm = 'WM';
        else
            trials(i).paradigm = 'Unknown';
        end
        
    end

end


% Go through trial_struct and remove any trials with fewer than min_spikes spikes
function trial_struct = remove_spikeless_trials( trial_struct, min_spikes )

    spike_sums = arrayfun(@(x) sum(x.spikes),trial_struct);
    spikeless_trials = find( spike_sums < min_spikes );
    trial_struct(spikeless_trials) = [];

end



% Go through all unique currents in a trial_struct; if there are fewer than
% 20 trials for any unique current, remove those trials from the
% trial_struct. If there are not two unique currents - one of which is -15,
% do not keep the trial_struct.
function [trial_struct, keep] = verify_drug_trialcounts( trial_struct, min_trial_num )
  
    % Find unique currents
    unique_currents = unique([trial_struct.drug]);
    keep = 1;
    toss_curr = [];
    
    % Go through each current and identify whether there are enough trials
    % or not. If there are not, remove those trials from trial_struct.
    for i = 1:length(unique_currents)        
        curr = unique_currents(i);
        curr_idxs = find( [trial_struct.drug] == curr );
        if length( curr_idxs ) < min_trial_num
            trial_struct( curr_idxs ) = [];
            toss_curr = [toss_curr i];
        end
        
        % For this current - go through all of the paradigms and ensure
        % that there are enough correct trials of a given paraidgm. If not;
        % eliminate those trials from the trial struct.
        unique_paradigms = unique( {trial_struct.paradigm} );
        for j = 1:length( unique_paradigms )       
            corr_paridigm_idxs = find( ([trial_struct.drug] == curr) &  ...
                                       (strcmp( {trial_struct.paradigm}, unique_paradigms(j))) & ...
                                       ([trial_struct.trial_error] == 0) );
            if length( corr_paridigm_idxs ) < min_trial_num
                paradigm_idxs = find( ([trial_struct.drug] == curr) &  ...
                                       (strcmp( {trial_struct.paradigm}, unique_paradigms(j))) );
                trial_struct( paradigm_idxs ) = []; 
            end
        end
        
    end
    unique_currents(toss_curr) = [];
    
    % If there are fewer than two currents and there is no retain current
    % info, report not to keep this trial_struct.
    now_unique_currents = unique([trial_struct.drug]);
    if (length(now_unique_currents) < 2) || ~ max(now_unique_currents == -15)
        keep = 0;
    end
    
end


% truncate_trials takes an ML bhv struct and the number of trials in length
% that it *should* be and removes trials from the end until it is the right
% length. This is necessary when the ML file outruns the PL file.
function bhv = truncate_trials( bhv, n_trials, starting_index )
  
    if starting_index <= 0
        starting_index_pos = - starting_index;
        selection = starting_index_pos +2 : (n_trials + starting_index_pos - 1);
    else
        selection = 1:n_trials;
    end
    
    bhv.TrialNumber = bhv.TrialNumber(selection,:);
    bhv.AbsoluteTrialStartTime = bhv.AbsoluteTrialStartTime(selection,:);
    bhv.BlockNumber = bhv.BlockNumber(selection,:);
    bhv.BlockIndex = bhv.BlockIndex(selection,:);
    bhv.ConditionNumber = bhv.ConditionNumber(selection,:);
    bhv.TrialError = bhv.TrialError(selection,:);
    bhv.CycleRate = bhv.CycleRate(selection,:);
	bhv.MinCycleRate = bhv.MinCycleRate(selection,:);
    bhv.NumCodes = bhv.NumCodes(selection,:);
    bhv.CodeNumbers = bhv.CodeNumbers(:,selection);
    bhv.CodeTimes = bhv.CodeTimes(:,selection);
    bhv.AnalogData = bhv.AnalogData(:,selection);
    bhv.ReactionTime = bhv.ReactionTime(:,selection);
    bhv.ObjectStatusRecord = bhv.ObjectStatusRecord(:,selection);
    bhv.RewardRecord = bhv.RewardRecord(:,selection);
    bhv.UserVars = bhv.UserVars(:,selection);
end

% Returns a struct with a subset of santized values from the original ML
% behavior struct. Again, all times have been aligned to Trial Start.
function ctd = get_clean_trial_data( bhv, idx, offset, bhv_code_times )
    
    ctd = struct;
    ctd.posx = bhv.AnalogData{1,idx}.EyeSignal(:,1);  %%% CONVERT TO DEG Vis Ang
    ctd.posy = bhv.AnalogData{1,idx}.EyeSignal(:,1);  %%% CONVERT TO DEG Vis Ang
    ctd.lever = bhv.AnalogData{1, idx}.General.Gen1;
    ctd.drugI = (bhv.AnalogData{1, idx}.General.Gen3) .* 100; % Converted to nA
    ctd.pupil = bhv.AnalogData{1, idx}.General.Gen2; %%% Units?
    ctd.millis = 1:length(ctd.posx); ctd.millis = ctd.millis' - offset;
    ctd.code_times = bhv_code_times;
    ctd.block = bhv.BlockNumber(idx);
    ctd.trial_error = bhv.TrialError(idx);
    ctd.event_codes = bhv.CodeNumbers{idx};
    
    reward_on = bhv.RewardRecord(idx).RewardOnTime;
    if reward_on
        ctd.reward_on = bhv.RewardRecord(idx).RewardOnTime - offset;
        ctd.reward_off = bhv.RewardRecord(idx).RewardOffTime - offset;
    else
        ctd.reward_on = [];
        ctd.reward_off = [];
    end
    
    ctd.comb = bhv.UserVars(idx).comb;
    ctd.nocue = bhv.UserVars(idx).nocue;
    ctd.theta = bhv.UserVars(idx).theta * (180/pi); % Converted to Degrees
    ctd.theta = adjust_theta(ctd.theta);
    ctd.radius = bhv.UserVars(idx).radius; % Still in MLUs. Must change.
    ctd.contrast = map_contrast(bhv.ConditionNumber(idx));
    
end

function contrast = map_contrast( condition_number)
    if condition_number <= 700
        contrast = 10;
    elseif condition_number <= 1400
        contrast = 15;
    elseif condition_number <= 2100
        contrast = 20;
    elseif condition_number <= 2900
        contrast = 25;
    end
end


% Change thetas that are in negative degrees into positive degrees.
function theta = adjust_theta( val )
    if val < 0
        theta = 360 + val;
    else
        theta = val;
    end
end


% adjust_drug does two things. 1) it removes trials in which the drug
% current changes by more than 5nA from the trial-struct. 2) It adds a new
% field which states what the actual (rounded, averaged) current applied
% was for the whole trial.
function trials = adjust_drug( trials )
    
    skip_idxs = [];

    for i = 1:length(trials)
        drug = trials(i).drugI;
        
        % Trials in which the current was being changed
        if range(drug) > 5
            skip_idxs = [skip_idxs i];
        end
        
        current = round(mean(drug));
        current = adjust_current(current);
        trials(i).drug = current;
    end
    
    trials(skip_idxs) = [];
end

% Adjust the current values so they are consistent within ~2nA error.
function current = adjust_current( current )
    current = round( current / 5 ) * 5;
end

% Returns a list of spike-times re-aligned to trial-onset.
function aligned_spike_times = get_aligned_spike_times( spike_times, trial_beg, trial_end, padding )

    STs_idxs = find( (spike_times >= (trial_beg - padding)) & (spike_times <= (trial_end + padding)));
    
    aligned_spike_times = spike_times(STs_idxs) - trial_beg;

end

% Makes a time series from a list of spike times for a trial with a given
% duration and padding. Returns not only the list of spike times [0 0 0 0 1 0], but also
% the associated trial-onset-aligned time values [ -2 -1 0 1 2 3 ];
function rslt = make_ts( times, duration, padding )

    length = duration + (2*padding);
    time_series = zeros(1,length);
    
    % idx when spikes occur
    ones_idxs = times + padding +1;

    if sum(times) % If there are actually spikes for this trial
        time_series(ones_idxs) = 1;
    end
    millis = (1:length) - padding;

    rslt = [time_series; millis];
end



