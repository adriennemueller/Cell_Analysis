% preprocess either takes a temporary file struct (or, if none is
% specified, loads the master file struct) and finds sessions of files  which have yet
% to be preprocessed. It then loads those files, finds the alignments - and
% which bhv files go with which units - generates set of new, clean structs
% for each bhv/unit pair, saves them out and (if no tmp struct is specified)
% updates the master file struct.
function preprocess(tmp_struct)

    % Load Master File Struct
    if nargin < 1
        load( 'master_file_struct', 'master_file_struct' );
    else
        master_file_struct = tmp_struct;
    end

    % Preprocess all sessions which have not yet been preprocessed.

    % Find indexes of sessions in master struct which haven't been preprocessed.
    proc_idxs = find( [master_file_struct.session.preprocessed] == 0 );
    disp(['Unprocessed Sessions = ' num2str( length(proc_idxs) )] );
    
    for i = 1:length(proc_idxs)
        
        % Print out which folder is currently being processed
        disp( strcat('Processing: ', master_file_struct.session( proc_idxs(i) ).sub_direc) );
        
        % Load data files into workspace
        raw_struct = load_data_files( master_file_struct.main_direc, master_file_struct.session( proc_idxs(i) ) );
        assignin( 'base', 'raw_struct', raw_struct);
        
        % Identify which unit files belong with which bhv files, and their alignment
        alignments = find_alignments( raw_struct );
        assignin( 'base', 'alignments', alignments);
      
        % Make 'clean' structs for future analyses
        session_struct = sanitize_structs(raw_struct, alignments);
        assignin('base', 'session_struct', session_struct);

        % Save Out the Clean Session Sub-Structs, Add them to the Master File Struct
        processed_files = export_sessions( master_file_struct.main_direc, master_file_struct.session(proc_idxs(i)), alignments, session_struct);
        master_file_struct.session(proc_idxs(i)).processed_files = processed_files;
        
        % Note that this session has now been preprocessed.
        master_file_struct.session( proc_idxs(i) ).preprocessed = 1;

        % Will only save out, and update, master_file_struct if told to
        % work on master_file_struct by lack of input argument
        if nargin < 1
        save('master_file_struct', 'master_file_struct');
        end
        
    end
    


end


% export_sessions loops through and exports clean structs for the found
% bhv/unit pairs for the session and updates the main_file_struct with the
% filenames for those output files.
function processed_files = export_sessions( main_direc, file_struct, alignments, session_struct )
    processed_files = {};

    for i = 1:length(alignments)
        unit_file_idx = alignments(i).unit_file;
        unit_file = file_struct.unit_files{unit_file_idx};
        
        save_fname = ['PROC_', unit_file];
        save_direc = [main_direc, file_struct.sub_direc];
        
        % Save out individual cleaned data struct for the bhv/unit
        % combinations in the session
        data_struct = session_struct{i};
        save( [save_direc, save_fname], 'data_struct');
        
        % Append the filename info to the file_struct that was passed in
        processed_files{i} = save_fname;
    end

end

% matches returns a matrix of alignment for behavioral and unit data. Each
% row represents a match, eg: [ 1 3 407 ] indicates that the first
% behavioral file in raw_struct aligns with the third unit in raw_struct at
% __________ index 407.
function alignments = find_alignments( raw_struct )

    alignments = struct;

    n_units = length( raw_struct.units );
    n_bhvs  = length( raw_struct.bhvs );
    
    for i = 1:n_bhvs

        % The times of some events recorded by Monkeylogic, relative to the trial beginning.
        ML_intra_trial_event_times = raw_struct.bhvs(i).CodeTimes;
        
        % The times that Plexon records events, relative to the beginning of
        % recording.
        PL_event_times = raw_struct.PL_events.TimeStamps;

        % Event numbers from Plexon (matching the timestamps in plexon_event_times)
        % are in Strobed. Find triples of the start and end event codes (which
        % appear to be recorded as 242 and 233, respectively).
        % PL_trial_idx has 2 columns, representing the start and end indexes of
        % trials from Plexon.
        PL_trial_idx = candidate_trial_indexes( raw_struct.PL_events.Strobed, 242, 233 );

        % Timestamps of the start and end events from Plexon
        PL_trial_times = PL_event_times(PL_trial_idx);
        % Duration of each trial from Plexon, in milliseconds
        PL_trial_lengths = PL_trial_times(:,2) - PL_trial_times(:,1);


        % For every trial from BHV, attempt to find the trial from Plexon
        % that most closely matches it in length. Store the index of the Plexon
        % trial into idx(i), matching the BHV trial i.
        idx = zeros(length(ML_intra_trial_event_times), 1);
        for j = 1:length(ML_intra_trial_event_times)
            times = ML_intra_trial_event_times{j};
            trial_length = (times(end) - times(1));
            differ = abs(PL_trial_lengths - trial_length);
            [~, k] = min(differ);
            idx(j) = k(1);
        end

        % Assuming that the BHV trials represent a contiguous slice of trials from
        % the Plexon data, find the index of the start of that run that would match
        % each Plexon trial to the one of closest length in BHV.
        % If the first BHV trial has been matched with the 180th Plexon trial, that
        % would indicate a start index of 180. If the second trial was matched with
        % Plexon trial 201, that would indicate a starting index of 200 (Plexon trial
        % 200 would line up with BHV trial 1). 
        candidate_starting_indexes = idx - (0:length(idx)-1)';
        % Find all suggested start indexes and their frequency of suggestion.
        bins = unique(candidate_starting_indexes);
        freq = histc(candidate_starting_indexes, bins);
        % Take the one that was suggested the most times by trial matches. In
        % practice, the highest fequency match is so by a significant margin.
        [total, j] = max(freq);
        probable_starting_index = bins(j);
        
        
        % The timestamps of the start and end trial events from Plexon, that we
        % think correspond to the trials from BHV. This is an array of two columns,
        % for the start and end timestamps, and the row indexes correspond to trial
        % indexes from BHV.
        ML_n_trials = length( ML_intra_trial_event_times );
        if probable_starting_index+ML_n_trials-1 > length(PL_trial_times)
            probable_end_index = length(PL_trial_times);
        else
            probable_end_index = probable_starting_index+ML_n_trials-1;
        end            
        PL_bhv_times = PL_trial_times(probable_starting_index:probable_end_index, :);
        PL_bhv_beg_time = PL_bhv_times(1,1);
        PL_bhv_end_time = PL_bhv_times(end,2);
        
                        
        
        if total > 10; % Should have at least 10 trials in any session being analyzed. 
            
            % Find associated unit files
            for j = 1:n_units
                unit = raw_struct.units(j);
                unit_beg_time = unit.SpikeTimes(1);
                unit_end_time = unit.SpikeTimes(end);
                
                % If unit starts between bhv_file start and end_times:
                % OR if bhv_file starts between unit start and end times:
                if ((unit_beg_time <= PL_bhv_beg_time) && (PL_bhv_beg_time <= unit_end_time)) ...
                    || ((PL_bhv_beg_time <= unit_beg_time) && (unit_beg_time <= PL_bhv_end_time))
                    
                    tmp_struct.bhv_file  = i;
                    tmp_struct.unit_file = j;
                    tmp_struct.unit_beg_time = unit_beg_time;
                    tmp_struct.unit_end_time = unit_end_time;
                    tmp_struct.PL_bhv_beg_time = PL_bhv_beg_time;
                    tmp_struct.PL_bhv_end_time = PL_bhv_end_time;
                    tmp_struct.PL_bhv_times = PL_bhv_times;
                    
                    % Append struct for this matched unit file and bhv file
                    if isempty(fieldnames(alignments))
                        alignments = tmp_struct;
                    else
                        alignments = [alignments tmp_struct];
                    end
                end
                
            end
            
        end
 
    end

end

% Finds triplets of start_event codes and end_event codes and returns the
% indices of the first of those triplets in an n by 2 matrix. 1 column for
% start event indices, 1 column for end event indices.
function rslt = candidate_trial_indexes( all_events, trial_start_event, trial_end_event )

    ntrials = floor(sum(all_events == trial_start_event) / 3);
    rslt = zeros(ntrials, 2);

    trial_idx = 1;
    event_idx = 1;
    while 1
        trial = find_trial( all_events, event_idx, trial_start_event, trial_end_event );
        if isempty(trial)
            break;
        end
        rslt(trial_idx,:) = trial;
        trial_idx = trial_idx + 1;
        event_idx = trial(2) + 1;
    end

    rslt = rslt(1:trial_idx-1,:);
end

function rslt = find_trial( all_events, start_index, trial_start_event, trial_end_event )

    rslt = [];

    if start_index + 2 > length(all_events)
        return;
    end

    starters = repmat(trial_start_event, 3, 1);
    enders = repmat(trial_end_event, 3, 1);
    start_idx = 0;
    for i = start_index:length(all_events)-2
        ev = all_events(i:i+2);
        if all(ev == starters)
            start_idx = i;
        elseif all(ev == enders)
            if start_idx > 0
                rslt = [start_idx, i+2];
                break;
            end
        end
    end

end




% load_data_files takes a session struct (sub struct of master struct) and
% loads the individual components into a raw structure. The output of this
% function has not be cleaned up, aligned, unified to the same time units,
% etc.
function raw_struct = load_data_files( main_direc, session_struct )

    % Load the plexon timestamped events
    event_file = strcat( main_direc, session_struct.sub_direc, session_struct.event_file );
    PL_events_tmp = load( event_file );
    PL_events.TimeStamps = PL_events_tmp.ans.Ts;
    PL_events.Strobed    = PL_events_tmp.ans.Strobed;
    PL_events.TimeStamps = round(PL_events.TimeStamps * 1000);

    
    % Load all of the monkeylogic bhv files for this session.
    for i = 1:length(session_struct.bhv_files)
        bhv_file = strcat( main_direc, session_struct.sub_direc, session_struct.bhv_files(i) );
        bhvs_struct(i) = bhv_read( bhv_file{1} ); % Accessing {1} to convert from cell array to string
    end
    
    % Load all of the plexon unit files for this session.
    for i = 1:length(session_struct.unit_files)
        unit_file = strcat( main_direc, session_struct.sub_direc, session_struct.unit_files(i) );
        unit_struct_tmp = load( unit_file{1} ); % Accessing {1} to convert from cell array to string
        units_struct(i).SpikeTimes = unit_struct_tmp.adc041(:,1);
        units_struct(i).SpikeTimes = round(units_struct(i).SpikeTimes * 1000);
        units_struct(i).SpikeVolts = unit_struct_tmp.adc041(:,[2:end]);
        
    end
    
    % Output everything into raw_struct
    raw_struct.PL_events = PL_events;
    raw_struct.bhvs = bhvs_struct;
    raw_struct.units = units_struct;
   
end
%Maybe rename raw_struct to session_struct

