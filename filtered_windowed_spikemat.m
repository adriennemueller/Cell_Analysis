% Make spike matrix of (currently only) attend window


%%% NEED DIRECTION AND INOUT FILTERING
%%% NEED TO FIX DRUG SVM CODE SO INCLUDES INOUT
function rslt_mat = filtered_windowed_spikemat( curr_data_mat, current, window, direction, inout )
    
    % Only keep correct attend trials in matrix
    correct_idxs = find( [curr_data_mat.trial_error] == 0 );
    corr_current_idxs = find( [curr_data_mat.drug] == current );
    attend_trial_idxs = find(cell2mat(cellfun( @(codes) ismember(126, codes), {curr_data_mat.event_codes}, 'Uniformoutput', 0 )));
    valid_idxs = correct_idxs( ismember(correct_idxs, attend_trial_idxs) );
    valid_idxs = valid_idxs( ismember( valid_idxs, corr_current_idxs ) );
    
    spike_mat = {curr_data_mat(valid_idxs).spikes};
    millis_mat = {curr_data_mat(valid_idxs).millis};
    rslt_mat = attend_only_window( spike_mat, millis_mat, {curr_data_mat(valid_idxs).event_codes}, {curr_data_mat(valid_idxs).code_times} );
    
end


% Make a smaller spike matrix of just the attend window from 300ms before
% blank on to blank on.
function attend_mat = attend_only_window( spike_mat, millis_mat, event_codes, code_times )

    % Loops over matching cells in event_codes and code_times and applies
    % the find (of eventcode 126 - blank off) to each pair.
    targ_off_times = cellfun(@(codes, times) times(codes == 126), event_codes, code_times);
        
    % Find the index of the milli value that corresponds to that time.
    targ_off_idxs = cellfun(@(millis, times) find(millis == times), millis_mat, num2cell(targ_off_times), 'uniformoutput', 0);
    
    % Grab the chunk of spikes from that index to -300 from that index
    attend_cellarray = cellfun(@(spikes, idxs) spikes( (idxs - 299):idxs ), spike_mat, targ_off_idxs, 'uniformoutput', 0);
    attend_mat = cell2mat(attend_cellarray);
end


