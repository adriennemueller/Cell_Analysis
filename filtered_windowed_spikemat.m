% Make spike matrix of (currently only) attend window

%%% NEED TO FIX DRUG SVM CODE SO INCLUDES INOUT
%%% Filter for Contrasts?     contrasts = [10 15 20 25];


% Beautiful way to identify attend trials - no longer necessary:
% attend_trial_idxs = find(cell2mat(cellfun( @(codes) ismember(126, codes), {curr_data_mat.event_codes}, 'Uniformoutput', 0 )));

function rslt_mat = filtered_windowed_spikemat( curr_data_mat, current, window, direction, inout )
    
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
    
    % Filter for window
    spike_mat = {curr_data_mat(valid_idxs).spikes};
    millis_mat = {curr_data_mat(valid_idxs).millis};
    rslt_mat = extract_window( window, spike_mat, millis_mat, {curr_data_mat(valid_idxs).event_codes}, {curr_data_mat(valid_idxs).code_times} );
    
end


% Make a smaller spike matrix of just the specified window
function windowed_mat = extract_window( window, spike_mat, millis_mat, event_codes, code_times )
    win_beg = window(1);
    win_end = window(2);
    
    % Loops over matching cells in event_codes and code_times and applies
    % the find of the win_beg and win_end event_codes to each pair.
    win_beg_times = cellfun(@(codes, times) times(codes == win_beg), event_codes, code_times);
    win_end_times = cellfun(@(codes, times) times(codes == win_end), event_codes, code_times);
        
    % Find the indices of the milli value that corresponds to those times.
    win_beg_idxs = cellfun(@(millis, times) find(millis == times), millis_mat, num2cell(win_beg_times), 'uniformoutput', 0);
    win_end_idxs = cellfun(@(millis, times) find(millis == times), millis_mat, num2cell(win_end_times), 'uniformoutput', 0);
    
    % Grab the chunk of spikes from beginning index to end index
    windowed_cellarray = cellfun(@(spikes, beg_idxs, end_idxs) spikes( beg_idxs:end_idxs ), spike_mat, win_beg_idxs, win_end_idxs, 'uniformoutput', 0);
    
    % Shave off any excess time_bins so can make a matrix. Shave off from
    % the front. This feels a little dodgy, but I think it's okay as long
    % as acknowledge in methods.
    windowed_cellarray = truncate( windowed_cellarray );
    
    % Convert to matrix
    windowed_mat = cell2mat(windowed_cellarray);
end

function out_array = truncate( in_array )
    lengths = cellfun(@(in_arr) length(in_arr), in_array );
    min_length = min(lengths);
    
    out_array = cellfun(@(in_arr) in_arr( (end-min_length+1):end ), in_array, 'uniformoutput', 0);
end


% reversed gives you the opposite direction (1-8) to the one you input. eg
% 0->180, 270->90 etc.
function rslt = reversed(direction)
    rslt = mod((direction + 135), 360) + 45;
    if rslt == 360, rslt = 0; end
end


