
function rslt = drug_svm( mfs )

    svm_mat = [];
    drug_labels = {};

    % Go through all units and make matrices / vectors of drug off and drug on
    % - for the different drugs individually
    % Use only attend window from targ blank start - 300ms (or - 500?)
    
    num_sessions = length(mfs.session);
    
    for i = 1:num_sessions
       
        num_units = length(mfs.session(i).data_mat);
        
        for j = 1:num_units
            % Make spike matrix of attend window
            curr_data_mat = mfs.session(i).data_mat(j).data_mat;
            
            % Only keep correct trials in matrix
            valid_idxs = find( [curr_data_mat.trial_error] == 0 );
            spike_mat = {curr_data_mat(valid_idxs).spikes};
            millis_mat = {curr_data_mat(valid_idxs).millis};
            attend_spike_mat = attend_only_window( spike_mat, millis_mat, {curr_data_mat(valid_idxs).event_codes}, {curr_data_mat(valid_idxs).code_times} );
            
            % Get list of drug vals for correct trials
            drug_vals = [curr_data_mat(valid_idxs).drug]; 
            curr_drug_labels = relabel_drug_vals( drug_vals, mfs.session(i).currents(j), mfs.session(i).drug );
            
            % Append spike matrix to svm_mat and drug_labels to
            % drug_labels.
            svm_mat = vertcat(svm_mat, attend_spike_mat');
            drug_labels = vertcat(drug_labels, curr_drug_labels');
        end
        
    end
    
    

    % Separate into training sets and test set
    
    % Train and crossvalidate
    
    % Test
    
    
    
    
end


% Take drug currents and convert them into a list of string labels
function drug_labels = relabel_drug_vals( drug_vals, currents, drug )
    currents = currents{1};
    drug_labels = cell(1,length(drug_vals));

    control_idxs = find(drug_vals == -15 );
    drug_labels(control_idxs) = {'Control'};
    
    for i = 2:length(currents)
        drug_curr = currents(i);
        current_idxs = find(drug_vals == drug_curr);
        drug_labels(current_idxs) = {strcat( drug, '_', num2str(drug_curr))};
    end
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







