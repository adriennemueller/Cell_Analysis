

% svm_parameter:
%  theta
%  attend_direction
%  decision
%  [rewarded] future plan
function svm_plot = svm_drug_effect( mfs, drug, current, epoch, crossval_n )

    if nargin < 5, crossval_n = 10; end

    % Get the subset of cells with relevant drug and current
    cell_struct = filter_cells( mfs, drug, current );
    
    
    % Loop through all relevant cells
    
    % For each epoch relevant to that svm_parameter
   
        % train and test an SVM
        
        % Get multiple estimates of SVM performance for a given cell

    % Once have gone through all cells; have a distribution of performances and
    % compare to chance. (With e.g. an anova? t-test)

end


%
function cell_struct = filter_cells( mfs, drug, epoch, current )

    cell_struct = struct;

    for i = 1:length(mfs.session)
        if strcmp(mfs.session(i).drug, drug)
            for j = 1:length(mfs.session(i).processed_files)
                if ismember(current_selection, mfs.session(i).currents{j})
                    
                    fname = mfs.session(i).processed_files{j};
                    disp(strcat('Adding: ', {' '},  fname));
        
                    data_struct = load_processed_file( mfs.session(i).sub_direc, fname );
                    
                    %%% MAKE THIS RETURN A MATRIX WITH RELEVANT IDENTIFIERS
                    %%% FOR THETA, ATTEND_DIRECTION AND DECISON 
                    control_trials = filtered_windowed_spikemat( data_struct, -15, epoch, [], [], 0 );
                    drug_trials    = filtered_windowed_spikemat( data_struct, current, epoch, [], [], 0 );
        
                    % Restrict matrix lengths to be the same widths (times)
                    % by removing columns (time points) from the front.
                    [svm_mat, control_trials, drug_trials] = shave_bins( svm_mat, control_trials, drug_trials );
                    
                    cell_struct(end+1).control_trials
                    cell_struct(end+1).drug_trials
                    

                end
            end
        end
    end


end


% Break data into training and test sets. This will probably have a struct
% output actually
function [training_data test_data] = partition_data_for_svm( svm_input_data, jackknife_flag )


end

% Train SVM and test it and then have it spit out its performance
function SVM_performance = run_test_SVM( training_data, test_data )

end
