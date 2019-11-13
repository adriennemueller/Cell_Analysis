

% svm_parameter:
%  theta
%  attend_direction
%  decision
%  [rewarded] future plan
function svm_plot = svm_drug_effect( mfs, drug, current, epoch, crossval_n )

    if nargin < 5, crossval_n = 10; end

    
    
    
    
    % For each epoch
     
    % For all cells
    
    window_str = 'attend_fixation';
    
    cell_struct = filter_cells( mfs, drug, current, window_str );

    % Loop through all relevant svm_parameters

    
    
   
        % train and test an SVM
        
        % Get multiple estimates of SVM performance for a given cell

    % Once have gone through all cells; have a distribution of performances and
    % compare to chance. (With e.g. an anova? t-test)

end


%
function cell_struct = filter_cells( mfs, drug, current, window_str )

    contrast_flag = 0; % Could be changed later.

    cell_struct = struct;

    for i = 1:length(mfs.session)
        if strcmp(mfs.session(i).drug, drug)
            for j = 1:length(mfs.session(i).processed_files)
                if ismember(current, mfs.session(i).currents{j})
                    
                    fname = mfs.session(i).processed_files{j};
                    disp(strcat('Adding: ', {' '},  fname));
        
                    data_struct = load_processed_file( mfs.session(i).sub_direc, fname );
                    
                    factored_mat = factored_data_mat( data_struct, -15, current, window_str, contrast_flag );
                    cell_struct(end+1).spikes = factored_mat.spikes;
                    cell_struct(end+1).factors = factored_mat.factors;
        
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
