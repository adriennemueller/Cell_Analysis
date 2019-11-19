% drug_svm takes:
% master_file_struct
% drug_selection (either SCH23390 or SKF81297, currently)
% current_selection: for that drug selection (usually either 20 or 50)
% crossval_n: The number of crossvalidation iterations to run
% jackknife_flag: 1 if data partitioned into 'leave-one-out' partitions for
%   jackknife crossvalidation; 0 if training and test trials drawn randomly
%   crossval_n times
% scramble_flag: 1 if data labels also scrambled for comparison

% It will create one large matrix of trials from cells which were exposed
% to the chosen drug at the chosen current and give each trial a label of
% either 'Control' or 'Drug'. It will then train 'straps'-many SVM
% classifiers on that data, with the matrix being broken into 'straps'-many
% sets of training and test matrices.

% The output struct includes the models, the labels, the scores and the percent of
% the labels that were correct for each model. The struct will also
% include, for each model, the percentage of the training matrix trials
% that were control trials.

% IMPORTANT: Currently will only use Attention Trials
function [CVSVMModel, signtest_pval] = drug_svm( mfs, drug_selection, current_selection, crossval_n, jackknife_flag, scramble_flag )

    if nargin < 6, scramble_flag = 1; end
    if nargin < 5, jackknife_flag = 0; end
    if nargin < 4, crossval_n = 10; end

    % Get giant matrix of trials for all cells with drug selection and
    % current selection
    [svm_mat, drug_labels] = multicell_mat( mfs, drug_selection, current_selection );
    assignin( 'base', 'svm_mat', svm_mat);
    assignin( 'base', 'drug_labels', drug_labels);

    if jackknife_flag
        % Separate data into crossvalidation partitions 
        partition_struct = partition_svm_mat( svm_mat, drug_labels, crossval_n );
    
        % Generating training and test sets from the crossval_n partions
        svm_struct = gen_jkCrossval_svm_struct( partition_struct );
    else
        % Generating training and test sets crossval_n times from the main
        % mat and labels
        svm_struct = gen_randCrossval_svm_struct( svm_mat, drug_labels, crossval_n );
    end

    CVSVMModel = train_test_svm( svm_struct, scramble_flag );
    
    % Run Signtest if data was also scrambled for comparison
    if scramble_flag
        signtest_pval = signtest( [CVSVMModel.perc_corr], [CVSVMModel.scr_perc_corr] );
        disp( strcat( 'Scrambled vs Unscrambled Signtest p-value = ', {' '}, num2str(signtest_pval) ) ); 
    else
        signtest_pval = NaN;
    end
    
end

% Iterates through the masterfilestruct and identifies processed files that
% contain data collected with the chosen drug and current. Returns a large
% matrix of all the trials appended together. Rows are trials. Columns are
% time points. The data currently being appended are the spike values (1 or
% 0) for the 300ms attend window.
function [svm_mat, drug_labels] = multicell_mat( mfs, drug_selection, current_selection )

    svm_mat = [];
    drug_labels = {};
    
    for i = 1:length(mfs.session)
        if strcmp(mfs.session(i).drug, drug_selection)
            for j = 1:length(mfs.session(i).processed_files)
                if ismember(current_selection, mfs.session(i).currents{j})
                    
                    fname = mfs.session(i).processed_files{j};
                    disp(strcat('Concatenating: ', {' '},  fname));
        
                    data_struct = load_processed_file( mfs.session(i).sub_direc, fname );
                    
                    control_trials = filtered_windowed_spikemat( data_struct, -15, 'attend_fixation', [], [], 0 );
                    drug_trials    = filtered_windowed_spikemat( data_struct, current_selection, 'attend_fixation', [], [], 0 );
        
                    % Restrict matrix lengths to be the same widths (times)
                    % by removing columns (time points) from the front.
                    [svm_mat, control_trials, drug_trials] = shave_bins( svm_mat, control_trials, drug_trials );
                    
                    % Append spike matrix to svm_mat and drug_labels to
                    % drug_labels.
                    svm_mat = horzcat(svm_mat, control_trials);
                    svm_mat = horzcat(svm_mat, drug_trials);
                    curr_drug_labels = get_drug_labels( size(control_trials,2), size(drug_trials,2) );
                    drug_labels = horzcat(drug_labels, curr_drug_labels);

                end
            end
        end
    end
end


% Take drug currents and convert them into a list of string labels
function drug_labels = get_drug_labels( control_length, drug_length )
    total_length = control_length+drug_length;
    drug_labels = cell( 1, total_length);

    drug_labels(1:control_length) = {'Control'};
    drug_labels(control_length+1:total_length) = {'Drug'};
end






