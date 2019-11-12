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
% sets of training and test matrices.8006977089 8002372767 8002413371

% The output struct includes the models, the labels, the scores and the percent of
% the labels that were correct for each model. The struct will also
% include, for each model, the percentage of the training matrix trials
% that were control trials.

% IMPORTANT: Currently will only use Attention Trials
function [CVSVMModel, signtest_pval] = drug_svm( mfs, drug_selection, current_selection, crossval_n, jackknife_flag, scramble_flag )

    if nargin < 6, scramble_flag = 1; end
    if nargin < 5, jackknife_flag = 1; end
    if nargin < 4, crossval_n = 10; end

    % Get giant matrix of trials for all cells with drug selection and
    % current selection
    [svm_mat, drug_labels] = multicell_mat( mfs, drug_selection, current_selection );
    assignin( 'base', 'svm_mat', svm_mat);
    assignin( 'base', 'drug_labels', drug_labels);

    if jackknife_flag
        % Separate data into crossvalidation partitions 
        partition_struct = partition_svm_mat( svm_mat, drug_labels, crossval_n );
    
        % Generating training and test sets from the 10 partions
        svm_struct = gen_jkCrossval_svm_struct( partition_struct );
    else
        % Generating training and test sets crossval_n times from the main
        % mat and labels
        svm_struct = gen_randCrossval_svm_struct( svm_mat, drug_labels, crossval_n );
    end

    CVSVMModel = struct;
    % Loop through svm_struct, generating an SVM for each train/test set
    % and predicting test values. Add predicted test labels to list
    % compared with real labels
    for i = 1:length(svm_struct)
        
        % Train SVM
        disp( strcat( 'Training SVM:', {' '}, num2str(i), '...' ) );
        train_svm_mat = svm_struct(i).train_svm_mat;
        train_drug_labels = svm_struct(i).train_drug_labels;
        test_svm_mat = svm_struct(i).test_svm_mat;
        test_drug_labels = svm_struct(i).test_drug_labels;
        SVMModel = fitcsvm(train_svm_mat, train_drug_labels, 'KernelFunction','rbf', 'Standardize',true,'ClassNames',{'Control','Drug'});
        
        % Test SVM
        [label,score] = predict( SVMModel, test_svm_mat );
        perc_corr = compare_labels( label, test_drug_labels );
        
        CVSVMModel(i).perc_control = sum(strcmp(train_drug_labels, 'Control')) / length(train_drug_labels);
        CVSVMModel(i).Model = SVMModel;
        CVSVMModel(i).label = label; 
        CVSVMModel(i).score = score;
        CVSVMModel(i).perc_corr = perc_corr;
        
        if scramble_flag
            scramble_train_drug_labels = train_drug_labels(randperm(length(train_drug_labels)));
            disp( strcat( 'Training Scrambled SVM:', {' '}, num2str(i), '...' ) );
            Scramble_SVMModel = fitcsvm(train_svm_mat, scramble_train_drug_labels, 'KernelFunction','rbf', 'Standardize',true,'ClassNames',{'Control','Drug'});
            [scr_label, scr_score] = predict( Scramble_SVMModel, test_svm_mat );
            scr_perc_corr = compare_labels( scr_label, test_drug_labels );

            CVSVMModel(i).Scamble_Model = Scramble_SVMModel;
            CVSVMModel(i).scr_label = scr_label; 
            CVSVMModel(i).scr_score = scr_score;
            CVSVMModel(i).scr_perc_corr = scr_perc_corr;
        end
        
    end 
    
    % Run Signtest if data was also scrambled for comparison
    if scramble_flag
        signtest_pval = signtest( [CVSVMModel.perc_corr], [CVSVMModel.scr_perc_corr] );
        disp( strcat( 'Scrambled vs Unscrambled Signtest p-value = ', {' '}, num2str(signtest_pval) ) ); 
    else
        signtest_pval = NaN;
    end
    
end

% Break a matrix of trials and a matrix of data labels into cross_val_n
% number of partitions. Return them in a struct.
function partition_struct = partition_svm_mat( svm_mat, drug_labels, crossval_n )

    partition_struct = struct;

    num_trials = size(svm_mat,1);

    partition_idxs = crossvalind('Kfold', num_trials, crossval_n);
    for i = 1:max(partition_idxs)
        partition_struct(i).svm_mat = svm_mat(partition_idxs == i,:);
        partition_struct(i).orig_labels = drug_labels(partition_idxs == i,:);
    end
end

% Great a structure of length crossval_n, that consists of training and
% test data with associated labels generated by random-draws (crossval_n
% times) from the full trial set
function svm_struct = gen_randCrossval_svm_struct( svm_mat, drug_labels, crossval_n )
    
    svm_struct = struct;
    for i = 1:crossval_n
        % Separate into training sets and test set
        [train_idxs, test_idxs] = separate_trials(size(svm_mat, 1), 0.9); % Train on 90%, leave 10% out
        svm_struct(i).train_svm_mat     = svm_mat( train_idxs, : );
        svm_struct(i).train_drug_labels = drug_labels( train_idxs );
        svm_struct(i).test_svm_mat      = svm_mat( test_idxs, : );
        svm_struct(i).test_drug_labels  = drug_labels(  test_idxs );
    end
end

% Separate a struct containing n trial matrices into n training sets and
% test sets
function svm_struct = gen_jkCrossval_svm_struct( partition_struct )

    svm_struct = struct;

    for left_out_idx = 1:length(partition_struct)
       
        curr_svm_struct = partition_struct;
        curr_svm_struct(left_out_idx) = [];
        
        svm_struct(left_out_idx).train_svm_mat = vertcat(curr_svm_struct.svm_mat);
        svm_struct(left_out_idx).train_drug_labels = vertcat(curr_svm_struct.orig_labels);
        
        svm_struct(left_out_idx).test_svm_mat = partition_struct(left_out_idx).svm_mat;
        svm_struct(left_out_idx).test_drug_labels = partition_struct(left_out_idx).orig_labels;

    end

end

% compare_labels determines how many of the svm labels in a test set were
% correctly assigned
function perc_corr = compare_labels( svm_labels, orig_labels )
    correct = strcmp( svm_labels, orig_labels );
    perc_corr = sum(correct) / length(svm_labels);
end

% separate_trials returns two sets of indexes to use to break a matrix of
% length mat_length into two new matrices. The first set of indices,
% train_idxs, will contain a fraction of the total indices equal to
% train_prop. The second set of indices, test_idxs, will contain the
% remainer.
function [train_idxs, test_idxs] = separate_trials( mat_length, train_prop )
    idx_vec = 1:mat_length;
    perm_idxs = randperm(numel( idx_vec ));

    train_length = round(mat_length * train_prop);

    train_idxs = perm_idxs(1:train_length);
    test_idxs  = perm_idxs(train_length+1:mat_length);
    
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
                    svm_mat = vertcat(svm_mat, control_trials');
                    svm_mat = vertcat(svm_mat, drug_trials');
                    curr_drug_labels = get_drug_labels( size(control_trials,2), size(drug_trials,2) );
                    drug_labels = vertcat(drug_labels, curr_drug_labels');

                end
            end
        end
    end
end

% Restrict matrix lengths to be the same widths (times)
% by removing columns (time points) from the front.
function [svm_mat, control_trials, drug_trials] = shave_bins( svm_mat, control_trials, drug_trials )
    [svm_mat_x svm_mat_y] = size( svm_mat );
    [ctrl_x ctrl_y]       = size( control_trials );
    [drug_x drug_y] = size( drug_trials );


    if svm_mat_x == 0
        shortest_val = min( [ctrl_x, drug_x] );
    else
        shortest_val = min( [svm_mat_y, ctrl_x, drug_x] );
        svm_discrepancy = svm_mat_y - shortest_val;
        svm_mat = svm_mat( :, svm_discrepancy+1:end ); 
    end
    
    ctrl_discrepancy = ctrl_x - shortest_val;
    drug_discrepancy = drug_x - shortest_val;
    
    control_trials = control_trials( ctrl_discrepancy+1:end,:);
    drug_trials    = drug_trials( drug_discrepancy+1:end,:);
    
end

% Take drug currents and convert them into a list of string labels
function drug_labels = get_drug_labels( control_length, drug_length )
    total_length = control_length+drug_length;
    drug_labels = cell( 1, total_length);

    drug_labels(1:control_length) = {'Control'};
    drug_labels(control_length+1:total_length) = {'Drug'};
end






