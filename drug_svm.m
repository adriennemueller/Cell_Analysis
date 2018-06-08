% drug_svm takes:
% master_file_struct
% drug_selection (either SCH23390 or SKF81297, currently)
% current: for that drug selection (usually either 20 or 50)
% straps: the number of crossvalidation iterations to run

% It will create one large matrix of trials from cells which were exposed
% to the chosen drug at the chosen current and give each trial a label of
% either 'Control' or 'Drug'. It will then train 'straps'-many SVM
% classifiers on that data, with the matrix being broken into 'straps'-many
% sets of training and test matrices.

% The output struct includes the models, the labels, the scores and the percent of
% the labels that were correct for each model. The struct will also
% include, for each model, the percentage of the training matrix trials
% that were control trials.
function rslt = drug_svm( mfs, drug_selection, current_selection, straps, scramble )
    
    if nargin < 5, scramble = 0; end
    if nargin < 4, straps = 100; end

    % Get giant matrix of trials for all cells with drug selection and
    % current selection
    
    [svm_mat, drug_labels] = multicell_mat( mfs, drug_selection, current_selection );
    
    assignin( 'base', 'svm_mat', svm_mat);
    assignin( 'base', 'drug_labels', drug_labels);

    
    % Alternative Jackknife method of bootstrapping - no replacement. Less
    % robust.
    % SVM_Model = fitcsvm( svm_mat, drug_labels );
    % CVSVMModel = crossval(SVMModel);

    CVSVMModel = struct;
     for i = 1:straps
         
        % Separate into training sets and test set
        [train_idxs, test_idxs] = separate_trials(size(svm_mat, 1), 0.9); % Train on 90%, leave 10% out
        train_svm_mat = svm_mat( train_idxs, : ); train_drug_labels = drug_labels( train_idxs );
        test_svm_mat  = svm_mat( test_idxs, : );  test_drug_labels = drug_labels(  test_idxs );

        % Training SVM
        disp( strcat( 'Training SVM:', {' '}, num2str(i), '...' ) );
        SVMModel = fitcsvm(train_svm_mat, train_drug_labels, 'KernelFunction','rbf', 'Standardize',true,'ClassNames',{'Control','Drug'});
        
        % Test SVM
        [label,score] = predict( SVMModel, test_svm_mat );
        perc_corr = compare_labels( label, test_drug_labels );
        
        CVSVMModel(i).perc_control = sum(strcmp(train_drug_labels, 'Control')) / length(train_drug_labels);
        CVSVMModel(i).Model = SVMModel;
        CVSVMModel(i).label = label; 
        CVSVMModel(i).score = score;
        CVSVMModel(i).perc_corr = perc_corr;
        
        if scramble
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

        
    rslt = CVSVMModel;
    
    signtest_pval = signtest( [CVSVMModel.perc_corr], [CVSVMModel.scr_perc_corr] );
    disp( strcat( 'Signtest p-value = ', {' '}, num2str(signtest_pval) ) ); 
    
    %SVMModel = fitclinear(train_svm_mat, train_drug_labels, 'ClassNames',{'Control', drug_label_string{1} });
    %classLoss = kfoldLoss(CVSVMModel) % The output of this is the generalization rate.
    %ScoreSVMModel = fitPosterior( SVMModel, train_svm_mat, train_drug_labels)
      
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
                    control_trials = filtered_windowed_spikemat( data_struct, -15, 'attend', [] );
                    drug_trials    = filtered_windowed_spikemat( data_struct, current_selection, 'attend', [] );
        
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

% Take drug currents and convert them into a list of string labels
function drug_labels = get_drug_labels( control_length, drug_length )
    total_length = control_length+drug_length;
    drug_labels = cell( 1, total_length);

    drug_labels(1:control_length) = {'Control'};
    drug_labels(control_length+1:total_length) = {'Drug'};
    %drug_labels(control_length+1:total_length) = {strcat( drug_selection, '_', num2str(current_selection))};    
end






