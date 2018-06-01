% drug_svm_jackknife takes:
% master_file_struct
% drug_selection (either SCH23390 or SKF81297, currently)
% current: for that drug selection (usually either 20 or 50)
% cross_val_n: the number of data partitions for cross-validation; default 10

% It will create one large matrix of trials from cells which were exposed
% to the chosen drug at the chosen current and give each trial a label of
% either 'Control' or 'Drug'. It will then partition the matrix into cross_val_n
% partitions. Cross_val_n - 1 partitions will be used as a trianing set and
% the remainder as the test set. Cross_val_n SVMs will be trained. The
% comparison of original labels and predicted labels (length of original
% matrix of trials) will be calculated and included in the output - as well
% as all of the trained models. 

% The output struct includes the models, the labels, the scores.
function rslt = drug_svm_jackknife( mfs, drug_selection, current_selection, crossval_n )
    
    if nargin < 4, crossval_n = 10; end

    % Get giant matrix of trials for all cells with drug selection and
    % current selection
    
    [svm_mat, drug_labels] = multicell_mat( mfs, drug_selection, current_selection );
    
    assignin( 'base', 'svm_mat', svm_mat);
    assignin( 'base', 'drug_labels', drug_labels);
    
    % Separate data into crossvalidation partitions 
    partition_struct = partition_svm_mat( svm_mat, drug_labels, cross_val_n );
    
    % Generating training and test sets from the 10 partions
    svm_struct = gen_crossval_svm_struct( partition_struct );
    
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

        svm_struct(i).predicted_labels = label;
        svm_struct(i).score = score;
        svm_struct(i).Model = SVMModel;
        
    end
    
    % Perform sign test on 'left out' trials
    rslt
      
end

% Separate a struct containing n trial matrices into n training sets and
% test sets
function svm_struct = gen_crossval_svm_struct( partition_struct )

    svm_struct = struct;

    for left_out_idx = 1:length(partition_struct)
       
        svm_struct(left_out_idx).train_svm_mat = vertcat(partition_struct( ~ left_out_idx ).svm_mat);
        svm_struct(left_out_idx).train_drug_labels = vertcat(partition_struct( ~ left_out_idx ).orig_labels);
        
        svm_struct(left_out_idx).test_svm_mat = vertcat(partition_struct( left_out_idx ).svm_mat);
        svm_struct(left_out_idx).test_drug_labels = vertcat(partition_struct( left_out_idx ).orig_labels);

    end

end


% Break a matrix of trials and a matrix of data labels into cross_val_n
% number of partitions. Return them in a struct.
function partition_struct = partition_svm_mat( svm_mat, drug_labels, cross_val_n )

    partition_struct = struct;

    num_trials = size(svm_mat,1);
    partition_idxs = crossvalind('Kfold', num_trials, crossval_n);

    for i = 1:max(partition_idxs)

        partition_struct(i).svm_mat = svm_mat(partition_idxs == i,:);
        partition_struct(i).orig_labels = drug_labels(partition_idxs == i,:);
        
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

% Loads the chosen processed file
function data_struct = load_processed_file( sub_direc, processed_fname )

    if strcmp( comp_mac_address, 'iMac'), save_direc = '/Users/Adrienne/Documents/MATLAB/Cell_Analysis/Processed';
    else save_direc = '/Users/eddi/Documents/MATLAB/Cell_Analysis/Processed';
    end
     
    ffp = fullfile( save_direc, sub_direc, processed_fname );
    load( ffp, 'data_struct' );
end

% Determines which computer the code is being run on - the iMac or the
% MacBook
function computer = comp_mac_address()
    localhost = java.net.InetAddress.getLocalHost;
    networkinterface = java.net.NetworkInterface.getByInetAddress(localhost);
    macaddress = typecast(networkinterface.getHardwareAddress, 'uint8');
    
    if min(macaddress == [56;201;134;2;215;75]) == 1
        computer = 'iMac';
    else, computer = 'MacBook';
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






