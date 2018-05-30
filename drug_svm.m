
function rslt = drug_svm( mfs, drug_selection, current_selection )


    % Get giant matrix of trials for all cells with drug selection and
    % current selection
    
    [svm_mat, drug_labels] = multicell_mat( mfs, drug_selection, current_selection );
    
    assignin( 'base', 'svm_mat', svm_mat);
    assignin( 'base', 'drug_labels', drug_labels);

    % Separate into training sets and test set
    [train_idxs, test_idxs] = separate_trials(size(svm_mat, 1), 0.8);
    
    train_svm_mat = svm_mat( train_idxs, : ); train_drug_labels = drug_labels( train_idxs );
    test_svm_mat  = svm_mat( test_idxs, : );  test_drug_labels = drug_labels(  test_idxs );
    
    
    
    % Train and crossvalidate
    disp( 'Training SVM...' );
    
    drug_label_string = strcat( drug_selection, {'_'}, num2str(current_selection) );
    SVMModel = fitcsvm(train_svm_mat, train_drug_labels, 'KernelFunction','rbf', 'Standardize',true,'ClassNames',{'Control',drug_label_string{1}});SVMModel = fitcsvm(train_svm_mat, train_drug_labels, 'KernelFunction','rbf', 'Standardize',true,'ClassNames',{'Control','SCH23390_50'});
    %SVMModel = fitclinear(train_svm_mat, train_drug_labels, 'ClassNames',{'Control', drug_label_string{1} });
    
    disp( 'Crossvalidating SVM...' );
    CVSVMModel = crossval(SVMModel);
    classLoss = kfoldLoss(CVSVMModel) % The output of this is the generalization rate.
 
    
    % Test
    [label,score] = predict( SVMModel, test_svm_mat );
    ScoreSVMModel = fitPosterior( SVMModel, train_svm_mat, train_drug_labels)
    ScoreTransform = ScoreSVMModel.ScoreTransform
    
    
    
end

function [train_idxs, test_idxs] = separate_trials( mat_length, train_prop )
    idx_vec = 1:mat_length;
    perm_idxs = randperm(numel( idx_vec ));

    train_length = round(mat_length * train_prop);

    train_idxs = perm_idxs(1:train_length);
    test_idxs  = perm_idxs(train_length+1:mat_length);
    
end


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


function data_struct = load_processed_file( sub_direc, processed_fname )

    if strcmp( comp_mac_address, 'iMac'), save_direc = '/Users/Adrienne/Documents/MATLAB/Cell_Analysis/Processed';
    else save_direc = '/Users/eddi/Documents/MATLAB/Cell_Analysis/Processed';
    end
     
    ffp = fullfile( save_direc, sub_direc, processed_fname );
    load( ffp, 'data_struct' );
end

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






