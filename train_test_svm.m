function CVSVMModel = train_test_svm( svm_struct, scramble_flag )

    CVSVMModel = struct;

    % Loop through svm_struct, generating an SVM for each train/test set
    % and predicting test values. Add predicted test labels to list
    % compared with real labels
    for i = 1:length(svm_struct)

        % Train SVM
        disp( strcat( 'Training SVM:', {' '}, num2str(i), '...' ) );
        train_svm_mat = svm_struct(i).train_svm_mat;
        train_labels = svm_struct(i).train_labels;
        test_svm_mat = svm_struct(i).test_svm_mat;
        test_labels = svm_struct(i).test_labels;
        
        if isa(train_labels, 'double')
            train_labels = string(train_labels);
            test_labels  = string(test_labels);
        end
        
        class_names = unique( train_labels );
        
        if length(class_names) == 2
            SVMModel = fitcsvm(train_svm_mat', train_labels, 'KernelFunction','rbf', 'Standardize',true,'ClassNames', class_names);
        else
            SVMModel = fitcecoc(train_svm_mat', train_labels, 'ClassNames', class_names);
        end
            
        % Test SVM
        [label,score] = predict( SVMModel, test_svm_mat' );
        perc_corr = compare_labels( label, test_labels' );

        CVSVMModel(i).perc_control = sum(strcmp(train_labels, 'Control')) / length(train_labels);
        CVSVMModel(i).Model = SVMModel;
        CVSVMModel(i).label = label; 
        CVSVMModel(i).score = score;
        CVSVMModel(i).perc_corr = perc_corr;

        if scramble_flag
            scramble_train_labels = train_labels(randperm(length(train_labels)));
            disp( strcat( 'Training Scrambled SVM:', {' '}, num2str(i), '...' ) );
            if length(class_names) == 2
                Scramble_SVMModel = fitcsvm(train_svm_mat', scramble_train_labels, 'KernelFunction','rbf', 'Standardize',true,'ClassNames', class_names);
            else
                Scramble_SVMModel = fitcecoc(train_svm_mat', scramble_train_labels, 'ClassNames', class_names);
            end    
            [scr_label, scr_score] = predict( Scramble_SVMModel, test_svm_mat' );
            scr_perc_corr = compare_labels( scr_label, test_labels' );

            CVSVMModel(i).Scamble_Model = Scramble_SVMModel;
            CVSVMModel(i).scr_label = scr_label; 
            CVSVMModel(i).scr_score = scr_score;
            CVSVMModel(i).scr_perc_corr = scr_perc_corr;
        end      
    end 
end

% compare_labels determines how many of the svm labels in a test set were
% correctly assigned
function perc_corr = compare_labels( svm_labels, orig_labels )
    correct = strcmp( svm_labels, orig_labels );
    perc_corr = sum(correct) / length(svm_labels);
end