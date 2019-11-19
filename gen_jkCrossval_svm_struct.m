
% Separate a struct containing n trial matrices into n training sets and
% test sets
function svm_struct = gen_jkCrossval_svm_struct( partition_struct )

    svm_struct = struct;

    for left_out_idx = 1:length(partition_struct)
       
        curr_svm_struct = partition_struct;
        curr_svm_struct(left_out_idx) = [];
        
        svm_struct(left_out_idx).train_svm_mat = horzcat(curr_svm_struct.svm_mat);
        svm_struct(left_out_idx).train_labels = horzcat(curr_svm_struct.orig_labels);
        
        svm_struct(left_out_idx).test_svm_mat = partition_struct(left_out_idx).svm_mat;
        svm_struct(left_out_idx).test_labels = partition_struct(left_out_idx).orig_labels;

    end

end