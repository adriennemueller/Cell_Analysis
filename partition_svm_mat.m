% Break a matrix of trials and a matrix of data labels into cross_val_n
% number of partitions. Return them in a struct.
function partition_struct = partition_svm_mat( svm_mat, labels, crossval_n )

    partition_struct = struct;

    num_trials = size(svm_mat,2);

    partition_idxs = crossvalind('Kfold', num_trials, crossval_n);
    for i = 1:max(partition_idxs)
        partition_struct(i).svm_mat = svm_mat(:, partition_idxs == i);
        partition_struct(i).orig_labels = labels( partition_idxs == i );
    end
end