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
