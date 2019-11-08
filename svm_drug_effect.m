
% HAVE a separate function combining what already have for just drug/nodrug
% comparison

% Function to create data formatted to do drug vs no drug svm comparison

% svm_parameter:
%  theta
%  attend_direction
%  decision
%  [rewarded] future plan
function svm_plot = svm_drug_effect( mfs, drug, current, svm_parameter ) % Include paradigm?

    % Get the subset of cells with relevant drug and current

    % Loop through all relevant cells
    
    % For each epoch relevant to that svm_parameter
   
        % train and test an SVM
        
        % Get multiple estimates of SVM performance for a given cell

    % Once have gone through all cells; have a distribution of performances and
    % compare to chance. (With e.g. an anova? t-test)

end


% Break data into training and test sets. This will probably have a struct
% output actually
function [training_data test_data] = partition_data_for_svm( svm_input_data, jackknife_flag )


end

% Train SVM and test it and then have it spit out its performance
function SVM_performance = run_test_SVM( training_data, test_data )

end
