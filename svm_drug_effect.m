
% paradigm: either 'Attention' or 'WM'

% svm_parameters will include
%  theta
%  attend_direction
%  decision
%  [rewarded] future plan
function svm_plot = svm_drug_effect( mfs, drug, current, paradigm, crossval_n )

    if nargin < 5, crossval_n = 3; end
    
    jackknife_flag = 0;
    
    % For a given paradigm
    if strcmp( paradigm, 'Attention' )
        window_str_list = { 'attend_visual', 'attend', 'blank', 'post_blank', 'attend_reward'};
    elseif strcmp( paradigm, 'WM' )
        window_str_list = { 'wm_visual', 'wm_delay', 'wm_response', 'wm_reward' };
    else 
        disp( 'Unknown paradigm' ); return;
    end
    
    epoch_struct = struct;
    
    % Loop through all epochs (window_str_list) and get the cells and the
    % appropriate labels for each of the svm parameters (factors)
    % associated with each trial
    for i = 1:length(window_str_list)
    
        % cell_struct contains the field 'spikes' and the field 'factors' and
        % is n units longs
        cell_struct = filter_cells( mfs, drug, current, window_str_list(i) );
                
        % Loop through all relevant svm_parameters 
        for j = 1:length(cell_struct) %%% CHECK THIS LENGTH 
        
            for k = 2:length( cell_struct(j).factors_strings )
                
                svm_parameter = cell_struct(j).factors_strings{k};
                svm_parameter_idx = k;
                drug_factor_idx = 1; % First factor is always drug
                factor_labels = cell_struct(j).factors{ 1, svm_parameter_idx };

                
                % Separate into drug and control trials
                control_idxs = find( cell_struct(j).factors{ drug_factor_idx } == -15 );
                control_trials = cell_struct(j).spikes( :, control_idxs );
                control_factor_labels = factor_labels( control_idxs );
                
                drug_idxs = find( cell_struct(j).factors{ drug_factor_idx } ~= -15 );
                drug_trials = cell_struct(j).spikes( :, drug_idxs );
                drug_factor_labels = factor_labels( drug_idxs );
                
                % train and test crossval_n number of SVMs each for drug trials 
                % and control trials separately
                
                if jackknife_flag
                    control_partition_struct = partition_svm_mat( control_trials, control_factor_labels, crossval_n );
                    contrl_struct = gen_jkCrossval_svm_struct( control_partition_struct );
                    
                    drug_partition_struct = partition_svm_mat( drug_trials, drug_factor_labels, crossval_n );
                    drug_struct = gen_jkCrossval_svm_struct( drug_partition_struct );
                else
                    contrl_struct = gen_randCrossval_svm_struct( control_trials, control_factor_labels, crossval_n );
                    drug_struct = gen_randCrossval_svm_struct( drug_trials, drug_factor_labels, crossval_n );    
                end

                control_CVSVMModel = train_test_svm( contrl_struct, 0 ); % 0 for not scrambling
                drug_CVSVMModel    = train_test_svm( drug_struct, 0 ); % 0 for not scrambling
           
                % Output the difference between the averages for this
                % cell
                epoch_struct = append_to_struct_field( ...
                    epoch_struct, ...
                    mean( [control_CVSVMModel.perc_corr] ), ...
                    window_str_list{i}, svm_parameter, 'ctrl_avg_perc_correct');
                epoch_struct = append_to_struct_field( ...
                    epoch_struct, ...
                    mean( [drug_CVSVMModel.perc_corr] ), ...
                    window_str_list{i}, svm_parameter, 'drug_avg_perc_correct');
                              
            end
        end
    end
    
    assignin( 'base', 'epoch_struct', epoch_struct );
    
    svm_plot = plot_svm_drug_effect_histos( epoch_struct );
end

%%% MOVE PLOT SVM DRUG EFFECT HISTOS BACK OVER HERE WHEN DONE TESTING


function S = append_to_struct_field( S, value, varargin )
    if isempty(varargin)
        S = [S value];
        return
    end
    if ~isfield(S, varargin{1})
        if nargin == 3
            S.(varargin{1}) = [];
        else
            S.(varargin{1}) = struct();
        end
    end
    S.(varargin{1}) = append_to_struct_field(S.(varargin{1}), value, varargin{2:end});
end

%
function cell_struct = filter_cells( mfs, drug, current, window_str )

    contrast_flag = 0; % Could be changed later.

    cell_struct = struct;

    for i = 1:length(mfs.session)
        if strcmp(mfs.session(i).drug, drug)
            for j = 1:length(mfs.session(i).processed_files)
                if ismember(current, mfs.session(i).currents{j})
                    
                    fname = mfs.session(i).processed_files{j};
                    disp(strcat('Adding: ', {' '},  fname));
        
                    data_struct = load_processed_file( mfs.session(i).sub_direc, fname );
                    
                    factored_mat = factored_data_mat( data_struct, -15, current, window_str, contrast_flag );
                    if ~isfield( cell_struct, 'spikes' ), idx = 0;
                    else, idx = size(cell_struct, 2);
                    end
                    
                    cell_struct(idx+1).spikes = factored_mat.spikes;
                    cell_struct(idx+1).factors = factored_mat.factors;
                    cell_struct(idx+1).factors_strings = factored_mat.factors_strings;
      
                end
            end
        end
    end


end


% Break data into training and test sets. This will probably have a struct
% output actually
function [training_data test_data] = partition_data_for_svm( svm_input_data, jackknife_flag )


end

% Train SVM and test it and then have it spit out its performance
function SVM_performance = run_test_SVM( training_data, test_data )

end
