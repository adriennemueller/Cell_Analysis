
function data_mat = factored_data_mat( data_struct, retain_current, eject_current, window_str, contrast_flag )
    
    currents = [retain_current, eject_current];
    spikes = []; theta = []; target_change = []; drug = []; contrast = [];

    for i = 1:length(currents)
    
        % Loop through datastruct, creating a matrix based on the epoch 
        [new_spikes, contrasts, events] = filtered_windowed_spikemat( data_struct, currents(i), window_str, [], [], contrast_flag );
        
        %%% SHAVE THE VINS HERE %%% CHECK THIS HORZCAT!
        [spikes, new_spikes, ~] = shave_bins( spikes, new_spikes, new_spikes );
        spikes = horzcat( spikes, new_spikes );
        
        %spikes        = horzcat( spikes, sum( ctrl_spike_mat, 1 ));
        theta         = horzcat( theta, [events.theta{:}] );
        target_change = horzcat( target_change, [events.target_change{:}] );
        drug          = horzcat( drug, repmat( currents(i), 1, length(events.theta)) );
        contrast      = horzcat( contrast, contrasts );
        %attend        =  
    end
    
    % Just drug as factor - 'Fixation' Window
    if ismember( window_str, {'fixation', 'attend_fixation', 'wm_fixation'} )
        merged_mat = {drug};
        factors_list = {'drug'};
    % If direction also a factor
    elseif sum( strcmp( window_str, {'visual', 'attend_visual', 'wm_visual', 'wm_delay', 'wm_response', 'attend_reward', 'wm_reward','reward'} ) )
        merged_mat = {drug, theta};
        factors_list = { 'drug', 'theta' };
    % IF ATTEND A FACTOR - DOESN'T WORK YET
    elseif sum( strcmp( window_str, {'attend', 'attContrast', 'blank'} ) )
%         merged_mat = {spikes, drug, theta, attend}; FIX THIS
%         factors_list = { 'drug', 'theta', 'attend' }; FIX THIS
        merged_mat = {drug, theta};
        factors_list = { 'drug', 'theta' };
    % If whether target flips is a factor
    elseif strcmp( window_str, 'post_blank' )
%         merged_mat = {spikes, drug, theta, attend, target_change}; FIX THIS
%         factors_list = { 'drug', 'theta', 'attend', 'target_change' }; FIX ThiS
        merged_mat = {drug, theta, target_change}; 
        factors_list = {'drug', 'theta', 'target_change'};
    end        

    % Add Contrast Flag if applicable
    if contrast_flag
        merged_mat{end+1} = cell2mat(contrast);
        factors_list{end+1} = 'contrast';
    end
    
    data_mat.spikes  = spikes; %merged_mat{ 1, :};
    data_mat.factors = merged_mat;
    data_mat.factors_strings = factors_list;
    
end