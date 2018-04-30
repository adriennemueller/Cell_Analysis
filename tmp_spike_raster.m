% 120 - Fixation Spot On
% 124 - Targets On



function rslt = tmp_spike_raster( full_trials, currents )

% If only one current was applied, return an empty rslt.
    if length(currents) == 1
        rslt = [];
    end
    
    % Go through all different ejected currents in the file. (currents(1)
    % will be the retain current.
    for i = 2:length(currents)
        
        % Get only retain trials and trials with that current
        eject_current = currents(i); retain_current = currents(1);
        trials = segregate_by_current( full_trials, retain_current, eject_current );

        % Remove error trials
        correct_trial_idx = find([trials.trial_error] == 0);
        correct_trials = trials(correct_trial_idx);

        % Adjust Theta
        if length(unique([correct_trials.theta])) > 9
            correct_trials = adjust_theta( correct_trials );
        end

        % Count spikes in attentional window for each trial. Get list of numbers
        % for attend in, list of numbers for attend out.
        correct_trials = remove_probe_trials( correct_trials );
        
        
        % Make a new Figure
        % Go through each direction
        
        % `[out,idx] = sort([14 8 91 19])
        drug_sort = [correct_trials.drug];
        [tmp, drug_sort_idxes ] = sort(drug_sort);
        correct_trials = correct_trials(drug_sort_idxes);

        direc_sort = [correct_trials.theta];
        [tmp, direc_idxs ] = sort(direc_sort);
        correct_trials = correct_trials(direc_idxs);

        
        
        
            % Sort into Attend In and Attend Out, Drug On, Drug Off
            
            % Plot Four Rasters in a Stack with different Colors
            
            % Plot Spike Dens underneath
        
        align_code = 124; % Targets On
        [spike_mat, aligned_millis] = align_trials( correct_trials, align_code );

        spikes_only = spike_mat(2:end, :);
        spikes_only = logical(spikes_only);
        figure();
        
        
        
        [xs, ys] = plotSpikeRaster( spikes_only, 'PlotType', 'scatter' );
        
        
%                    plot([501:1000], C25attNOut, 'r:');
%            plot([501:1000], C25attNIn, 'r');
 
        % Should maybe use text to make labels centered.
        %Redo this properly
        offset = 300;
        hold_time = 300; 
        cue_time = 300; %%% SOMETIMES 500!!!
        blank_time = 300;

        set(gca,'Xtick',[offset, offset+hold_time, offset+hold_time+cue_time, offset+hold_time+cue_time+blank_time],'XTickLabel',{ ['Targets \newline On'], ['Cue \newline On'], ['Targ \newline Flip \newline Off'], ['Targ \newline Flip \newline On'] });
        
        
        
        % For each direction - get idxs for attend in and attend out, when drug
        % is present and when drug is absent
 %       idx_struct = segregate( correct_trials );
 %       rslt = idx_struct; %%% TMP TMP TMP %%%

    end
    
    
    
    
end


% Make a new trials structure with only the retain trials and trials with
% ejection of current of interest
function trials = segregate_by_current( trials, retain_c, eject_c )
    valid_idxs = find( ([trials.drug] == retain_c) | ([trials.drug] == eject_c));
    trials = trials(valid_idxs);
end

% For sessions that have probe trials, this removes the probe trials.
function correct_trials = remove_probe_trials( correct_trials )
    indexes = [];

    % Can I do this without a for loop?
    for i = 1:length(correct_trials)
        if find( correct_trials(i).event_codes == 126 ) % Contains a blanked cue; therefore not a probe trial
            indexes = [indexes i];
        end
    end

    correct_trials = correct_trials(indexes);
end


% Event codes for attentional window (Cue on, before targets blank and flip);
% Need a function for this because I changed the event codes in Jan/Mar 2016   
function window = get_attend_window( correct_trials )
    if find(correct_trials(1).event_codes == 121);
        window = [121 126];
    else
        window = [133 126];
    end
end
        
% Align all of them so the zeros are at cue_onset and crop to appropriate lengths.
% They need to be of uniform length so that the average spike density
% isn't calculated from just a subset of trials at the tails.
% Put into a matrix of set length in ms, and n trials long

%%% THE MILLIS ARE IN THE FIRST ROW!!!!!
function [spike_mat, aligned_millis] = align_trials( trials, align_code )
    window = [-300, 1300];

    aligned_millis = window(1) : window(2);
    spike_mat = zeros(length(trials) + 1, window(2) - window(1) +1);
    spike_mat(1,:) = aligned_millis;
    
    for i = 1:length(trials)
        align_time = trials(i).code_times( [trials(i).event_codes] == align_code );
        millis = trials(i).spike_millis;
        spikes = trials(i).spikes;
        
        align_idx = find(millis == align_time);
        spike_mat(i+1,:) = spikes( (align_idx + window(1)) : (align_idx + window(2)) );
    end
end

% adjust_theta takes a trial-struct and changes the thetas to the closest
% '8-direction' theta. Non-ideal, but RFs are large enough that grouping
% them shouldn't change the firing too much. Also, collapsing 3-4
% directions down to 1 direction, if the RFs are small, will only make my
% signal worse.
% Crude method.
function adjusted_trials = adjust_theta( trials )

    for i = 1:length(trials)
        theta = round(trials(i).theta);
        if (theta <= 20) || (340 <= theta);
            trials(i).theta = 0;
        elseif (30 <= theta) && (theta <= 60)
            trials(i).theta = 45;
        elseif (70 <= theta) && (theta <= 110)
            trials(i).theta = 90;
        elseif (120 <= theta) && (theta <= 150)
            trials(i).theta = 135;
        elseif (160 <= theta) && (theta <= 200)
            trials(i).theta = 180;
        elseif (210 <= theta) && (theta <= 240)
            trials(i).theta = 225;
        elseif (250 <= theta) && (theta <= 290)
            trials(i).theta = 270;
        elseif (300 <= theta) && (theta <= 330)
            trials(i).theta = 315;
        end
    end
    
    adjusted_trials = trials;
end
