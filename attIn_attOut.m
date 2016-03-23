
% Attend In - Attend Out

function rslt = attIn_attOut( trials )

    % Remove error trials
    correct_trial_idx = find([trials.trial_error] == 0);
    correct_trials = trials(correct_trial_idx);
    
    % Adjust Theta
    if length(unique([correct_trials.theta])) > 9
        correct_trials = adjust_theta( correct_trials );
    end
       
    % Count spikes in attentional window for each trial. Get list of numbers
    % for attend in, list of numbers for attend out.
    window = [121 126]; % Event codes for attentional window (Cue on, before targets blank and flip);
    correct_trials = count_spikes( correct_trials, window );
    
    % For each direction - get idxs for attend in and attend out, when drug
    % is present and when drug is absent
    idx_struct = segregate( correct_trials );
    
    % Get D' for result of this
    dmat = gen_dprime_struct( idx_struct, correct_trials );
    
    % Get Attend In/Out Drug On/Off SDen averages and SEs
    sdens = get_attend_sdens( idx_struct, correct_trials, 121 );
    sden_summs =  get_trial_sum( sdens );
    
    rslt.dmat = dmat;
    rslt.sdens = sdens;
    rslt.sden_summs = sden_summs;
    
    % MAKE A FUNCTION TO RUN THIS FUNCTION AND ADD THE OUTPUT TO THE
    % SESSION_STRUCT AND SAVE IT OUT?

end

function rslt = get_trial_sum( trials )
    
    rslt = struct;

    for i = 1:length(trials); % eight directions
        rslt(i).dOaI_millis = trials(i).drugOff_attIn(1,:);
        dOaI_sdens  = trials(i).drugOff_attIn(2:end,:);
        [ts ms] = size(dOaI_sdens);
        rslt(i).dOaI_avg    = mean(dOaI_sdens, 1);
        rslt(i).dOaI_ste    = std(dOaI_sdens) / sqrt(ts); % Make sure std is in right dimension
        
        rslt(i).dOaO_millis = trials(i).drugOff_attOut(1,:);
        dOaO_sdens  = trials(i).drugOff_attOut(2:end,:);
        [ts ms] = size(dOaO_sdens);
        rslt(i).dOaO_avg    = mean(dOaO_sdens, 1);
        rslt(i).dOaO_ste    = std(dOaO_sdens) / sqrt(ts); % Make sure std is in right dimension

        rslt(i).dNaI_millis = trials(i).drugOn_attIn(1,:);
        dNaI_sdens  = trials(i).drugOn_attIn(2:end,:);
        [ts ms] = size(dNaI_sdens);
        rslt(i).dNaI_avg    = mean(dNaI_sdens, 1);
        rslt(i).dNaI_ste    = std(dNaI_sdens) / sqrt(ts); % Make sure std is in right dimension

        rslt(i).dNaO_millis = trials(i).drugOn_attOut(1,:);
        dNaO_sdens  = trials(i).drugOn_attOut(2:end,:);
        [ts ms] = size(dNaO_sdens);
        rslt(i).dNaO_avg    = mean(dNaO_sdens, 1);
        rslt(i).dNaO_ste    = std(dNaO_sdens) / sqrt(ts); % Make sure std is in right dimension
        
    end
end


function sdens = get_attend_sdens( idx_struct, trials, align_idx )

    drugOff_attIn = [];
    drugOff_attOut = [];
    drugOn_attIn = [];
    drugOn_attOut = [];

   for i = 1:length(unique([idx_struct.theta])); % For the 8 directions
        
        theta = map_direction(i);
        
        % Find the Drug Off, Attend In trials and
        row = find(([idx_struct.theta] == theta) & ([idx_struct.drug] == 0) & ([idx_struct.attend] == 1));
        sub_trials = trials(idx_struct(row).idxs);
        sdens(i).drugOff_attIn = align_trials( sub_trials, align_idx );
        
        % Same for Drug Off, Attend Out
        row = find(([idx_struct.theta] == theta) & ([idx_struct.drug] == 0) & ([idx_struct.attend] == 0));
        sub_trials = trials(idx_struct(row).idxs);
        sdens(i).drugOff_attOut = align_trials( sub_trials, align_idx );
        
        % Find the Drug On, Attend In trials and
        row = find(([idx_struct.theta] == theta) & ([idx_struct.drug] == 1) & ([idx_struct.attend] == 1));
        sub_trials = trials(idx_struct(row).idxs);
        sdens(i).drugOn_attIn = align_trials( sub_trials, align_idx );
        
        % Same for Drug On, Attend Out
        row = find(([idx_struct.theta] == theta) & ([idx_struct.drug] == 1) & ([idx_struct.attend] == 0));
        sub_trials = trials(idx_struct(row).idxs);
        sdens(i).drugOn_attOut = align_trials( sub_trials, align_idx );
      
        sdens(i).theta = theta;
    end
    
    % Get averaged sdens, + Ses for the four matrices, for the 8
    % directions
    
    % Save out a struct with the trial subset sdens, ses, (+ spikes? )

end

% Align all of them so the zeros are at cue_onset and crop to appropriate lengths.
% They need to be of uniform length so that the average spike density
% isn't calculated from just a subset of trials at the tails.
% Put into a matrix of set length in ms, and n trials long

%%% THE MILLIS ARE IN THE FIRST ROW!!!!!
function [spike_mat, aligned_millis] = align_trials( trials, align_code )
    window = [-500, 750];

    aligned_millis = window(1) : window(2);
    spike_mat = zeros(length(trials) + 1, window(2) - window(1) +1);
    spike_mat(1,:) = aligned_millis;
    
    for i = 1:length(trials)
        align_time = trials(i).code_times( [trials(i).event_codes] == align_code );
        millis = trials(i).spike_millis;
        sdens  = trials(i).sden150;
        
        align_idx = find(millis == align_time);
        spike_mat(i+1,:) = sdens( (align_idx + window(1)) : (align_idx + window(2)) );
    end
end


% gen_dprime_struct takes an index-struct (containing lists of trial
% indexes for different conditions) and a trial struct (containing the
% behavioral and spike info for the conditions themselves. It outputs a
% matrix of the AttendIn vs AttendOut d' values for drug off, and drug on trials
% separately; for each direction. So the first 2 lines of output would be:
% 0  0.4 0.2
% 45 1.0 0.5
% For Direction DrugOff DrugOn
function rslt = gen_dprime_struct( idx_struct, trials )   
    rslt = [];

    for i = 1:length(unique([idx_struct.theta])); % For the 8 directions
        
        theta = map_direction(i);
        
        % Find the Drug Off, Attend In trials and make a vector of the
        % number of spikes for them.
        row = find(([idx_struct.theta] == theta) & ([idx_struct.drug] == 0) & ([idx_struct.attend] == 1));
        drugOff_AttIn_n_spikes = [trials(idx_struct(row).idxs).n_spikes];
        
        % Same for Drug Off, Attend Out
        row = find(([idx_struct.theta] == theta) & ([idx_struct.drug] == 0) & ([idx_struct.attend] == 0));
        drugOff_AttOut_n_spikes = [trials(idx_struct(row).idxs).n_spikes];
        
        % Calculate the D' for Drug Off Attend In vs Attend Out.
        drugOff_dprime = d_prime( drugOff_AttIn_n_spikes, drugOff_AttOut_n_spikes );
        
        % Find the Drug On, Attend In trials and make a vector of the
        % number of spikes for them.        
        row = find(([idx_struct.theta] == theta) & ([idx_struct.drug] == 1) & ([idx_struct.attend] == 1));
        drugOn_AttIn_n_spikes = [trials(idx_struct(row).idxs).n_spikes];
        
        % Same for Drug On, Attend Out
        row = find(([idx_struct.theta] == theta) & ([idx_struct.drug] == 1) & ([idx_struct.attend] == 0));
        drugOn_AttOut_n_spikes = [trials(idx_struct(row).idxs).n_spikes];
        
        % Calculate the D' for Drug On Attend In vs Attend Out.
        drugOn_dprime = d_prime( drugOn_AttIn_n_spikes, drugOn_AttOut_n_spikes );
        
        % Append the result.
        rslt = [rslt; theta drugOff_dprime drugOn_dprime];
    end
end



% count_spikes counts the number of spikes in a given window (eg cue,
% blank, etc) for each trial and appends it to the trial_structure that was
% input. The window is given by two behavioral codes. eg:
% [121 126] - Event codes for attentional window [Cue On, Targets Blank];
function trials = count_spikes( trials, window )

    for i = 1:length(trials)
        
        % Find the times during the trial at which the window of interest (eg cue on, blank, etc) 
        % starts and ends.
        start_idx = trials(i).code_times( [trials(i).event_codes] == window(1) );
        end_idx   = trials(i).code_times( [trials(i).event_codes] == window(2) );
        
        % Get the start and end indexes for the window we want to count the
        % spikes in, but finding the index at which they occur in the
        % spike_millis vector for the trial. (The spikes vector is always
        % longer than the behavioral data vectors, but they are both
        % aligned to trial onset == 0.)
        spike_start_idx = find(trials(i).spike_millis == start_idx );
        spike_end_idx   = find(trials(i).spike_millis == end_idx );
        
        % Add up the number of spikes in the window
        n_spikes = sum(trials(i).spikes(spike_start_idx:spike_end_idx));
        
        % This will be ambiguous if there are ever other spike counts in the
        % struct.
        trials(i).n_spikes = n_spikes;
    end
    

end


% segregate takes a structure of trials and makes a new index_struct which
% lists the indexes of trials matching particular conditions, eg:
% Direction 45�, Drug On, Attend Out
% Direction 180�, Drug Off, Attend In
% Currently allows for only two drug conditions: drug off (-), drug on (+).
function idx_struct = segregate( trials )
    idx_struct = struct;
    tmp_struct = struct;
    
    % For each possible of 8 directions
    for direction = 1:8
        
        % Get what theta that direction matches up to.
        tmp_struct.theta = map_direction(direction);
        
        % And for conditions whether the drug is on or off.
        for drug = 0:1
            
            tmp_struct.drug = drug;
            
            % Currently assuming retaining a drug is a negative sign and
            % ejecting it is a positive sign. Also not allowing for
            % different values of injected drug. Should amend to allow.
            % (Unique?)
            if drug == 0;
                drug_sign = -1;
            else
                drug_sign = 1;
            end
            
            % And for either attend-in or attend-out conditions.
            for attend = 0:1
                
                tmp_struct.attend = attend;
                
                % Attend-out conditions are when the attention is in the
                % oposite direction.
                if attend == 0
                    attend_direction = reversed(direction);
                else
                    attend_direction = direction;
                end
                
                % Get all the trials that match
                tmp_struct.idxs = find( (sign([trials.drug]) == drug_sign) & ...
                                                   ([trials.theta] == map_direction(attend_direction)));
                                               
                if isempty(fieldnames(idx_struct))
                    idx_struct = tmp_struct;
                else
                    idx_struct = [idx_struct tmp_struct];
                end
            end
        end
    end
end



% reversed gives you the opposite direction (1-8) to the one you input. eg
% 1->5, 3->7, 6->2, etc.
function rslt = reversed(direction)
    rslt = mod((direction + 3), 8)+1;
end


% map_direction takes a value of 1-8 and maps it to a direction from 0-315.
function rslt = map_direction( dir_theta )
    theta_vec = [0 45 90 135 180 225 270 315];
    rslt = theta_vec(dir_theta); % TEST ME
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

%     NGroups = 8;
%     groupWidth = 2 * pi / NGroups;
%     thetas = [trials.theta];
% 
%     for group = 1:NGroups
%         theta = group * groupWidth;
% 
%         
%         if close_angles(thetas(i), theta, groupWidth)
% 
%         end
%     end
% 
% function rslt = close_angles(a, b, t)
%     rslt = 0;
%     d = abs(a - b);
%     if d < t
%         rslt = 1;
%      elseif 2 * pi - d < t
%          rslt = 1;
%     end
% end
