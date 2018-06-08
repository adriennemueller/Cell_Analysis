
% Attend In - Attend Out

function rslt = attIn_attOut( full_trials, currents )

    

        % MAKE A FUNCTION TO RUN THIS FUNCTION AND ADD THE OUTPUT TO THE
        % SESSION_STRUCT AND SAVE IT OUT?
    end
end

%
function event_durs = get_event_durations( trial )

    e_times = trial.code_times;
    codes = trial.event_codes;
    
    fix_on   = get_time_matching_code( e_times, codes, 120);
    targ_on  = get_time_matching_code( e_times, codes, 124);
    targ_off = get_time_matching_code( e_times, codes, 126);
    targ_on2 = get_time_matching_code( e_times, codes, 128);
    
    if find(trial.event_codes == 121)
        cue_on = get_time_matching_code( e_times, codes, 121);
    else
        cue_on = get_time_matching_code( e_times, codes, 133);
    end

    event_durs.pre_cue = round((cue_on - targ_on) / 50) * 50;
    event_durs.cue     = round((targ_off - cue_on) / 50) * 50;
    event_durs.blank   = round((targ_on2 - targ_off) / 50) * 50;
end

function e_time = get_time_matching_code(times, codes, code)
    idx = find(codes == code);
    e_time = times(idx);
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
    if find(correct_trials(1).event_codes == 121)
        window = [121 126];
    else
        window = [133 126];
    end
end


function vis_pval = gen_vis_pval( idx_struct, trials )

    % Can combine AttendIn and AttendOut trials because identical during period
    % of initial stimulus presentation

    vis_pval = [];

    for i = 1:length(unique([idx_struct.theta])) % For the 8 directions
        
        theta = map_direction(i);
        
        % Find the Drug Off, Attend In trials and make a vector of the
        % number of spikes for them.
        row = find(([idx_struct.theta] == theta) & ([idx_struct.drug] == 0) & ([idx_struct.attend] == 1));
        drugOff_AttIn_n_spikes = [trials(idx_struct(row).idxs).visual_n_spikes];
        drugOff_AttIn_n_spikes_baseline = [trials(idx_struct(row).idxs).baseline_n_spikes];
        
        % Same for Drug Off, Attend Out
        row = find(([idx_struct.theta] == theta) & ([idx_struct.drug] == 0) & ([idx_struct.attend] == 0));
        drugOff_AttOut_n_spikes = [trials(idx_struct(row).idxs).visual_n_spikes];
        drugOff_AttOut_n_spikes_baseline = [trials(idx_struct(row).idxs).baseline_n_spikes];
        
        
        % Combine lists
        drugOff_vis_n_spikes = [drugOff_AttIn_n_spikes drugOff_AttOut_n_spikes];
        drugOff_base_n_spikes = [drugOff_AttIn_n_spikes_baseline drugOff_AttOut_n_spikes_baseline];
        
        % Calculate the Wilcox Rank Sum for Drug Off Visual Response vs Baseline.
        if ~ isempty(drugOff_vis_n_spikes) &&  ~ isempty(drugOff_base_n_spikes)
            drugOff_ranksum = ranksum( drugOff_vis_n_spikes, drugOff_base_n_spikes );
        else
            drugOff_ranksum = nan;
        end
        
        % Find the Drug On, Attend In trials and make a vector of the
        % number of spikes for them.        
        row = find(([idx_struct.theta] == theta) & ([idx_struct.drug] == 1) & ([idx_struct.attend] == 1));
        drugOn_AttIn_n_spikes = [trials(idx_struct(row).idxs).visual_n_spikes];
        drugOn_AttIn_n_spikes_baseline = [trials(idx_struct(row).idxs).baseline_n_spikes];
        
        % Same for Drug On, Attend Out
        row = find(([idx_struct.theta] == theta) & ([idx_struct.drug] == 1) & ([idx_struct.attend] == 0));
        drugOn_AttOut_n_spikes = [trials(idx_struct(row).idxs).visual_n_spikes];
        drugOn_AttOut_n_spikes_baseline = [trials(idx_struct(row).idxs).baseline_n_spikes];
        
        drugOn_vis_n_spikes = [drugOn_AttIn_n_spikes drugOn_AttOut_n_spikes];
        drugOn_base_n_spikes = [drugOn_AttIn_n_spikes_baseline drugOn_AttOut_n_spikes_baseline];
        
        % Calculate the D' and Wilcox Rank Sum for Drug On Attend In vs Attend Out.
        if ~ isempty(drugOn_vis_n_spikes) &&  ~ isempty(drugOn_base_n_spikes)
            drugOn_ranksum = ranksum( drugOn_vis_n_spikes, drugOn_base_n_spikes );
        else
            drugOn_ranksum = nan;
        end
        
        % Get Firing Rate for Drug Off and Drug On
        drugOff_VisAvg = mean(drugOff_vis_n_spikes)/0.25; % So bad. So, so bad.
        drugOn_VisAvg  = mean(drugOn_vis_n_spikes)/0.25;
        
        
        % Append the result.
        vis_pval = [vis_pval; theta drugOff_ranksum drugOn_ranksum drugOff_VisAvg drugOn_VisAvg];
    end
    
    
    
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

        %%% ATT In and ATT Out Combined %%%
        rslt(i).dN_millis = trials(i).drugOn_attIn(1,:);
        dN_sdens = vertcat(dNaI_sdens, dNaO_sdens);
        [ts ms] = size(dN_sdens);
        rslt(i).dN_avg    = mean(dN_sdens, 1);
        rslt(i).dN_ste    = std(dN_sdens) / sqrt(ts); % Make sure std is in right dimension

        rslt(i).dO_millis = trials(i).drugOn_attOut(1,:);
        dO_sdens  = vertcat(dOaI_sdens, dOaO_sdens);
        [ts ms] = size(dO_sdens);
        rslt(i).dO_avg    = mean(dO_sdens, 1);
        rslt(i).dO_ste    = std(dO_sdens) / sqrt(ts); % Make sure std is in right dimension
        
        
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
        drugOff_AttIn_n_spikes = [trials(idx_struct(row).idxs).attend_n_spikes];
        
        % Same for Drug Off, Attend Out
        row = find(([idx_struct.theta] == theta) & ([idx_struct.drug] == 0) & ([idx_struct.attend] == 0));
        drugOff_AttOut_n_spikes = [trials(idx_struct(row).idxs).attend_n_spikes];
        
        % Calculate the D' and Wilcox Rank Sum for Drug Off Attend In vs Attend Out.
        drugOff_dprime  = d_prime( drugOff_AttIn_n_spikes, drugOff_AttOut_n_spikes );
        if ~ isempty(drugOff_AttIn_n_spikes) &&  ~ isempty(drugOff_AttOut_n_spikes)
            drugOff_ranksum = ranksum( drugOff_AttIn_n_spikes, drugOff_AttOut_n_spikes );
        else
            drugOff_ranksum = nan;
        
        end
        
        % Find the Drug On, Attend In trials and make a vector of the
        % number of spikes for them.        
        row = find(([idx_struct.theta] == theta) & ([idx_struct.drug] == 1) & ([idx_struct.attend] == 1));
        drugOn_AttIn_n_spikes = [trials(idx_struct(row).idxs).attend_n_spikes];
        
        % Same for Drug On, Attend Out
        row = find(([idx_struct.theta] == theta) & ([idx_struct.drug] == 1) & ([idx_struct.attend] == 0));
        drugOn_AttOut_n_spikes = [trials(idx_struct(row).idxs).attend_n_spikes];
        
        % Calculate the D' and Wilcox Rank Sum for Drug On Attend In vs Attend Out.
        drugOn_dprime  = d_prime( drugOn_AttIn_n_spikes, drugOn_AttOut_n_spikes );
        if ~ isempty(drugOn_AttIn_n_spikes) &&  ~ isempty(drugOn_AttOut_n_spikes)
            drugOn_ranksum = ranksum( drugOn_AttIn_n_spikes, drugOn_AttOut_n_spikes );
        else
            drugOn_ranksum = nan;
        end
        
        drugOff_avgFR = mean(horzcat(drugOff_AttIn_n_spikes, drugOff_AttOut_n_spikes)) / 0.3; % 300ms lat cue attend window - to get FR
        drugOn_avgFR =  mean(horzcat(drugOn_AttIn_n_spikes, drugOn_AttOut_n_spikes)) / 0.3;
        
        
        % Append the result.
        rslt = [rslt; theta drugOff_dprime drugOn_dprime drugOff_ranksum drugOn_ranksum drugOff_avgFR drugOn_avgFR];
    end
end


function rslt = gen_anova_struct( idx_struct, trials )
    data_vec = [];
    direction = [];
    drug = [];
    attend = [];
    
    for i = 1:length(unique([idx_struct.theta])); % For the 8 directions % Only 4 dir b/c copied att in/out.
        
        theta = map_direction(i);
        if (theta >= 180), continue, end
        
        
        % Find the Drug Off, Attend In trials and make a vector of the
        % number of spikes for them.
        row = find(([idx_struct.theta] == theta) & ([idx_struct.drug] == 0) & ([idx_struct.attend] == 1));
        drugOff_AttIn_n_spikes = [trials(idx_struct(row).idxs).attend_n_spikes];
        
        % Same for Drug Off, Attend Out
        row = find(([idx_struct.theta] == theta) & ([idx_struct.drug] == 0) & ([idx_struct.attend] == 0));
        drugOff_AttOut_n_spikes = [trials(idx_struct(row).idxs).attend_n_spikes];
        
        
        % Find the Drug On, Attend In trials and make a vector of the
        % number of spikes for them.        
        row = find(([idx_struct.theta] == theta) & ([idx_struct.drug] == 1) & ([idx_struct.attend] == 1));
        drugOn_AttIn_n_spikes = [trials(idx_struct(row).idxs).attend_n_spikes];
        
        % Same for Drug On, Attend Out
        row = find(([idx_struct.theta] == theta) & ([idx_struct.drug] == 1) & ([idx_struct.attend] == 0));
        drugOn_AttOut_n_spikes = [trials(idx_struct(row).idxs).attend_n_spikes];
        
        % Append the result.
        
        l_dOff_AI = length(drugOff_AttIn_n_spikes);
        l_dOff_AO = length(drugOff_AttOut_n_spikes);
        l_dOn_AI  = length(drugOn_AttIn_n_spikes);
        l_dOn_AO  = length(drugOn_AttOut_n_spikes);
        
        data_vec = [data_vec drugOff_AttIn_n_spikes drugOff_AttOut_n_spikes drugOn_AttIn_n_spikes drugOn_AttOut_n_spikes];
        direction = [ direction repmat(theta, 1, l_dOff_AI) repmat(theta, 1, l_dOff_AO) repmat(theta, 1, l_dOn_AI) repmat(theta, 1, l_dOn_AO) ];
        drug = [ drug zeros(1, l_dOff_AI) zeros(1, l_dOff_AO) ones(1, l_dOn_AI) ones(1, l_dOn_AO) ];
        attend = [ attend ones(1, l_dOff_AI) zeros(1, l_dOff_AO) ones(1, l_dOn_AI) zeros(1, l_dOn_AO) ];
    end
    
    [p,tbl] = anovan( data_vec, {direction, drug, attend}, 'model','full','varnames',{'direction','drug','attend'}, 'display','off', 'sstype', 1 );
    
    rslt.p = p;
    rslt.tbl = tbl;
end


% count_spikes counts the number of spikes in a given window (eg cue,
% blank, etc) for each trial and appends it to the trial_structure that was
% input. The window is given by two behavioral codes. eg:
% [121 126] - Event codes for attentional window [Cue On, Targets Blank];
function trials = count_spikes( trials, window, type )

    for i = 1:length(trials)
        
        % Find the times during the trial at which the window of interest (eg cue on, blank, etc) 
        % starts and ends.
        start_idx = trials(i).code_times( [trials(i).event_codes] == window(1) );
        end_idx   = trials(i).code_times( [trials(i).event_codes] == window(2) );
        
        % For cue period (attention period) - ditch the beginning of the
        % post-cue period and only take the last 300ms before the target
        % changes.
        if strcmp(type, 'attend')
            start_idx = (end_idx - 300);
        end
        if strcmp(type, 'visual')
            start_idx = (end_idx - 250);
        end
        
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
        if strcmp(type, 'attend')
            trials(i).attend_n_spikes = n_spikes;
        elseif strcmp(type, 'visual')
            trials(i).visual_n_spikes = n_spikes;
            
            % baseline
            dur = spike_end_idx - spike_start_idx;
            base_start_idx = spike_start_idx - dur;
            base_end_idx = spike_end_idx - dur;
            
            trials(i).baseline_n_spikes = sum(trials(i).spikes(base_start_idx:base_end_idx));
        end
    end
    

end


% segregate takes a structure of trials and makes a new index_struct which
% lists the indexes of trials matching particular conditions, eg:
% Direction 45º, Drug On, Attend Out
% Direction 180º, Drug Off, Attend In
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
