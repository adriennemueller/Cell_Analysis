
% Attend In - Attend Out

function rslt = attIn_attOut_Contrasts( full_trials, currents )

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
        window = get_attend_window( correct_trials );
        correct_trials = count_spikes( correct_trials, window, 'attend' );

        % For each direction - get idxs for attend in and attend out, when drug
        % is present and when drug is absent
        idx_struct = segregate( correct_trials );
        rslt = idx_struct; %%% TMP TMP TMP %%%

        
        event_durations = get_event_durations( correct_trials(1) );

        
        % Get D' for result of this
%         dmat = gen_dprime_struct( idx_struct, correct_trials );
%         
%         % Get Anova for this
         anova_mat = gen_anova_struct( idx_struct, correct_trials );
% 
%         % Get Attend In/Out Drug On/Off SDen averages and SEs
         sdens = get_attend_sdens( idx_struct, correct_trials, window(1) );
         sden_summs =  get_trial_sum( sdens );
% 
%         % Get Visual Reponse p-values %%% VERIFY VERIFY VERIFY %%%
%         vis_window = [124 window(1)]; %THIS WINDOW 1 THING IS CONFUSING. FIX.
%         correct_trials = count_spikes( correct_trials, vis_window, 'visual' );
%         vis_pval = gen_vis_pval( idx_struct, correct_trials );
%     
%         % Return results for each current separately
         rslt(i-1).current = eject_current;
%         rslt(i-1).dmat = dmat;
         rslt(i-1).anova_mat = anova_mat;
         rslt(i-1).sdens = sdens;
         rslt(i-1).sden_summs = sden_summs;
         rslt(i-1).event_durations = event_durations;

%         rslt(i-1).vis_pval = vis_pval;
% 
%         % MAKE A FUNCTION TO RUN THIS FUNCTION AND ADD THE OUTPUT TO THE
%         % SESSION_STRUCT AND SAVE IT OUT?
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
    if find(correct_trials(1).event_codes == 121);
        window = [121 126];
    else
        window = [133 126];
    end
end


function vis_pval = gen_vis_pval( idx_struct, trials )

    % Can combine AttendIn and AttendOut trials because identical during period
    % of initial stimulus presentation

    vis_pval = [];

    for i = 1:length(unique([idx_struct.theta])); % For the 8 directions
        
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
        
        
        %%% THIS IS AWFUL. FIX IT LATER 
        
        %%% C10 %%%
        rslt(i).C10dOaI_millis = trials(i).C10drugOff_attIn(1,:);
        C10dOaI_sdens  = trials(i).C10drugOff_attIn(2:end,:);
        [ts ms] = size(C10dOaI_sdens);
        rslt(i).C10dOaI_avg    = mean(C10dOaI_sdens, 1);
        rslt(i).C10dOaI_ste    = std(C10dOaI_sdens) / sqrt(ts); % Make sure std is in right dimension
        
        rslt(i).C10dOaO_millis = trials(i).C10drugOff_attOut(1,:);
        C10dOaO_sdens  = trials(i).C10drugOff_attOut(2:end,:);
        [ts ms] = size(C10dOaO_sdens);
        rslt(i).C10dOaO_avg    = mean(C10dOaO_sdens, 1);
        rslt(i).C10dOaO_ste    = std(C10dOaO_sdens) / sqrt(ts); % Make sure std is in right dimension

        rslt(i).C10dNaI_millis = trials(i).C10drugOn_attIn(1,:);
        C10dNaI_sdens  = trials(i).C10drugOn_attIn(2:end,:);
        [ts ms] = size(C10dNaI_sdens);
        rslt(i).C10dNaI_avg    = mean(C10dNaI_sdens, 1);
        rslt(i).C10dNaI_ste    = std(C10dNaI_sdens) / sqrt(ts); % Make sure std is in right dimension

        rslt(i).C10dNaO_millis = trials(i).C10drugOn_attOut(1,:);
        C10dNaO_sdens  = trials(i).C10drugOn_attOut(2:end,:);
        [ts ms] = size(C10dNaO_sdens);
        rslt(i).C10dNaO_avg    = mean(C10dNaO_sdens, 1);
        rslt(i).C10dNaO_ste    = std(C10dNaO_sdens) / sqrt(ts); % Make sure std is in right dimension

        %%% ATT In and ATT Out Combined %%%
        rslt(i).C10dN_millis = trials(i).C10drugOn_attIn(1,:);
        C10dN_sdens = vertcat(C10dNaI_sdens, C10dNaO_sdens);
        [ts ms] = size(C10dN_sdens);
        rslt(i).C10dN_avg    = mean(C10dN_sdens, 1);
        rslt(i).C10dN_ste    = std(C10dN_sdens) / sqrt(ts); % Make sure std is in right dimension

        rslt(i).C10dO_millis = trials(i).C10drugOn_attOut(1,:);
        C10dO_sdens  = vertcat(C10dOaI_sdens, C10dOaO_sdens);
        [ts ms] = size(C10dO_sdens);
        rslt(i).C10dO_avg    = mean(C10dO_sdens, 1);
        rslt(i).C10dO_ste    = std(C10dO_sdens) / sqrt(ts); % Make sure std is in right dimension
               
        
        
        %%% C15 %%%
        rslt(i).C15dOaI_millis = trials(i).C15drugOff_attIn(1,:);
        C15dOaI_sdens  = trials(i).C15drugOff_attIn(2:end,:);
        [ts ms] = size(C15dOaI_sdens);
        rslt(i).C15dOaI_avg    = mean(C15dOaI_sdens, 1);
        rslt(i).C15dOaI_ste    = std(C15dOaI_sdens) / sqrt(ts); % Make sure std is in right dimension
        
        rslt(i).C15dOaO_millis = trials(i).C15drugOff_attOut(1,:);
        C15dOaO_sdens  = trials(i).C15drugOff_attOut(2:end,:);
        [ts ms] = size(C15dOaO_sdens);
        rslt(i).C15dOaO_avg    = mean(C15dOaO_sdens, 1);
        rslt(i).C15dOaO_ste    = std(C15dOaO_sdens) / sqrt(ts); % Make sure std is in right dimension

        rslt(i).C15dNaI_millis = trials(i).C15drugOn_attIn(1,:);
        C15dNaI_sdens  = trials(i).C15drugOn_attIn(2:end,:);
        [ts ms] = size(C15dNaI_sdens);
        rslt(i).C15dNaI_avg    = mean(C15dNaI_sdens, 1);
        rslt(i).C15dNaI_ste    = std(C15dNaI_sdens) / sqrt(ts); % Make sure std is in right dimension

        rslt(i).C15dNaO_millis = trials(i).C15drugOn_attOut(1,:);
        C15dNaO_sdens  = trials(i).C15drugOn_attOut(2:end,:);
        [ts ms] = size(C15dNaO_sdens);
        rslt(i).C15dNaO_avg    = mean(C15dNaO_sdens, 1);
        rslt(i).C15dNaO_ste    = std(C15dNaO_sdens) / sqrt(ts); % Make sure std is in right dimension

        %%% ATT In and ATT Out Combined %%%
        rslt(i).C15dN_millis = trials(i).C15drugOn_attIn(1,:);
        C15dN_sdens = vertcat(C15dNaI_sdens, C15dNaO_sdens);
        [ts ms] = size(C15dN_sdens);
        rslt(i).C15dN_avg    = mean(C15dN_sdens, 1);
        rslt(i).C15dN_ste    = std(C15dN_sdens) / sqrt(ts); % Make sure std is in right dimension

        rslt(i).C15dO_millis = trials(i).C15drugOn_attOut(1,:);
        C15dO_sdens  = vertcat(C15dOaI_sdens, C15dOaO_sdens);
        [ts ms] = size(C15dO_sdens);
        rslt(i).C15dO_avg    = mean(C15dO_sdens, 1);
        rslt(i).C15dO_ste    = std(C15dO_sdens) / sqrt(ts); % Make sure std is in right dimension
               
    
        %%% 20 %%%
        rslt(i).C20dOaI_millis = trials(i).C20drugOff_attIn(1,:);
        C20dOaI_sdens  = trials(i).C20drugOff_attIn(2:end,:);
        [ts ms] = size(C20dOaI_sdens);
        rslt(i).C20dOaI_avg    = mean(C20dOaI_sdens, 1);
        rslt(i).C20dOaI_ste    = std(C20dOaI_sdens) / sqrt(ts); % Make sure std is in right dimension
        
        rslt(i).C20dOaO_millis = trials(i).C20drugOff_attOut(1,:);
        C20dOaO_sdens  = trials(i).C20drugOff_attOut(2:end,:);
        [ts ms] = size(C20dOaO_sdens);
        rslt(i).C20dOaO_avg    = mean(C20dOaO_sdens, 1);
        rslt(i).C20dOaO_ste    = std(C20dOaO_sdens) / sqrt(ts); % Make sure std is in right dimension

        rslt(i).C20dNaI_millis = trials(i).C20drugOn_attIn(1,:);
        C20dNaI_sdens  = trials(i).C20drugOn_attIn(2:end,:);
        [ts ms] = size(C20dNaI_sdens);
        rslt(i).C20dNaI_avg    = mean(C20dNaI_sdens, 1);
        rslt(i).C20dNaI_ste    = std(C20dNaI_sdens) / sqrt(ts); % Make sure std is in right dimension

        rslt(i).C20dNaO_millis = trials(i).C20drugOn_attOut(1,:);
        C20dNaO_sdens  = trials(i).C20drugOn_attOut(2:end,:);
        [ts ms] = size(C20dNaO_sdens);
        rslt(i).C20dNaO_avg    = mean(C20dNaO_sdens, 1);
        rslt(i).C20dNaO_ste    = std(C20dNaO_sdens) / sqrt(ts); % Make sure std is in right dimension

        %%% ATT In and ATT Out Combined %%%
        rslt(i).C20dN_millis = trials(i).C20drugOn_attIn(1,:);
        C20dN_sdens = vertcat(C20dNaI_sdens, C20dNaO_sdens);
        [ts ms] = size(C20dN_sdens);
        rslt(i).C20dN_avg    = mean(C20dN_sdens, 1);
        rslt(i).C20dN_ste    = std(C20dN_sdens) / sqrt(ts); % Make sure std is in right dimension

        rslt(i).C20dO_millis = trials(i).C20drugOn_attOut(1,:);
        C20dO_sdens  = vertcat(C20dOaI_sdens, C20dOaO_sdens);
        [ts ms] = size(C20dO_sdens);
        rslt(i).C20dO_avg    = mean(C20dO_sdens, 1);
        rslt(i).C20dO_ste    = std(C20dO_sdens) / sqrt(ts); % Make sure std is in right dimension        
        
        
        
        
        
        %%% C25 %%%
        rslt(i).C25dOaI_millis = trials(i).C25drugOff_attIn(1,:);
        C25dOaI_sdens  = trials(i).C25drugOff_attIn(2:end,:);
        [ts ms] = size(C25dOaI_sdens);
        rslt(i).C25dOaI_avg    = mean(C25dOaI_sdens, 1);
        rslt(i).C25dOaI_ste    = std(C25dOaI_sdens) / sqrt(ts); % Make sure std is in right dimension
        
        rslt(i).C25dOaO_millis = trials(i).C25drugOff_attOut(1,:);
        C25dOaO_sdens  = trials(i).C25drugOff_attOut(2:end,:);
        [ts ms] = size(C25dOaO_sdens);
        rslt(i).C25dOaO_avg    = mean(C25dOaO_sdens, 1);
        rslt(i).C25dOaO_ste    = std(C25dOaO_sdens) / sqrt(ts); % Make sure std is in right dimension

        rslt(i).C25dNaI_millis = trials(i).C25drugOn_attIn(1,:);
        C25dNaI_sdens  = trials(i).C25drugOn_attIn(2:end,:);
        [ts ms] = size(C25dNaI_sdens);
        rslt(i).C25dNaI_avg    = mean(C25dNaI_sdens, 1);
        rslt(i).C25dNaI_ste    = std(C25dNaI_sdens) / sqrt(ts); % Make sure std is in right dimension

        rslt(i).C25dNaO_millis = trials(i).C25drugOn_attOut(1,:);
        C25dNaO_sdens  = trials(i).C25drugOn_attOut(2:end,:);
        [ts ms] = size(C25dNaO_sdens);
        rslt(i).C25dNaO_avg    = mean(C25dNaO_sdens, 1);
        rslt(i).C25dNaO_ste    = std(C25dNaO_sdens) / sqrt(ts); % Make sure std is in right dimension

        %%% ATT In and ATT Out Combined %%%
        rslt(i).C25dN_millis = trials(i).C25drugOn_attIn(1,:);
        C25dN_sdens = vertcat(C25dNaI_sdens, C25dNaO_sdens);
        [ts ms] = size(C25dN_sdens);
        rslt(i).C25dN_avg    = mean(C25dN_sdens, 1);
        rslt(i).C25dN_ste    = std(C25dN_sdens) / sqrt(ts); % Make sure std is in right dimension

        rslt(i).C25dO_millis = trials(i).C25drugOn_attOut(1,:);
        C25dO_sdens  = vertcat(C25dOaI_sdens, C25dOaO_sdens);
        [ts ms] = size(C25dO_sdens);
        rslt(i).C25dO_avg    = mean(C25dO_sdens, 1);
        rslt(i).C25dO_ste    = std(C25dO_sdens) / sqrt(ts); % Make sure std is in right dimension
        
    end
end


function sdens = get_attend_sdens( idx_struct, trials, align_idx )

    drugOff_attIn = [];
    drugOff_attOut = [];
    drugOn_attIn = [];
    drugOn_attOut = [];

   for i = 1:length(unique([idx_struct.theta])); % For the 8 directions
        
        theta = map_direction(i);
        
        %%% THIS IS AWFUL. FIX IT LATER.
        
        %%% C10 %%%
        % Find the Drug Off, Attend In trials and
        row = find(([idx_struct.theta] == theta) & ([idx_struct.drug] == 0) & ([idx_struct.attend] == 1) & ([idx_struct.contrast] == 10) );
        sub_trials = trials(idx_struct(row).idxs);
        sdens(i).C10drugOff_attIn = align_trials( sub_trials, align_idx );
        
        % Same for Drug Off, Attend Out
        row = find(([idx_struct.theta] == theta) & ([idx_struct.drug] == 0) & ([idx_struct.attend] == 0) & ([idx_struct.contrast] == 10));
        sub_trials = trials(idx_struct(row).idxs);
        sdens(i).C10drugOff_attOut = align_trials( sub_trials, align_idx );
        
        % Find the Drug On, Attend In trials and
        row = find(([idx_struct.theta] == theta) & ([idx_struct.drug] == 1) & ([idx_struct.attend] == 1) & ([idx_struct.contrast] == 10));
        sub_trials = trials(idx_struct(row).idxs);
        sdens(i).C10drugOn_attIn = align_trials( sub_trials, align_idx );
        
        % Same for Drug On, Attend Out
        row = find(([idx_struct.theta] == theta) & ([idx_struct.drug] == 1) & ([idx_struct.attend] == 0) & ([idx_struct.contrast] == 10));
        sub_trials = trials(idx_struct(row).idxs);
        sdens(i).C10drugOn_attOut = align_trials( sub_trials, align_idx );

        
        %%% C15 %%%
        % Find the Drug Off, Attend In trials and
        row = find(([idx_struct.theta] == theta) & ([idx_struct.drug] == 0) & ([idx_struct.attend] == 1) & ([idx_struct.contrast] == 15));
        sub_trials = trials(idx_struct(row).idxs);
        sdens(i).C15drugOff_attIn = align_trials( sub_trials, align_idx );
        
        % Same for Drug Off, Attend Out
        row = find(([idx_struct.theta] == theta) & ([idx_struct.drug] == 0) & ([idx_struct.attend] == 0) & ([idx_struct.contrast] == 15));
        sub_trials = trials(idx_struct(row).idxs);
        sdens(i).C15drugOff_attOut = align_trials( sub_trials, align_idx );
        
        % Find the Drug On, Attend In trials and
        row = find(([idx_struct.theta] == theta) & ([idx_struct.drug] == 1) & ([idx_struct.attend] == 1) & ([idx_struct.contrast] == 15));
        sub_trials = trials(idx_struct(row).idxs);
        sdens(i).C15drugOn_attIn = align_trials( sub_trials, align_idx );
        
        % Same for Drug On, Attend Out
        row = find(([idx_struct.theta] == theta) & ([idx_struct.drug] == 1) & ([idx_struct.attend] == 0) & ([idx_struct.contrast] == 15));
        sub_trials = trials(idx_struct(row).idxs);
        sdens(i).C15drugOn_attOut = align_trials( sub_trials, align_idx );
        
        
        %%% C20 %%%
        % Find the Drug Off, Attend In trials and
        row = find(([idx_struct.theta] == theta) & ([idx_struct.drug] == 0) & ([idx_struct.attend] == 1) & ([idx_struct.contrast] == 20));
        sub_trials = trials(idx_struct(row).idxs);
        sdens(i).C20drugOff_attIn = align_trials( sub_trials, align_idx );
        
        % Same for Drug Off, Attend Out
        row = find(([idx_struct.theta] == theta) & ([idx_struct.drug] == 0) & ([idx_struct.attend] == 0) & ([idx_struct.contrast] == 20));
        sub_trials = trials(idx_struct(row).idxs);
        sdens(i).C20drugOff_attOut = align_trials( sub_trials, align_idx );
        
        % Find the Drug On, Attend In trials and
        row = find(([idx_struct.theta] == theta) & ([idx_struct.drug] == 1) & ([idx_struct.attend] == 1) & ([idx_struct.contrast] == 20));
        sub_trials = trials(idx_struct(row).idxs);
        sdens(i).C20drugOn_attIn = align_trials( sub_trials, align_idx );
        
        % Same for Drug On, Attend Out
        row = find(([idx_struct.theta] == theta) & ([idx_struct.drug] == 1) & ([idx_struct.attend] == 0) & ([idx_struct.contrast] == 20));
        sub_trials = trials(idx_struct(row).idxs);
        sdens(i).C20drugOn_attOut = align_trials( sub_trials, align_idx );
        
        
        
        %%% C25 %%%
        % Find the Drug Off, Attend In trials and
        row = find(([idx_struct.theta] == theta) & ([idx_struct.drug] == 0) & ([idx_struct.attend] == 1) & ([idx_struct.contrast] == 25));
        sub_trials = trials(idx_struct(row).idxs);
        sdens(i).C25drugOff_attIn = align_trials( sub_trials, align_idx );
        
        % Same for Drug Off, Attend Out
        row = find(([idx_struct.theta] == theta) & ([idx_struct.drug] == 0) & ([idx_struct.attend] == 0) & ([idx_struct.contrast] == 25));
        sub_trials = trials(idx_struct(row).idxs);
        sdens(i).C25drugOff_attOut = align_trials( sub_trials, align_idx );
        
        % Find the Drug On, Attend In trials and
        row = find(([idx_struct.theta] == theta) & ([idx_struct.drug] == 1) & ([idx_struct.attend] == 1) & ([idx_struct.contrast] == 25));
        sub_trials = trials(idx_struct(row).idxs);
        sdens(i).C25drugOn_attIn = align_trials( sub_trials, align_idx );
        
        % Same for Drug On, Attend Out
        row = find(([idx_struct.theta] == theta) & ([idx_struct.drug] == 1) & ([idx_struct.attend] == 0) & ([idx_struct.contrast] == 25));
        sub_trials = trials(idx_struct(row).idxs);
        sdens(i).C25drugOn_attOut = align_trials( sub_trials, align_idx );
        
        
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
    contrast = [];
    
    for i = 1:length(idx_struct)
        
        for j = 1:length(trials)
            
            if sum((j == idx_struct(i).idxs))
                
                theta_i = idx_struct(i).theta;
                if (theta_i >= 180), continue, end
                direction = [direction theta_i];
                
                drug_i = idx_struct(i).drug;
                drug = [drug drug_i];
                
                attend_i = idx_struct(i).attend;
                attend = [attend attend_i];
                
                contrast_i = idx_struct(i).contrast;
                contrast = [contrast contrast_i];
                
                data_vec_i = trials(j).attend_n_spikes;
                data_vec = [data_vec data_vec_i];
            else
                continue
            end
        end
    
    end
    
    [p,tbl] = anovan( data_vec, {direction, drug, attend, contrast}, 'model','full','varnames',{'direction','drug','attend', 'contrast'}, 'display','off', 'sstype', 1 );
     
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
    contrasts = [10 15 20 25];
    
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
         
            for contrast_i = 1:length(contrasts) 
                
            tmp_struct.contrast = contrasts(contrast_i);
            
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
                                            ([trials.theta] == map_direction(attend_direction)) & ...
                                             [trials.contrast] == tmp_struct.contrast );

                    if isempty(fieldnames(idx_struct))
                        idx_struct = tmp_struct;
                    else
                        idx_struct = [idx_struct tmp_struct];
                    end
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
