function rslt = drug_vs_control_AttIn_AttOut( bhv_data, cell_data, plexon_tstamps )

% load bhv data
BHV = bhv_read(bhv_data);
cell_data = load(cell_data);
plexon_tstamps = load(plexon_tstamps);

assignin('base', 'BHV', BHV);
assignin('base', 'cell_data', cell_data);
assignin('base', 'plexon_tstamps', plexon_tstamps);

% find alignment

% The times of some events recorded by Monkeylogic, relative to the trial
% beginning.
intra_trial_event_times = BHV.CodeTimes;
% The times that Plexon records events, relative to the beginning of
% recording, converted to milliseconds.
plexon_event_times = plexon_tstamps.ans.Ts * 1000;

% Event numbers from Plexon (matching the timestamps in plexon_event_times)
% are in Strobed. Find triples of the start and end event codes (which
% appear to be recorded as 242 and 233, respectively).
% all_trial_idx has 2 columns, representing the start and end indexes of
% trials from Plexon.
all_trial_idx = candidate_trial_indexes( plexon_tstamps.ans.Strobed, 242, 233 );
% Timestamps of the start and end events from Plexon
all_trial_times = plexon_event_times(all_trial_idx);
% Duration of each trial from Plexon, in milliseconds
all_trial_lengths = all_trial_times(:,2) - all_trial_times(:,1);


% For every trial from BHV, attempt to find the trial from Plexon
% that most closely matches it in length. Store the index of the Plexon
% trial into idx(i), matching the BHV trial i.
idx = [];
for i = 1:length(intra_trial_event_times)
    times = intra_trial_event_times{i};
    trial_length = (times(end) - times(1));
    differ = abs(all_trial_lengths - trial_length);
    [~, j] = min(differ);
    idx = [idx ; j(1)];
end

% Assuming that the BHV trials represent a contiguous slice of trials from
% the Plexon data, find the index of the start of that run that would match
% each Plexon trial to the one of closest length in BHV.
% If the first BVH trial has been matched with the 180th Plexon trial, that
% would indicate a start index of 180. If the second trial was matched with
% Plexon trial 201, that would indicate a starting index of 200 (Plexon trial
% 200 would line up with BHV trial 1).
candidate_starting_indexes = idx - (0:length(idx)-1)';
% Find all suggested start indexes and their frequency of suggestion.
bins = unique(candidate_starting_indexes);
freq = histc(candidate_starting_indexes, bins);
% Take the one that was suggested the most times by trial matches. In
% practice, the highest fequency match is so by a significant margin.
[~, i] = max(freq);
probable_starting_index = bins(i);
assignin('base', 'probable_starting_index', probable_starting_index);

% chunk out cell data and add it to bhv_struct, at appropriate times

ntrials = length(intra_trial_event_times);
% The timestamps of the start and end trial events from Plexon, that we
% think correspond to the trials from BHV. This is an array of two columns,
% for the start and end timestamps, and the row indexes correspond to trial
% indexes from BHV.
trial_times = all_trial_times(probable_starting_index:probable_starting_index+ntrials-1, :);
assignin('base', 'trial_times', trial_times);

% The Plexon spike timestamps, in milliseconds.
spike_times = cell_data.adc041(:,1) * 1000;
assignin('base', 'spike_times', spike_times);

spikes = cell(ntrials, 1);
for i = 1:ntrials
    % Get only the spikes that happened between the start and end times for
    % each relevant trial.
    start_time = trial_times(i,1);
    end_time = trial_times(i,2);
    
    % Take some spikes either side of the trial so that we have some buffer
    % when convolving.
    
    window_start_time = start_time - 1000;
    window_end_time = end_time + 1000;
    
    trial_spike_times = spike_times(spike_times >= window_start_time & spike_times <= window_end_time);

    % Align spikes with behaviour events
    if ~isempty(trial_spike_times)
        % Event times from BHV, relative to trial start
        event_times = intra_trial_event_times{i};
        % Subtract the start time from Plexon and add the start
        % event time from BHV, so that spikes are
        % relative to the same time as the events from BHV.
        trial_spike_times = (trial_spike_times - start_time) + event_times(1);
        
        % Convert again so spike times are relative to the 5th event from
        % BHV (the second event after the start triple).
        trial_spike_times = trial_spike_times - event_times(5);

        spikes{i} = trial_spike_times;
    end
end

% Save the spike times in the BHV struct.
BHV.Spikes = spikes;
assignin('base', 'BHV2', BHV);

% plot sd for DrugOn and DrugOff for each direction

% break trials out into correct trials
correct_trial_idx = find(BHV.TrialError == 0);

% break trials out into different directions
NGroups = 8;
groupWidth = 2 * pi / NGroups;
Thetas = [BHV.UserVars.theta];
% figure();

averages_on = {};
averages_off = {};

    for group = 1:NGroups
        theta = group * groupWidth;

        % break trials into DrugOn and DrugOff
        drug_on_spikes = {};
        drug_off_spikes = {};
        drug_on_trials = 0;
        drug_off_trials = 0;
        for i = correct_trial_idx'
            trial_spike_times = BHV.Spikes{i};
            if close_angles(Thetas(i), theta, groupWidth)
                % Drug is on when analog output is high?
                drug = mean(BHV.AnalogData{i}.General.Gen3);
                drugIsOn = drug > 0;
                if drugIsOn
                    drug_on_spikes = [drug_on_spikes ; trial_spike_times];
                    drug_on_trials = drug_on_trials + 1;
                else
                    drug_off_spikes = [drug_off_spikes ; trial_spike_times];
                    drug_off_trials = drug_off_trials + 1;
                end
            end
        end

        assignin('base', 'drug_on_spikes', drug_on_spikes);
        assignin('base', 'drug_off_spikes', drug_off_spikes);

        % Convert to Matrix
        drug_on_full  = spikelist_fill( drug_on_spikes,  -1500, 4500 ); %WANT real trial end time instead of 3500
        drug_off_full = spikelist_fill( drug_off_spikes, -1500, 4500 ); %SAME SAME

        assignin('base', 'drug_on_full', drug_on_full);
        assignin('base', 'drug_off_full', drug_off_full);

        % Align on Target Onset (SEE ABOVE - USE DIFFERENT event_times(idx)?

        % Generate Mean Line + SE Lines
        drug_on_avgs  = gen_avgd_traces( drug_on_full );
        drug_off_avgs = gen_avgd_traces( drug_off_full );

        assignin('base', 'drug_on_avgs', drug_on_avgs);
        assignin('base', 'drug_off_avgs', drug_off_avgs);
        
        averages_on = [averages_on, drug_on_avgs];
        averages_off = [averages_off, drug_off_avgs];

        

        % Plot Nice Alphaed Overlays

%         subplot(3,3, get_subplotidx(group));
%         hold on;
%         title(['Group ' num2str(group) ' (Theta ' num2str(theta) ')']);
% 
%         if ~isempty(drug_on_spikes)
%             %[Y, X] = spike_hist(drug_on_spikes);
%             %plot(X, Y ./ drug_on_trials, 'red');
%             shadedErrorBar(1:length( drug_on_avgs), drug_on_avgs(1,:), drug_on_avgs(2,:), 'r', 1 );
% 
%         end
% 
%         if ~isempty(drug_off_spikes)
%             %[Y, X] = spike_hist(drug_off_spikes);
%             %plot(X, Y ./ drug_off_trials, 'black');
% 
%             shadedErrorBar(1:length( drug_off_avgs), drug_off_avgs(1,:), drug_off_avgs(2,:), 'k', 1 );
% 
%         end
%         
%         % HACK: Cut off at the window padding.
%         xlim([500 3500]);
%         
%         hold off;
    end
    
    assignin('base', 'averages_on', averages_on);
    assignin('base', 'averages_off', averages_off);
    
    averages_on_AttOut = rotate_cells( averages_on );
    averages_off_AttOut = rotate_cells( averages_off );
    
    
    
    figure();
    hold on;
    for i = 1:length(averages_off)
        % Get the averages out of the first row of each cell
        ON = averages_on{i}(1,:) - averages_on_AttOut{i}(1,:);
        OFF = averages_off{i}(1,:) - averages_off_AttOut{i}(1,:);
        
        plot(ON{i}(1:,) - OFF{i}(1,:));
    end
    xlim([500 3500]);
    line([1500 1500], [-100 100], 'YLimInclude', 'off', 'XLimInclude', 'off');
    
    ax = gca;
    ax.XTick = [1000 1500 2000 2500 3000 3500];
    ax.XTickLabel = {'-500','0','500','1000','1500','2000'};
    hold off;
    

%     subplot(3, 3, 5);
%     
%     [m, n] = size(cell_data.adc041);
%     hold on;
%     for i = 1:100:m
%         plot(cell_data.adc041(i,2:n), 'k');
%     end
%     hold off;
%     
%     xlabel( 'Time (AU)' );
%     ylabel( 'Voltage (uv)' );
    
end

function rslt = rotate_cells( cells )
    rslt = [cells(5:8), cells(1:4)];
end

function rslt = get_subplotidx( group );

    if group == 1
        rslt = 6;
    elseif group == 2
        rslt = 3;
    elseif group == 3
        rslt = 2;
    elseif group == 4
        rslt = 1;
    elseif group == 5
        rslt = 4;
    elseif group == 6
        rslt = 7;
    elseif group == 7
        rslt = 8;
    elseif group == 8
        rslt = 9;
    end

end

function rslt = gen_avgd_traces( spike_array )
    %max_idx = max(size(spike_array)); % NO?

    array_sizes = cellfun(@size, spike_array, 'UniformOutput', 0);
    max_idx = max(cell2mat(array_sizes));
    
    avg_array = {};
    % Make the Spike Arrays Spike Densities instead
  
    % Run through each index
    for i = 1:max_idx
        curr_list = [];
        for j = 1:length(spike_array)
            curr_sden = cell2mat(spike_array(j));
            % See if that index exists in a given spike array
            if i <= length(curr_sden) 
                % If it does, add its value to a new cell array, with that index
                curr_list = [curr_list curr_sden(i)];
            end
        end
        avg_array{i} = curr_list;

    end
    
    assignin('base', 'avg_array', avg_array);

    
    % Get the Avg and SE for the new array.
    avgs = cellfun(@mean,avg_array);
    stds = cellfun(@std, avg_array);
    ns   = cellfun(@length, avg_array);
   
    stes = stds ./ sqrt(ns);
    
    upper = avgs + stes;
    lower = avgs - stes;
    
    rslt = [avgs; stes];


end

% spike_times is a cell array of vectors of times in milliseconds, relative
% to trial starts, or some other aligned point.
function rslt_cell = spikelist_fill( spike_times, min_spike_time, max_spike_time )
    nTrials = length(spike_times);

    spike_times_rounded = cellfun(@round, spike_times, 'UniformOutput', 0);

    % Get the minimum and maximum relative time from across all cells, fi
    % we want to get a min and max spike time that are guaranteed to fit
    % our data, rather than get them from outside.
%     min_spike_time = min(cellfun(@min, spike_times_rounded));
%     max_spike_time = max(cellfun(@max, spike_times_rounded));

    rslt = zeros(nTrials, max_spike_time - min_spike_time);
    
    for i = 1:nTrials
        % Start with zeros and fill in a 1 for every millisecond containing
        % a spike. The first entry is at time min_spike_time, and the last
        % at time max_spike_time.
        spikes = zeros(1, max_spike_time - min_spike_time);
        spike_indexes = spike_times_rounded{i} - min_spike_time;
        % If we crash here, it's because min/max spike times don't
        % encompass the whole spike data. If we want to truncate the data,
        % we can do that here, or perhaps rearrange to do it after the
        % spike density function.
        spikes(spike_indexes) = 1;
        
        sden = spike_density( spikes, 200 ); % Gaussian
        %sden = spike_density_psp( spikes ); % PSP
        
        rslt(i,:) = sden;
    end
  
    
    % We want to return a cell array just because other code expects it. We
    % just split the array into one row per cell.
    % TODO: Fix the other code.
    rslt_cell = {};
    for i = 1:size(rslt,1)
        rslt_cell = [ rslt_cell, rslt(i,:) ];
    end
end

function [rslt, millis] = spike_hist( spike_times )
    A = floor(min(spike_times));
    B = floor(max(spike_times));
    millis = (A:B);
    rslt = histc( spike_times, millis );
    rslt = spike_density( rslt, 200 );
end

function rslt = close_angles(a, b, t)
    rslt = 0;
    d = abs(a - b);
    if d < t
        rslt = 1;
     elseif 2 * pi - d < t
         rslt = 1;
    end
end

% Finds triplets of start_event codes and end_event codes and returns the
% indices of the first of those triplets in an n by 2 matrix. 1 column for
% start event indices, 1 column for end event indices.
function rslt = candidate_trial_indexes( all_events, trial_start_event, trial_end_event )

    ntrials = floor(sum(all_events == trial_start_event) / 3);
    rslt = zeros(ntrials, 2);

    trial_idx = 1;
    event_idx = 1;
    while 1
        trial = find_trial( all_events, event_idx, trial_start_event, trial_end_event );
        if isempty(trial)
            break;
        end
        rslt(trial_idx,:) = trial;
        trial_idx = trial_idx + 1;
        event_idx = trial(2) + 1;
    end

    rslt = rslt(1:trial_idx-1,:);
end

function rslt = find_trial( all_events, start_index, trial_start_event, trial_end_event )

    rslt = [];

    if start_index + 2 > length(all_events)
        return;
    end

    starters = repmat(trial_start_event, 3, 1);
    enders = repmat(trial_end_event, 3, 1);
    start_idx = 0;
    for i = start_index:length(all_events)-2
        ev = all_events(i:i+2);
        if all(ev == starters)
            start_idx = i;
        elseif all(ev == enders)
            if start_idx > 0
                rslt = [start_idx, i+2];
                break;
            end
        end
    end

end
