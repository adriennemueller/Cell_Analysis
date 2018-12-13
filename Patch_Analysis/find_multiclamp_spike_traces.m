% threshold of 100 is good for monky3 and monky4, use 50 for the rest.

% Identify a good threshold value automatically?
function find_multiclamp_spike_traces(threshold)

    % Open Igor File
    current_direc = cd;
    [data, digitalchannels, time] = readIgor_withDigital2();
    cd(current_direc);
    
    time_dur = time(end) - time(1);
    sample_rate = (length(time) / time_dur) * 1000; % In Hz
    high_pass_cutoff_freq = 300; %Hz
            
    % Rows will be time, columns will be traces
    data = squeeze( data(1,:,:) );
 
    % Correction for Voltage Offset
    if data(1,1) < -300, data = data + 375; end 
    
    % High pass filter the data at 300Hz
    filtered_data = highpass(data, high_pass_cutoff_freq, sample_rate);
    
    % Identify whether these traces do or do not contain a pulse
    sample_trace = data(:,1);
    
    pulse_test_beg = 100;
    pulse_test_end = 200;

    test_vals = sample_trace(pulse_test_beg:pulse_test_end);
    pulse_test = sum(test_vals >= threshold);
    if pulse_test > 100, pulse = 1;
    else, pulse = 0;close all;close
    end
        
    % Search for spikes - voltages above or below a specific threshold
    if pulse
        end_idx = length(time) - 50;
        [potential_rows, potential_cols] = find( filtered_data(450:end_idx,:) >= threshold );
        filtered_data = filtered_data(475:end,:);
        time = time(475:end);
    else % If there is a pulse, do not search in the pulse range
        end_idx = length(time) - 200;
        [potential_rows, potential_cols] = find( filtered_data(200:end_idx,:) >= threshold );
        filtered_data = filtered_data(200:end_idx,:);
        time = time(200:end_idx);
        time = time - time(1); % Shift so starts at 0
    end
    
    potential_trace_cols = unique(potential_cols);
    potential_traces = filtered_data( :, potential_trace_cols );
    
    % Plot the traces that putatively contain spikes
    num_traces = size(potential_traces,2);
    for i = 1:num_traces
        figure();
        plot(time,potential_traces(:,i), 'LineWidth', 2);
        xlabel( 'Time (ms)', 'FontSize', 16, 'FontWeight', 'bold');
        ylabel( 'Voltage (\muV)', 'FontSize', 16, 'FontWeight', 'bold');
        set(gca,'FontSize',14)
        xlim( [time(1) time(end)] );
    end
    
%     for i = 1:size(filtered_data,2)
%         figure();
%         plot(filtered_data(i,:));
%     end
    
end