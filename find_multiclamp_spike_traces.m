function find_multiclamp_spike_traces()

    [data, digitalchannels, time] = readIgor_withDigital2();

    time_dur = time(end) - time(1);
    sample_rate = (length(time) / time_dur) * 1000; % In Hz
    high_pass_cutoff_freq = 300; %Hz
    
    threshold = 100; %mV
    
    
    % Rows will be time, columns will be traces
    data = squeeze( data(1,:,:) );
    
    % High pass filter the data at 300Hz
    filtered_data = highpass(data, high_pass_cutoff_freq, sample_rate);
    
    % Identify whether these traces do or do not contain a pulse
    sample_trace = data(:,1);
    
    thresh_idxs = (sample_trace >= threshold);
    
    diff(thresh_idxs)== 1;
    
    
    % If there is a pulse, chop it off
    
    
    % Identify a good threshold value automatically?
    
    % Search for spikes - voltages above or below a specific threshold
    
    % Plot the traces that putatively contain spikes
    
        
    figure();
    plot(trace1);
    
end