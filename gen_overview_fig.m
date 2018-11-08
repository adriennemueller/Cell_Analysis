% Will make and saveout a summary of att-in/att-out, drug-on/drug-off by
% direction for a given file.
function overview_fig = gen_overview_fig( data_struct, currents )

    all_paradigms = {'Attention', 'WM', 'Attention_Contrast' };
    file_paradigms = unique({data_struct.paradigm});
    usable_paradigms = file_paradigms(contains( file_paradigms, all_paradigms ) ); % TEST ME

    % Set up superfig
    overview_fig = figure('visible', 'off');
    base_overview_fig_numplots = 5; % 4 Drug off/on raster+sden + 1, overlapping d' (fix, vis, att/wm)
    numplots =  base_overview_fig_numplots * length(currents) * length(usable_paradigms);
    
    % Filter Data by current and paradigm
    % Loop through all paradigms in that file an
    for i = 2:length(currents) 
        for j = 1:length(usable_paradigms)

            % Filtered data_struct by paradigm
            data_struct = data_struct( strcmp( {data_struct.paradigm}, usable_paradigms{j} ) );
            
            if strcmp( usable_paradigms{j}, 'Attention_Contrast' )
                contrast_flag = 1; else, contrast_flag = 0;
            end
            
            retain_current = currents(1);
            eject_current = currents(i);
            
            window_str = 'fullNoMotor';
            corr_idx = find( [data_struct.trial_error] == 0 );
            
%             if isempty( corr_idx ) %%%HACK FOR TESTING - CAN BE REMOVED!
%                 disp()
%             end
%             
            sample_correct_trial = data_struct(corr_idx(1)); 
            %window = get_window( sample_correct_trial, window_str );
            
            % Filter by direction
            control_spikemat_in  = get_directional_spikemat( data_struct, retain_current, window_str, 'in', contrast_flag );
            control_spikemat_out = get_directional_spikemat( data_struct, retain_current, window_str, 'out', contrast_flag );
            drug_spikemat_in     = get_directional_spikemat( data_struct, eject_current, window_str, 'in', contrast_flag );
            drug_spikemat_out    = get_directional_spikemat( data_struct, eject_current, window_str, 'out', contrast_flag );

            % Plot Histograms and SDen Overlay for the subsets
            for k = 1:(length( control_spikemat_in ) /2 ) % Number of Directions / 2; because data the same in other four - just flipped in and out
                plot_data.ctrl_in  = logical(control_spikemat_in(k).spikes');
                plot_data.ctrl_out = logical(control_spikemat_out(k).spikes');
                plot_data.drug_in  = logical(drug_spikemat_in(k).spikes');
                plot_data.drug_out = logical(drug_spikemat_out(k).spikes');
             
                spike_sden_subplot = raster_sden_plot( plot_data, sample_correct_trial, window_str );

            end

            % Plot d' plot for each sub-window (Fix, Vis, Attend/WM)

            % Plot Anova -pval plot?
        end
    end
    
    % Make some sort of composite overview_fig for contrast data?
    
end

function output_plot = raster_sden_plot( plot_data, sample_correct_trial, window_str  )

    %output_plot = figure(); % Maybe?


    % Plot the four rasters with appropriate colors
    output_plot(1) = subplot(6,1,1);
    [ctrlOut_xs, ctrlOut_ys] = plotSpikeRaster( plot_data.ctrl_out, 'PlotType', 'scatter' );
    
    output_plot(2) = subplot(6,1,2);
    [ctrlIn_xs, ctrlIn_ys]   = plotSpikeRaster( plot_data.ctrl_in, 'PlotType', 'scatter' );
    
    output_plot(3) = subplot(6,1,3);
    [drugOut_xs, drugOut_ys] = plotSpikeRaster( plot_data.drug_out, 'PlotType', 'scatter' );

    output_plot(4) = subplot(6,1,4);
    [drugIn_xs, drugIn_ys]   = plotSpikeRaster( plot_data.drug_in, 'PlotType', 'scatter' );

    % If any of these plots is void of data, skip the rest of the code
    % here. HACK. FIX LATER.
    if isempty( plot_data.ctrl_out) || isempty( plot_data.ctrl_in ) || isempty( plot_data.drug_out ) || isempty( plot_data.drug_in )
        output_plot(5) = subplot(6,1,5); plot( 0,0);
        output_plot(6) = subplot(6,1,6); plot( 0,0);
        return
    end
    
    
    % Add the two Sden Overlay plots.

    % Calculate SDens, and SEs
    sden_kernelwidth = 50; %ms Gaussian
    [ctrl_out_sden, ctrl_out_se] = gen_sden_data( plot_data.ctrl_out, sden_kernelwidth );
    [ctrl_in_sden, ctrl_in_se]   = gen_sden_data( plot_data.ctrl_in, sden_kernelwidth );
    [drug_out_sden, drug_out_se] = gen_sden_data( plot_data.drug_out, sden_kernelwidth );
    [drug_in_sden, drug_in_se]   = gen_sden_data( plot_data.drug_in, sden_kernelwidth );
    
    % Control In vs Out
    output_plot(5) = subplot(6,1,5);
    x = [ 1 :(2 * length(ctrl_out_sden)) ]; % x, forwards and backwards
    yy = [ctrl_out_sden, fliplr(max(ctrl_out_sden, ctrl_in_sden))]; % Draw where attIn > attOut
    fill(x,yy,'k', 'FaceAlpha', 0.3, 'LineStyle', 'none');
    yy = [ctrl_out_sden, fliplr(min(ctrl_out_sden, ctrl_in_sden))]; % Draw where attIn < attOut
    fill(x,yy,'k', 'FaceAlpha', 0.1, 'LineStyle', 'none');
    
    % Drug In vs Out
    output_plot(6) = subplot(6,1,6);
    x = [ 1 :(2 * length(drug_out_sden)) ]; % x, forwards and backwards
    yy = [drug_out_sden, fliplr(max(drug_out_sden, drug_in_sden))]; % Draw where attIn > attOut
    fill(x,yy,'r', 'FaceAlpha', 0.3, 'LineStyle', 'none');
    yy = [drug_out_sden, fliplr(min(drug_out_sden, drug_in_sden))]; % Draw where attIn < attOut
    fill(x,yy,'r', 'FaceAlpha', 0.1, 'LineStyle', 'none');
    
    % Find the times for the different events during the trial
    event_struct = find_event_times( sample_correct_trial, window_str );
    for i = 1:4 % number of plots
        for j = 1:length(event_struct)
        
            curr_plot = output_plot(i);
            
            % Add the lines for the events 
            ylimits = ylim( curr_plot ); ylength = ylimits(2) - ylimits(1) + 1;
            line( output_plot(i), repmat(event_struct(j).e_time, ylength), ylimits(1):ylimits(2), 'Color','r' ); %%% TODO
    
            set(output_plot(i), 'xticklabel', {});
        end
    end

    % Add Event Strings to X Axis
    e_strings = {event_struct.e_string}; xlabel_times = [event_struct.e_time];
    set(output_plot(4), 'xtick', xlabel_times, 'xticklabel', e_strings);
    
    
    % Make plot pretty
end

% Generates a small struct with
% Currently convolves with gaussian kernel.
function [sden_avg sden_se] = gen_sden_data( spike_mat, kernel_width )

    % Generage sdens for all trials in spike_mat
    sden_mat = arrayfun( @(row_idx) spike_density( spike_mat( row_idx, :), kernel_width), (1:size(spike_mat,1) ), 'UniformOutput', 0 );
    sden_mat = cell2mat( sden_mat' );
    
    % Average them
    sden_avg = mean( sden_mat, 1 );
    
    % Generate standard error
    sden_std = std( sden_mat );
    num_trials = size( sden_mat, 1);
    sden_se = sden_std ./ sqrt( num_trials );
end


% This will only work for fullNoMotor trials
function event_struct = find_event_times( corr_trial, window_str )
    
    if ~strcmp( window_str, 'fullNoMotor' )
        disp( 'Inappropriate trial window selection. Only fullNoMotor case currently handled.' );
    end

    event_struct = struct;

    e_codes = corr_trial.event_codes;
    e_times = corr_trial.code_times;
  
    attend_earlysession_flag = ~isempty( find(e_codes == 121) );
    attend_latesession_flag  = ~isempty( find(e_codes == 133) );
    wm_flag                  = ~isempty( find(e_codes == 153) );
    
    offset = e_times(e_codes == 120) - 1; 
    
    % Fixation Onset the same for all trials.
    event_struct(1).e_string = 'Fix'; 
    event_struct(1).e_time = e_times(e_codes == 120) - offset; 
    
    % Visual Onset
    event_struct(2).e_string = 'Target'; 
    if attend_earlysession_flag || attend_latesession_flag
        event_struct(2).e_time = e_times(e_codes == 124) - offset; % Attend Conditions
    else
        event_struct(2).e_time = e_times(e_codes == 153) - offset; % WM Condition
    end
    
    if attend_earlysession_flag
        event_struct(3).e_string = 'Cue';
        event_struct(3).e_time   = e_times(e_codes == 121) - offset;
    elseif attend_latesession_flag
        event_struct(3).e_string = 'Cue';
        event_struct(3).e_time   = e_times(e_codes == 133) - offset;
    else % WM Condition
        event_struct(3).e_string = 'Delay';
        event_struct(3).e_time   = e_times(e_cdoes == 155) - offset;
    end
    
end
