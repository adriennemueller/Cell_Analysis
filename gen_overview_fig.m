% Will make and saveout a summary of att-in/att-out, drug-on/drug-off by
% direction for a given file.
function overview_fig = gen_overview_fig( data_struct_in, currents )

    all_paradigms = {'Attention', 'WM', 'Attention_Contrast' };
    file_paradigms = unique({data_struct_in.paradigm});
    usable_paradigms = file_paradigms(contains( file_paradigms, all_paradigms ) );

    % Set up superfig
    overview_fig = figure();%'visible', 'off');

    % Filter Data by current and paradigm
    % Loop through all paradigms and currents in that file and make a separate
    % overview_fig for each one.
    
    direc_plot_num = 1; % Starting count for number of plots
    total_num_direc_plots = (length(currents) - 1) * length(usable_paradigms);
    
    % Establish figure size
    overview_fig_pos = get(overview_fig, 'Position'); %// gives x left, y bottom, width, height
    overview_fig_width = overview_fig_pos(3) * total_num_direc_plots;  overview_fig_height = overview_fig_pos(4);
    set(overview_fig, 'Position', [10 10 overview_fig_width overview_fig_height]);
    
    for i = 2:length(currents) 
        for j = 1:length(usable_paradigms)

            curr_paradigm = usable_paradigms{j};
            
            % Filtered data_struct by paradigm
            data_struct = data_struct_in( strcmp( {data_struct_in.paradigm}, curr_paradigm ) );
            
            retain_current = currents(1);
            eject_current  = currents(i);
            
            % Different Window for WM
                %%%% TMP %%%
    window_str = 'fullNoMotor';
    %%%%%
    
    
            corr_idx = find( [data_struct.trial_error] == 0 );
            sample_correct_trial = data_struct(corr_idx(1)); 
            
            % Filter by direction
            control_spikemat_in  = get_combined_spikemat( data_struct, retain_current, curr_paradigm, 'in' );
            control_spikemat_out = get_combined_spikemat( data_struct, retain_current, curr_paradigm, 'out' );
            drug_spikemat_in     = get_combined_spikemat( data_struct, eject_current, curr_paradigm, 'in' );
            drug_spikemat_out    = get_combined_spikemat( data_struct, eject_current, curr_paradigm, 'out' );

            % Plot Histograms and SDen Overlay for the subsets
            total_num_directions = length( control_spikemat_in ) / 2;
            direc_fig = figure();
            for k = 1:total_num_directions % Number of Directions / 2; because data the same in other four - just flipped in and out
                plot_data.ctrl_in  = logical(control_spikemat_in(k).spikes');
                plot_data.ctrl_out = logical(control_spikemat_out(k).spikes');
                plot_data.drug_in  = logical(drug_spikemat_in(k).spikes');
                plot_data.drug_out = logical(drug_spikemat_out(k).spikes');
                
                spike_sden_subplot = raster_sden_plot( plot_data, sample_correct_trial, window_str );
                
                direc_fig = insert_subpanel( direc_fig, spike_sden_subplot, k, total_num_directions );
            end
        
            overview_fig = append_direc_fig( overview_fig, direc_fig, curr_paradigm, currents(i), total_num_direc_plots, direc_plot_num );
            direc_plot_num = direc_plot_num + 1; % Suboptimal?
            
            % Plot d' plot for each sub-window (Fix, Vis, Attend/WM)

            % Plot Anova -pval plot?
        end
    end
    
    % Make some sort of composite overview_fig for contrast data?
    
end

function combined_spikemat = get_combined_spikemat( data_struct, current, paradigm, attend_type )

    if strcmp( paradigm, 'Attention_Contrast' )
        contrast_flag = 1; else, contrast_flag = 0;
    end

    blank_period = 100; % 100ms break in between sections
    
    % WM Trials
    if strcmp( paradigm, 'WM' )
        fix_spikemat  = get_directional_spikemat( data_struct, current, 'fixation', attend_type, contrast_flag );
        vis_spikemat  = get_directional_spikemat( data_struct, current, 'visual', attend_type, contrast_flag );
        main_spikemat = get_directional_spikemat( data_struct, current, 'wm_last500', attend_type, contrast_flag );
        rwd_spikemat  = get_directional_spikemat( data_struct, current, 'reward', attend_type, contrast_flag );
   
        num_trials = size( main_spikemat(1).spikes, 2 );
        blank_mat = zeros( blank_period, num_trials );
          
        combined_spikemat_spikes = horzcat( fix_spikemat(1).spikes, vis_spikemat(1).spikes, blank_mat, main_spikemat(1).spikes, blank_mat, rwd_spikemat(1).spikes );
        
        
    % ATTEND TRIALS
    else
        main_spikemat = get_directional_spikemat( data_struct, current, 'fullNoMotor', attend_type, contrast_flag );
        rwd_spikemat  = get_directional_spikemat( data_struct, current,  'reward', attend_type, contrast_flag );
        
        num_trials = size( main_spikemat(1).spikes, 2 );
        blank_mat = zeros( blank_period, num_trials );
        combined_spikemat_spikes = vertcat(  main_spikemat(1).spikes, blank_mat, rwd_spikemat(1).spikes );
    end
    
    combined_spikemat = main_spikemat;
    combined_spikemat(1).spikes = combined_spikemat_spikes;
    
end


function overview_fig = append_direc_fig( overview_fig, direc_fig, paradigm, current, total_num_direc_plots, direc_plot_num )

    figure( overview_fig ); % Make this figure the current figure
    
    % Get number of subfig panels in the current direc_fig
    sub_fig_N_x = 4; % 4 Positions
    sub_fig_N_y = 6; % 6 Subplots: Attend in/out for Control/Drug, + Control SDen + Drug SDen
    direc_fig_subplot_num = sub_fig_N_x * sub_fig_N_y;
    
    if total_num_direc_plots == 1
        overview_fig = direc_fig;
    else
        % Loop through all of the subfigs in this direc_fig
        hFigIAxes = findobj('Parent', direc_fig, 'Type','axes');
        hFigIAxes = remap_subplot_positions(hFigIAxes);
        row_width = sub_fig_N_x * total_num_direc_plots;
        for i = 1:direc_fig_subplot_num %24
            
            subfig_row = (ceil(i/sub_fig_N_x) - 1);
            subfig_col_in_panel = mod(i-1,sub_fig_N_x) + 1;
            curr_panel_num = subfig_row * row_width + subfig_col_in_panel;
            curr_panel_num = curr_panel_num + sub_fig_N_x * (direc_plot_num - 1);
            
            % Copy this subfig from the direc_fig into the appropriate position in the
            % overview_fig
            hAxes = hFigIAxes(i);
            subplot_child = get(hAxes, 'children');
            
            %subplot_child = get(direc_fig(i), 'children');
            curr_subplot = subplot( sub_fig_N_y, sub_fig_N_x * total_num_direc_plots, curr_panel_num );
            copyobj(subplot_child, curr_subplot);
        end
        subtitle_str = strcat(paradigm, " ", num2str(current), "nA");
        subtitle_spacing = 1 / (total_num_direc_plots+1);
        dim  = [subtitle_spacing * direc_plot_num, 0.96, 0.3, 0.01]; 
        annotation('textbox', dim, 'string', subtitle_str, 'FitBoxToText', 'on', 'LineStyle', 'none');
    end
    
end

function axes_list = remap_subplot_positions( axes_list )
    x = reshape(1:24,4,6);
    x = x';
    tmp_indices = flip(x(:),1);
    [~,new_indices] = sort(tmp_indices);
    axes_list = axes_list( new_indices );
end


function direc_fig = insert_subpanel( direc_fig, ss_subplot, direc_num, total_num_directions )
    
    num_subfig_panels = length(ss_subplot);
    figure(direc_fig); % Make this figure the current figure

    for i = 1:num_subfig_panels
        curr_panel_num = sub2ind( [total_num_directions, num_subfig_panels], direc_num, i );
        subplot_child = get(ss_subplot(i), 'children');
        curr_subplot = subplot( num_subfig_panels, total_num_directions, curr_panel_num );
        copyobj(subplot_child, curr_subplot);
    end
end


function output_plot = raster_sden_plot( plot_data, sample_correct_trial, window_str  )

    output_plot = figure();%'visible', 'off');

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
    
    
    %%% Add the two Sden Overlay plots %%%

    % Calculate SDens, and SEs
    sden_kernelwidth = 50; %ms Gaussian
    [ctrl_out_sden, ctrl_out_se] = gen_sden_data( plot_data.ctrl_out, sden_kernelwidth );
    [ctrl_in_sden, ctrl_in_se]   = gen_sden_data( plot_data.ctrl_in, sden_kernelwidth );
    [drug_out_sden, drug_out_se] = gen_sden_data( plot_data.drug_out, sden_kernelwidth );
    [drug_in_sden, drug_in_se]   = gen_sden_data( plot_data.drug_in, sden_kernelwidth );
    
    % Control In vs Out
    output_plot(5) = subplot(6,1,5); hold on;
    x = [ 1 : length(ctrl_out_sden), length(ctrl_out_sden) : -1 : 1];
    yy = [ctrl_out_sden, fliplr(max(ctrl_out_sden, ctrl_in_sden))]; % Draw where attIn > attOut
    fill(x,yy,'k', 'FaceAlpha', 0.3, 'LineStyle', 'none');
    yy = [ctrl_out_sden, fliplr(min(ctrl_out_sden, ctrl_in_sden))]; % Draw where attIn < attOut
    fill(x,yy,'k', 'FaceAlpha', 0.1, 'LineStyle', 'none'); hold off;
    
    % Drug In vs Out
    output_plot(6) = subplot(6,1,6); hold on;
    x = [ 1 : length(drug_out_sden), length(drug_out_sden) : -1 : 1];
    yy = [drug_out_sden, fliplr(max(drug_out_sden, drug_in_sden))]; % Draw where attIn > attOut
    fill(x,yy,'r', 'FaceAlpha', 0.3, 'LineStyle', 'none');
    yy = [drug_out_sden, fliplr(min(drug_out_sden, drug_in_sden))]; % Draw where attIn < attOut
    fill(x,yy,'r', 'FaceAlpha', 0.1, 'LineStyle', 'none'); hold off;
    
    % Find the times for the different events during the trial

    event_struct = find_event_times( sample_correct_trial, window_str );
    
    
    %%% Debug code
    for i = 1:length(event_struct)
        disp( strcat( event_struct(i).e_string,{': '}, num2str(event_struct(i).e_time) ) );
    end
    disp( 'Durations:' );
    disp( [event_struct(2:end).e_time] - [event_struct(1:end-1).e_time] );
    %%%
    
    
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
    
    event_struct = struct;

    e_codes = corr_trial.event_codes;
    e_times = corr_trial.code_times;
  
    attend_earlysession_flag = ~isempty( find(e_codes == 121) );
    attend_latesession_flag  = ~isempty( find(e_codes == 133) );
    wm_flag                  = ~isempty( find(e_codes == 153) );
    
    offset = e_times(e_codes == 120) - 1; 
    
    
    if strcmp( window_str, 'fullNoMotor' )
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
        end
    elseif strcmp( window_str, 'wm_last500' )    % WM Condition  %%% TEST TEST TEST TEST
            event_struct(1).e_string = 'Delay End';
            event_struct(1).e_time   = e_times(e_codes == 155) - offset;
        
    end
    
    
end
