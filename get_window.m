% get_window identifies the event codes that delimit the specified window.
% A correct trial from the data_struct is passed in, as well as string
% identifying what window is desired:
% 'attend'      - cue onset to blank onset
% 'attContrast' - cue onset to blank onset
% 'wm'          - target offset to fixation offset % VERIFY
% 'visual'      - targets onset to either cue onset (attend trials) or target offset (wm trials)
% 'fixation'    - fixation onset to targets onset 
% 'fullNoMotor' - fixation onset to blank onset (attend trials) or to fixation offset (wm trials)
% 
% Will return a pair of numbers:
% end_code of the window
% win_length - how many samples should be counted back from that window
function win_info = get_window( correct_trial, window_string )
    attend_earlysession_flag = find(correct_trial.event_codes == 121);
    attend_latesession_flag  = find(correct_trial.event_codes == 133);
    wm_flag                  = find(correct_trial.event_codes == 153);

    if strcmp( window_string, 'attend' ) || strcmp( window_string, 'attContrast' )
        if attend_earlysession_flag, trial_window = [121 126]; % Need a function for this because I changed the event codes in Jan/Mar 2016 
        else, trial_window = [133 126];
        end
    elseif strcmp( window_string, 'wm' ) || strcmp( window_string, 'wm_last500' )% These are variable durations.
        trial_window = [155 161]; 
    elseif strcmp( window_string, 'visual' )
        if attend_earlysession_flag, trial_window = [124 121]; % Attend Trials Early Sessions
        elseif attend_latesession_flag, trial_window = [124 133]; % Attend Trials Late Sessions
        else, trial_window = [153 155]; % WM Trials
        end
    elseif strcmp( window_string, 'fixation' )
        if wm_flag, trial_window = [120 153]; % WM Trials
        else, trial_window = [120 124]; % Attend Trials
        end
    elseif strcmp( window_string, 'fullNoMotor' )
        if wm_flag, trial_window = [120 161];
        else, trial_window = [120 126]; 
        end
    end
    
    end_code = trial_window(2);
    win_length = get_win_length( correct_trial, trial_window, window_string );
    
    if strcmp( window_string, 'wm_last500' ) % Useful because WM ranges are variable duration.
        win_length = 500;
    end
    
    win_info = [end_code win_length];
end


function win_length = get_win_length( correct_trial, trial_window, window_string )

    %attend_earlysession_flag = find(correct_trial.event_codes == 121);
    %attend_latesession_flag  = find(correct_trial.event_codes == 133);
    wm_flag                  = find(correct_trial.event_codes == 153);

    fix_win_length = 300; % Always take the last 300ms of the fixation window
    wm_win_length  = 750; % Always take the last 750ms of the delay window in the WM task
    
    % If fixation window 
    if strcmp(window_string, 'fixation' )
        win_length = fix_win_length; 
        
    % fullNoMotor window
    elseif strcmp( window_string, 'fullNoMotor' )
        if wm_flag, fix_end_code = 153;
        else, fix_end_code = 124; 
            end_time = correct_trial.code_times( correct_trial.event_codes == trial_window(2) );
            fix_end_time = correct_trial.code_times( correct_trial.event_codes == fix_end_code );
            win_length = (end_time - fix_end_time) + fix_win_length;
        end
        
    % If visual or attend or working memory windows - consistent times
    else
        beg_time = correct_trial.code_times( correct_trial.event_codes == trial_window(1) );
        end_time = correct_trial.code_times( correct_trial.event_codes == trial_window(2) );
        win_length = end_time - beg_time;
    end
    
end